# =============================================================================
# esmc — Equivalent Soil Mass
# Unified comparison function: esm_compare()
#
# Version: 0.3.0
#
# Changes from v0.2.0:
#   - smm -> mass in all internal and output column names
#   - Extrapolation now always performed (no more NA for TARGET_NOT_REACHED)
#   - Column `extrapolated` added to $layer, $cumulative, $diagnostics
#   - warning() issued when any profile requires extrapolation
# =============================================================================


#' Equivalent Soil Mass Comparison
#'
#' Standardizes soil organic carbon (SOC) stocks — or any other soil property
#' that varies with depth — to a common mineral soil mass reference, using
#' the Equivalent Soil Mass (ESM) framework.
#'
#' Soil samples collected at fixed depths contain different amounts of mineral
#' soil mass when bulk density (BD) differs among treatments or over time.
#' Comparing SOC stocks across unequal soil masses can produce artifacts
#' (false positives, biased estimates). \code{esm_compare} corrects for this
#' by evaluating a continuous C(M) curve — cumulative SOC stock as a function
#' of cumulative mineral soil mass — for each profile, and interpolating to a
#' common reference mass.
#'
#' When the reference mass exceeds the maximum observed mass of a profile,
#' the function extrapolates beyond the observed data range. Extrapolated
#' values are flagged with \code{extrapolated = TRUE} and may be biased;
#' use them with caution.
#'
#' @param data data.frame. Must contain columns: \code{upper}, \code{lower}
#'   (depth interval boundaries in cm), \code{bd} (bulk density in g/cm3),
#'   \code{trt} (treatment), the column named by \code{soil_prop}, and all
#'   columns listed in \code{profile_by} and \code{reference_by}.
#' @param soil_prop Character. Name of the soil property column to analyze
#'   (e.g. \code{"soc"} for soil organic carbon concentration).
#' @param profile_by Character vector. Columns that uniquely identify a soil
#'   profile/core (e.g. \code{c("site", "year", "trt", "rep")}).
#' @param output Character. What to return: \code{"both"} (layers and
#'   cumulative, default), \code{"layers"} (ESM stock per depth layer only),
#'   or \code{"cumulative"} (cumulative ESM stock at the reference mass only).
#'   Note: \code{esm_plot_spline()} requires \code{output = "both"} or
#'   \code{output = "cumulative"}.
#' @param reference_mode Character. How to define the reference mineral soil
#'   mass: \code{"min"} uses the minimum observed mass across profiles (avoids
#'   extrapolation), \code{"ref_trt"} uses the mean mass of a reference
#'   treatment, \code{"custom"} uses a user-supplied value.
#' @param reference_by Character vector or NULL. Columns defining strata
#'   within which the reference mass is computed independently (e.g.
#'   \code{c("site", "year")}). If NULL, a single global reference mass is
#'   used across all profiles.
#' @param ref_trt Value in column \code{trt} to use as the reference treatment.
#'   Required when \code{reference_mode = "ref_trt"}.
#' @param custom_mass Numeric. Reference mineral soil mass in kg/m2. Required
#'   when \code{reference_mode = "custom"}.
#' @param interp Character. Interpolation method for the C(M) curve:
#'   \code{"spline"} (cubic spline, recommended when >= 3 depth intervals are
#'   available), \code{"linear"} (recommended for 2 intervals), or
#'   \code{"step"} (constant within interval).
#' @param soil_prop_unit Character. Unit of \code{soil_prop}: \code{"pct"}
#'   (percent, 0-100) or \code{"g_kg"} (g per kg).
#' @param som_col Character or NULL. Name of a column containing measured
#'   soil organic matter (SOM) content, used to compute mineral soil mass.
#'   If NULL (default), mineral mass is approximated as total soil mass
#'   (i.e. SOM correction is skipped unless \code{som_from_soc = TRUE}).
#' @param som_unit Character. Unit of \code{som_col}: \code{"pct"},
#'   \code{"g_kg"}, or \code{"frac"} (0-1).
#' @param som_from_soc Logical. If TRUE, estimate SOM from SOC using
#'   \code{SOM = SOC * som_factor}.
#' @param soc_col Character or NULL. SOC column used to estimate SOM when
#'   \code{som_from_soc = TRUE} and SOC is in a different column than
#'   \code{soil_prop}. Defaults to \code{soil_prop}.
#' @param soc_unit Character. Unit of \code{soc_col} when it differs from
#'   \code{soil_prop}: \code{"pct"} or \code{"g_kg"}.
#' @param som_factor Numeric > 0. Conversion factor from SOC to SOM
#'   (default 2, i.e. SOM = 2 * SOC).
#' @param output_unit Character. Unit for stock outputs: \code{"kg_m2"}
#'   (default) or \code{"Mg_ha"} (adds convenience columns with suffix
#'   \code{_Mg_ha}; 1 kg/m2 = 10 Mg/ha).
#'
#' @return A list with the following elements (some depend on \code{output}):
#' \describe{
#'   \item{\code{interval}}{Per-interval computations: soil mass, mineral mass,
#'     cumulative mineral mass, and stock for each depth interval and profile.}
#'   \item{\code{profile_mass}}{C(M) curve nodes per profile: cumulative
#'     mineral mass and cumulative stock, including the (0, 0) anchor.
#'     Useful for diagnostic plotting with \code{esm_plot_spline}.}
#'   \item{\code{layer}}{(\code{output = "layers"} or \code{"both"})
#'     ESM-adjusted stock per depth layer and profile, computed as
#'     C(M_lower) - C(M_upper). Includes column \code{extrapolated}.}
#'   \item{\code{targets}}{(\code{output = "layers"} or \code{"both"})
#'     Reference mineral mass at each depth boundary.}
#'   \item{\code{cumulative}}{(\code{output = "cumulative"} or \code{"both"})
#'     Cumulative ESM stock per profile at the reference mass. Includes
#'     column \code{extrapolated}.}
#'   \item{\code{diagnostics}}{(\code{output = "cumulative"} or \code{"both"})
#'     Per-profile diagnostics: target mass, maximum observed mass, whether
#'     extrapolation was required, and any notes.}
#'   \item{\code{meta}}{List of settings used in the call.}
#' }
#'
#' @examples
#' \dontrun{
#' # Layer-wise ESM stocks (most common use case)
#' result <- esm_compare(
#'   data           = my_data,
#'   soil_prop      = "soc",
#'   profile_by     = c("site", "year", "trt", "rep"),
#'   output         = "layers",
#'   reference_mode = "min",
#'   interp         = "spline",
#'   soil_prop_unit = "pct",
#'   som_from_soc   = TRUE
#' )
#' result$layer
#'
#' # Cumulative ESM stock at reference treatment mass
#' result2 <- esm_compare(
#'   data           = my_data,
#'   soil_prop      = "soc",
#'   profile_by     = c("site", "year", "trt", "rep"),
#'   output         = "cumulative",
#'   reference_mode = "ref_trt",
#'   ref_trt        = "control",
#'   reference_by   = c("site", "year"),
#'   interp         = "spline",
#'   soil_prop_unit = "pct",
#'   som_from_soc   = TRUE
#' )
#' result2$cumulative
#' }
#'
#' @export
esm_compare <- function(
    data,
    soil_prop,
    profile_by,
    output         = c("both", "layers", "cumulative"),
    reference_mode = c("min", "ref_trt", "custom"),
    reference_by   = NULL,
    ref_trt        = NULL,
    custom_mass    = NULL,
    interp         = c("spline", "linear", "step"),
    soil_prop_unit = c("pct", "g_kg"),

    som_col        = NULL,
    som_unit       = c("pct", "g_kg", "frac"),
    som_from_soc   = FALSE,
    soc_col        = NULL,
    soc_unit       = c("pct", "g_kg"),
    som_factor     = 2,

    output_unit    = c("kg_m2", "Mg_ha")
) {

  # ---- Match arguments -------------------------------------------------------
  output         <- match.arg(output)
  reference_mode <- match.arg(reference_mode)
  interp         <- match.arg(interp)
  soil_prop_unit <- match.arg(soil_prop_unit)
  som_unit       <- match.arg(som_unit)
  soc_unit       <- match.arg(soc_unit)
  output_unit    <- match.arg(output_unit)

  # ---- Validate required columns ---------------------------------------------
  req_cols <- unique(c("upper", "lower", "bd", "trt",
                       soil_prop, profile_by, reference_by))
  missing  <- setdiff(req_cols, names(data))
  if (length(missing) > 0)
    stop("Missing required columns: ", paste(missing, collapse = ", "))

  if (reference_mode == "ref_trt" && is.null(ref_trt))
    stop("reference_mode = 'ref_trt' requires ref_trt.")

  if (reference_mode == "custom" && (is.null(custom_mass) || !is.finite(custom_mass)))
    stop("reference_mode = 'custom' requires custom_mass (kg/m2).")

  if (!is.null(som_col) && !som_col %in% names(data))
    stop("som_col not found in data: ", som_col)

  if (isTRUE(som_from_soc)) {
    if (is.null(soc_col)) soc_col <- soil_prop
    if (!soc_col %in% names(data))
      stop("soc_col not found in data: ", soc_col)
    if (!is.finite(som_factor) || som_factor <= 0)
      stop("som_factor must be a positive number.")
  } else {
    soc_col <- NULL
  }

  # ---- Drop rows with missing essential values --------------------------------
  df <- data[
    is.finite(data$upper)          &
    is.finite(data$lower)          &
    is.finite(data$bd)             &
    is.finite(data[[soil_prop]]),
    , drop = FALSE
  ]
  if (nrow(df) == 0)
    stop("No valid rows after removing NAs in upper / lower / bd / soil_prop.")

  # ---- Interval computations (full, for internal use) ------------------------
  soc_unit_eff <- if (isTRUE(som_from_soc))
                    (if (identical(soc_col, soil_prop)) soil_prop_unit else soc_unit)
                  else soc_unit

  interval_full <- .esm_compute_intervals(
    df             = df,
    soil_prop      = soil_prop,
    profile_by     = profile_by,
    soil_prop_unit = soil_prop_unit,
    som_col        = som_col,
    som_unit       = som_unit,
    som_from_soc   = som_from_soc,
    soc_col        = soc_col,
    soc_unit       = soc_unit_eff,
    som_factor     = if (isTRUE(som_from_soc)) som_factor else NULL
  )

  # Slim for user-facing $interval output
  interval_tbl <- .esm_slim_interval(interval_full)

  # ---- C(M) curve per profile (includes (0,0) anchor) -----------------------
  profile_mass <- .esm_profile_mass(interval_full, profile_by)

  # ---- Branch by output type ------------------------------------------------
  do_layers     <- output %in% c("layers", "both")
  do_cumulative <- output %in% c("cumulative", "both")

  layer_tbl   <- NULL
  targets_tbl <- NULL
  cum_out     <- NULL
  diag_out    <- NULL

  # --- Layers branch ----------------------------------------------------------
  if (do_layers) {
    layers_tbl  <- .esm_layers_from_data(interval_full, reference_by)

    targets_tbl <- .esm_targets_by_boundary(
      interval_tbl   = interval_full,
      profile_by     = profile_by,
      reference_by   = reference_by,
      reference_mode = reference_mode,
      ref_trt        = ref_trt,
      custom_mass    = custom_mass
    )

    layer_tbl <- .esm_layer_stocks_from_targets(
      profile_mass = profile_mass,
      layers_tbl   = layers_tbl,
      targets_tbl  = targets_tbl,
      profile_by   = profile_by,
      reference_by = reference_by,
      interp       = interp
    )
  }

  # --- Cumulative branch ------------------------------------------------------
  if (do_cumulative) {
    ref_tbl <- .esm_define_reference_mass(
      interval_tbl   = interval_full,
      reference_mode = reference_mode,
      reference_by   = reference_by,
      ref_trt        = ref_trt,
      custom_mass    = custom_mass,
      profile_by     = profile_by
    )

    fit_out  <- .esm_fit_and_evaluate(
      profile_mass = profile_mass,
      ref_tbl      = ref_tbl,
      profile_by   = profile_by,
      reference_by = reference_by,
      interp       = interp
    )

    cum_out  <- fit_out$cumulative
    diag_out <- fit_out$diagnostics
  }

  # ---- Warning if any extrapolation occurred ---------------------------------
  n_extrap_layers <- if (!is.null(layer_tbl))
    sum(layer_tbl$extrapolated, na.rm = TRUE) else 0L

  n_extrap_cum <- if (!is.null(cum_out))
    sum(cum_out$extrapolated, na.rm = TRUE) else 0L

  if (n_extrap_layers > 0 || n_extrap_cum > 0) {
    msg_parts <- character(0)
    if (n_extrap_layers > 0)
      msg_parts <- c(msg_parts,
        paste0(n_extrap_layers, " layer estimate(s) in $layer"))
    if (n_extrap_cum > 0)
      msg_parts <- c(msg_parts,
        paste0(n_extrap_cum, " profile(s) in $cumulative"))
    warning(
      "ESM extrapolation beyond observed mineral mass range detected in: ",
      paste(msg_parts, collapse = " and "), ".\n",
      "  These estimates required extrapolation because the reference mass ",
      "exceeded the maximum observed mineral mass of the profile.\n",
      "  Extrapolated values are flagged with extrapolated = TRUE and ",
      "may be biased. Use with caution.\n",
      "  Consider using reference_mode = 'min' to avoid extrapolation.",
      call. = FALSE
    )
  }

  # ---- Optional unit conversion ----------------------------------------------
  if (output_unit == "Mg_ha") {
    interval_tbl <- .esm_add_Mg_ha(interval_tbl)
    profile_mass <- .esm_add_Mg_ha(profile_mass)
    if (!is.null(targets_tbl)) targets_tbl <- .esm_add_Mg_ha(targets_tbl)
    if (!is.null(layer_tbl))   layer_tbl   <- .esm_add_Mg_ha(layer_tbl)
    if (!is.null(cum_out))     cum_out     <- .esm_add_Mg_ha(cum_out)
    if (!is.null(diag_out))    diag_out    <- .esm_add_Mg_ha(diag_out)
  }

  # ---- Assemble output list --------------------------------------------------
  out <- list(
    interval     = interval_tbl,
    profile_mass = profile_mass,
    meta = list(
      package        = "esmc",
      version        = "0.3.0",
      fn             = "esm_compare",
      principle      = "Depth is sampling; mineral soil mass is comparison.",
      output         = output,
      reference_mode = reference_mode,
      reference_by   = reference_by,
      ref_trt        = ref_trt,
      custom_mass    = custom_mass,
      interp         = interp,
      soil_prop      = soil_prop,
      soil_prop_unit = soil_prop_unit,
      output_unit    = output_unit,
      som = list(
        som_col      = som_col,
        som_unit     = som_unit,
        som_from_soc = som_from_soc,
        soc_col      = soc_col,
        soc_unit     = if (isTRUE(som_from_soc)) soc_unit_eff else NULL,
        som_factor   = if (isTRUE(som_from_soc)) som_factor else NULL
      )
    )
  )

  if (do_layers) {
    out$targets <- targets_tbl
    out$layer   <- layer_tbl
  }
  if (do_cumulative) {
    out$cumulative  <- cum_out
    out$diagnostics <- diag_out
  }

  out
}


# =============================================================================
# Internal helpers
# =============================================================================

# Convert kg/m2 columns to Mg/ha
.esm_add_Mg_ha <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  kg_cols  <- num_cols[grepl("_kg_m2$", num_cols)]
  for (cc in kg_cols) {
    new_name <- sub("_kg_m2$", "_Mg_ha", cc)
    if (!new_name %in% names(df)) df[[new_name]] <- df[[cc]] * 10
  }
  df
}

# Keep only user-facing columns in $interval
.esm_slim_interval <- function(interval_tbl) {
  keep <- c(
    "upper", "lower", "bd",
    "soil_mass_kg_m2",
    "mineral_mass_kg_m2", "cum_mineral_mass_kg_m2",
    "stock_kg_m2", "cum_stock_kg_m2"
  )
  keep <- intersect(keep, names(interval_tbl))

  id_cols <- setdiff(names(interval_tbl), c(
    "thick_cm", "valid_thick", "soil_prop_frac",
    "smm_g_cm2", "cum_smm_g_cm2", "soil_mass_g_cm2",
    "som_frac", "som_source", "som_factor",
    keep
  ))
  interval_tbl[, c(id_cols, keep), drop = FALSE]
}

# Build C(M) nodes per profile with (0,0) anchor
.esm_profile_mass <- function(interval_tbl, profile_by) {
  gid    <- interaction(interval_tbl[profile_by], drop = TRUE)
  groups <- split(interval_tbl, gid)

  out <- lapply(groups, function(g) {
    ok <- is.finite(g$cum_mineral_mass_kg_m2) & is.finite(g$cum_stock_kg_m2)
    g2 <- g[ok, c(profile_by, "cum_mineral_mass_kg_m2", "cum_stock_kg_m2"),
             drop = FALSE]
    g2 <- g2[order(g2$cum_mineral_mass_kg_m2), , drop = FALSE]

    if (nrow(g2) > 0 && g2$cum_mineral_mass_kg_m2[1] > 0) {
      z                          <- g2[1, , drop = FALSE]
      z$cum_mineral_mass_kg_m2   <- 0
      z$cum_stock_kg_m2          <- 0
      g2 <- rbind(z, g2)
    }
    g2
  })

  out <- out[!vapply(out, is.null, logical(1)) &
               vapply(out, nrow, integer(1)) > 0]
  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

# Per-interval computations: soil mass, mineral mass, stocks
.esm_compute_intervals <- function(
    df, soil_prop, profile_by, soil_prop_unit,
    som_col, som_unit, som_from_soc, soc_col, soc_unit, som_factor
) {
  gid    <- interaction(df[profile_by], drop = TRUE)
  groups <- split(df, gid)

  out <- lapply(groups, function(g) {
    g <- g[order(g$upper, g$lower), , drop = FALSE]

    g$thick_cm    <- g$lower - g$upper
    g$valid_thick <- is.finite(g$thick_cm) & g$thick_cm > 0

    # Total soil mass per interval (g/cm2)
    g$soil_mass_g_cm2 <- ifelse(g$valid_thick, g$bd * g$thick_cm, NA_real_)

    # SOM fraction (0-1)
    g$som_source <- "none"
    g$som_frac   <- 0
    g$som_factor <- NA_real_

    if (!is.null(som_col)) {
      g$som_source <- "measured"
      som_vals     <- g[[som_col]]
      g$som_frac   <- switch(som_unit,
        pct  = som_vals / 100,
        g_kg = som_vals / 1000,
        frac = som_vals
      )
    } else if (isTRUE(som_from_soc)) {
      g$som_source <- paste0("estimated_from_", soc_col)
      g$som_factor <- som_factor
      soc_vals     <- g[[soc_col]]
      soc_frac     <- switch(soc_unit,
        pct  = soc_vals / 100,
        g_kg = soc_vals / 1000
      )
      g$som_frac <- soc_frac * som_factor
    }

    som_clamped <- pmin(pmax(g$som_frac, 0), 1)

    # Mineral mass (g/cm2) and cumulative, then converted to kg/m2
    g$mineral_mass_g_cm2     <- g$soil_mass_g_cm2 * (1 - som_clamped)
    g$cum_mineral_mass_g_cm2 <- cumsum(
      replace(g$mineral_mass_g_cm2, is.na(g$mineral_mass_g_cm2), 0))

    g$soil_mass_kg_m2        <- g$soil_mass_g_cm2        * 10
    g$mineral_mass_kg_m2     <- g$mineral_mass_g_cm2     * 10
    g$cum_mineral_mass_kg_m2 <- g$cum_mineral_mass_g_cm2 * 10

    # Soil property as fraction -> interval stock (kg/m2)
    g$soil_prop_frac  <- if (soil_prop_unit == "pct") g[[soil_prop]] / 100
                         else                          g[[soil_prop]] / 1000
    g$stock_kg_m2     <- g$soil_prop_frac * g$mineral_mass_kg_m2
    g$cum_stock_kg_m2 <- cumsum(
      replace(g$stock_kg_m2, is.na(g$stock_kg_m2), 0))
    g
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# Reference mass per stratum — cumulative branch
.esm_define_reference_mass <- function(
    interval_tbl, reference_mode, reference_by, ref_trt, custom_mass, profile_by
) {
  by_cols  <- unique(c(profile_by, reference_by))
  prof_max <- stats::aggregate(
    interval_tbl$cum_mineral_mass_kg_m2,
    by  = interval_tbl[by_cols],
    FUN = function(x) max(x, na.rm = TRUE)
  )
  names(prof_max)[ncol(prof_max)] <- "max_cum_mass_kg_m2"

  if (is.null(reference_by)) {
    prof_max$.__stratum__ <- "ALL"
    strata_cols <- ".__stratum__"
  } else {
    strata_cols <- reference_by
  }

  if (reference_mode == "custom") {
    ref_tbl               <- unique(prof_max[strata_cols])
    ref_tbl$mass_ref_kg_m2 <- custom_mass
    return(ref_tbl)
  }

  if (reference_mode == "min") {
    ref_tbl <- stats::aggregate(
      prof_max$max_cum_mass_kg_m2,
      by  = prof_max[strata_cols],
      FUN = function(x) min(x, na.rm = TRUE)
    )
    names(ref_tbl)[ncol(ref_tbl)] <- "mass_ref_kg_m2"
    return(ref_tbl)
  }

  if (reference_mode == "ref_trt") {
    if (!"trt" %in% names(prof_max))
      stop("reference_mode = 'ref_trt' requires column 'trt' in profile_by / data.")
    prof_ref <- prof_max[prof_max$trt == ref_trt, , drop = FALSE]
    if (nrow(prof_ref) == 0)
      stop("No profiles found for ref_trt = ", ref_trt)
    ref_tbl <- stats::aggregate(
      prof_ref$max_cum_mass_kg_m2,
      by  = prof_ref[strata_cols],
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    names(ref_tbl)[ncol(ref_tbl)] <- "mass_ref_kg_m2"
    return(ref_tbl)
  }

  stop("Unknown reference_mode: ", reference_mode)
}

# Evaluate C(M_ref) per profile — cumulative branch
# Extrapolates when target exceeds observed range; flags with extrapolated = TRUE
.esm_fit_and_evaluate <- function(
    profile_mass, ref_tbl, profile_by, reference_by, interp
) {
  safe_rbind <- function(x) {
    x <- x[!vapply(x, is.null, logical(1))]
    if (length(x) == 0) return(data.frame())
    cols <- unique(unlist(lapply(x, names)))
    x2   <- lapply(x, function(d) {
      for (cc in cols) if (!cc %in% names(d)) d[[cc]] <- NA
      d[, cols, drop = FALSE]
    })
    do.call(rbind, x2)
  }

  if (is.null(reference_by)) {
    profile_mass$.__stratum__ <- "ALL"
    ref_tbl$.__stratum__      <- "ALL"
    strata_cols <- ".__stratum__"
  } else {
    strata_cols <- reference_by
  }

  gid    <- interaction(profile_mass[profile_by], drop = TRUE)
  groups <- split(profile_mass, gid)

  cum_list  <- list()
  diag_list <- list()

  for (nm in names(groups)) {
    pm  <- groups[[nm]]
    pm  <- pm[order(pm$cum_mineral_mass_kg_m2), , drop = FALSE]

    stratum_vals <- pm[1, strata_cols, drop = FALSE]
    idx          <- .esm_match_row(ref_tbl[strata_cols], stratum_vals)
    if (is.na(idx))
      stop("Could not match profile stratum to reference mass table.")
    mass_ref <- ref_tbl$mass_ref_kg_m2[idx]

    x     <- pm$cum_mineral_mass_kg_m2
    y     <- pm$cum_stock_kg_m2
    npts  <- length(x)
    max_x <- if (npts > 0) max(x, na.rm = TRUE) else NA_real_

    # Extrapolation flag: target exceeds observed range
    extrapolated <- is.finite(mass_ref) && is.finite(max_x) && mass_ref > max_x

    esm_cum <- NA_real_
    note    <- NULL

    if (interp == "step") {
      esm_cum <- .esm_eval_step(x, y, mass_ref)
    } else if (interp == "linear") {
      if (npts < 2) {
        note <- "INSUFFICIENT_POINTS_LINEAR"
      } else {
        # rule = 2: extrapolates using boundary value for linear;
        # for true extrapolation we use the slope of the last segment
        esm_cum <- .esm_eval_linear(x, y, mass_ref)
      }
    } else if (interp == "spline") {
      if (npts < 3) {
        note <- "INSUFFICIENT_POINTS_SPLINE"
        if (npts >= 2) esm_cum <- .esm_eval_linear(x, y, mass_ref)
      } else {
        f       <- stats::splinefun(x, y, method = "natural")
        esm_cum <- f(mass_ref)
      }
    }

    base_row <- pm[1, c(profile_by, strata_cols), drop = FALSE]

    cum_row                  <- base_row
    cum_row$target_mass_kg_m2 <- mass_ref
    cum_row$cum_stock_kg_m2  <- esm_cum
    cum_row$extrapolated     <- extrapolated
    cum_list[[nm]]           <- cum_row

    d                        <- base_row
    d$interp                 <- interp
    d$target_mass_kg_m2      <- mass_ref
    d$max_cum_mass_kg_m2     <- max_x
    d$extrapolated           <- extrapolated
    d$note                   <- if (!is.null(note)) note else NA_character_
    diag_list[[nm]]          <- d
  }

  list(
    cumulative  = safe_rbind(cum_list),
    diagnostics = safe_rbind(diag_list)
  )
}

# Step interpolation (and extrapolation: holds last value)
.esm_eval_step <- function(x, y, xout) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  if (xout <= x[1])          return(y[1])
  if (xout >= x[length(x)]) return(y[length(y)])
  y[max(which(x <= xout))]
}

# Linear interpolation with true extrapolation beyond range
.esm_eval_linear <- function(x, y, xout) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  n <- length(x)
  if (xout <= x[1]) {
    # Extrapolate left using slope of first segment
    slope <- if (n >= 2) (y[2] - y[1]) / (x[2] - x[1]) else 0
    return(y[1] + slope * (xout - x[1]))
  }
  if (xout >= x[n]) {
    # Extrapolate right using slope of last segment
    slope <- if (n >= 2) (y[n] - y[n-1]) / (x[n] - x[n-1]) else 0
    return(y[n] + slope * (xout - x[n]))
  }
  stats::approx(x, y, xout = xout, method = "linear", rule = 1)$y
}

# Row matching helper
.esm_match_row <- function(tbl, rowvals) {
  if (nrow(tbl) == 0) return(NA_integer_)
  ok <- rep(TRUE, nrow(tbl))
  for (cn in names(rowvals)) ok <- ok & (tbl[[cn]] == rowvals[[cn]])
  idx <- which(ok)
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

# Layer definitions from sampled depth intervals
.esm_layers_from_data <- function(interval_tbl, reference_by) {
  if (is.null(reference_by)) {
    out <- unique(interval_tbl[, c("upper", "lower"), drop = FALSE])
    out <- out[order(out$upper, out$lower), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  key    <- interaction(interval_tbl[reference_by], drop = TRUE)
  groups <- split(interval_tbl, key)
  out    <- lapply(groups, function(g)
    unique(g[, c(reference_by, "upper", "lower"), drop = FALSE]))
  out    <- do.call(rbind, out)
  out    <- out[do.call(order, out[, c(reference_by, "upper", "lower"),
                                   drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Reference mass at each depth boundary — layers branch
.esm_targets_by_boundary <- function(
    interval_tbl, profile_by, reference_by = NULL,
    reference_mode = c("min", "ref_trt", "custom"),
    ref_trt = NULL, custom_mass = NULL
) {
  reference_mode <- match.arg(reference_mode)
  bnds <- sort(unique(interval_tbl$lower))

  if (reference_mode == "custom") {
    # Case 1: single scalar -> same target mass for all boundaries
    if (is.numeric(custom_mass) && length(custom_mass) == 1 && is.null(names(custom_mass))) {
      out <- data.frame(lower = bnds, target_mass_kg_m2 = custom_mass,
                        stringsAsFactors = FALSE)
      return(out[order(out$lower), , drop = FALSE])
    }
    # Case 2: named vector -> names are lower boundaries, values are target masses
    if (is.numeric(custom_mass) && !is.null(names(custom_mass))) {
      out <- data.frame(lower = as.numeric(names(custom_mass)),
                        target_mass_kg_m2 = as.numeric(custom_mass),
                        stringsAsFactors = FALSE)
      return(out[order(out$lower), , drop = FALSE])
    }
    # Case 3: data.frame with columns lower + target_mass_kg_m2
    if (is.data.frame(custom_mass)) {
      needed <- c("lower", "target_mass_kg_m2")
      miss   <- setdiff(needed, names(custom_mass))
      if (length(miss) > 0)
        stop("custom_mass data.frame is missing columns: ",
             paste(miss, collapse = ", "))
      if (!is.null(reference_by)) {
        miss2 <- setdiff(reference_by, names(custom_mass))
        if (length(miss2) > 0)
          stop("custom_mass data.frame is missing reference_by columns: ",
               paste(miss2, collapse = ", "))
      }
      out       <- custom_mass
      out$lower <- as.numeric(out$lower)
      out$target_mass_kg_m2 <- as.numeric(out$target_mass_kg_m2)
      return(out[do.call(order, out[, c(reference_by, "lower"),
                                    drop = FALSE]), , drop = FALSE])
    }
    stop(
      "custom_mass must be one of:\n",
      "  - a single numeric value (same target mass for all layers)\n",
      "  - a named numeric vector (names = lower boundaries, values = target masses)\n",
      "  - a data.frame with columns 'lower' and 'target_mass_kg_m2'"
    )
  }

  build_global <- function(fun_bnd) {
    do.call(rbind, lapply(bnds, function(b) {
      sub <- interval_tbl[interval_tbl$lower == b, , drop = FALSE]
      data.frame(lower = b, target_mass_kg_m2 = fun_bnd(sub),
                 stringsAsFactors = FALSE)
    }))
  }

  build_by_group <- function(gcols, fun_bnd) {
    key    <- interaction(interval_tbl[gcols], drop = TRUE)
    groups <- split(interval_tbl, key)
    do.call(rbind, lapply(groups, function(g) {
      base     <- unique(g[, gcols, drop = FALSE])[1, , drop = FALSE]
      bnds_grp <- sort(unique(g$lower))
      tt   <- do.call(rbind, lapply(bnds_grp, function(b) {
        sub <- g[g$lower == b, , drop = FALSE]
        data.frame(lower = b, target_mass_kg_m2 = fun_bnd(sub),
                   stringsAsFactors = FALSE)
      }))
      cbind(base, tt)
    }))
  }

  fun_bnd <- if (reference_mode == "min") {
    function(sub) min(sub$cum_mineral_mass_kg_m2, na.rm = TRUE)
  } else {
    if (is.null(ref_trt))
      stop("reference_mode = 'ref_trt' requires ref_trt.")
    function(sub) {
      sub2 <- sub[sub$trt == ref_trt, , drop = FALSE]
      mean(sub2$cum_mineral_mass_kg_m2, na.rm = TRUE)
    }
  }

  out <- if (is.null(reference_by)) build_global(fun_bnd) else
    build_by_group(reference_by, fun_bnd)

  out <- out[is.finite(out$target_mass_kg_m2) & out$target_mass_kg_m2 >= 0, ,
             drop = FALSE]
  rownames(out) <- NULL
  out
}

# Evaluate C(M) at each boundary and compute layer stocks — layers branch
# Extrapolates when target exceeds observed range; flags with extrapolated = TRUE
.esm_layer_stocks_from_targets <- function(
    profile_mass, layers_tbl, targets_tbl,
    profile_by, reference_by = NULL,
    interp = c("spline", "linear", "step")
) {
  interp <- match.arg(interp)

  # Add M = 0 boundary
  if (is.null(reference_by)) {
    if (!any(targets_tbl$lower == 0))
      targets_tbl <- rbind(
        data.frame(lower = 0, target_mass_kg_m2 = 0), targets_tbl)
  } else {
    strata  <- unique(targets_tbl[, reference_by, drop = FALSE])
    add0    <- cbind(strata, data.frame(lower = 0, target_mass_kg_m2 = 0))
    targets_tbl <- rbind(add0, targets_tbl)
  }

  # Evaluate C(M) at a single mass value; returns list(value, extrapolated)
  eval_C <- function(pm, M) {
    M  <- as.numeric(M)
    if (!is.finite(M)) return(list(value = NA_real_, extrapolated = FALSE))

    pm <- pm[order(pm$cum_mineral_mass_kg_m2), , drop = FALSE]
    x  <- pm$cum_mineral_mass_kg_m2
    y  <- pm$cum_stock_kg_m2
    if (length(x) == 0) return(list(value = NA_real_, extrapolated = FALSE))

    max_x        <- x[length(x)]
    extrapolated <- M > max_x

    value <- if (interp == "step") {
      .esm_eval_step(x, y, M)
    } else if (interp == "linear") {
      .esm_eval_linear(x, y, M)
    } else {
      # spline
      if (length(unique(x)) < 3) {
        .esm_eval_linear(x, y, M)
      } else {
        stats::splinefun(x, y, method = "natural")(M)
      }
    }
    list(value = value, extrapolated = extrapolated)
  }

  gid       <- interaction(profile_mass[profile_by], drop = TRUE)
  pm_groups <- split(profile_mass, gid)

  out_list <- vector("list", length(pm_groups))
  names(out_list) <- names(pm_groups)

  for (i in seq_along(pm_groups)) {
    pm      <- pm_groups[[i]]
    prof_id <- pm[1, profile_by, drop = FALSE]

    if (is.null(reference_by)) {
      tsub <- targets_tbl
      lyr  <- layers_tbl
    } else {
      key  <- prof_id[, reference_by, drop = FALSE]
      tsub <- targets_tbl
      lyr  <- layers_tbl
      for (cc in reference_by) {
        tsub <- tsub[tsub[[cc]] == key[[cc]], , drop = FALSE]
        lyr  <- lyr[lyr[[cc]]  == key[[cc]], , drop = FALSE]
      }
    }

    m_map <- setNames(tsub$target_mass_kg_m2, as.character(tsub$lower))

    res <- lapply(seq_len(nrow(lyr)), function(j) {
      u   <- lyr$upper[j]
      l   <- lyr$lower[j]
      M_u <- if (as.character(u) %in% names(m_map)) m_map[[as.character(u)]] else NA_real_
      M_l <- if (as.character(l) %in% names(m_map)) m_map[[as.character(l)]] else NA_real_

      res_u <- eval_C(pm, M_u)
      res_l <- eval_C(pm, M_l)

      data.frame(
        prof_id,
        upper                    = u,
        lower                    = l,
        target_mass_upper_kg_m2  = M_u,
        target_mass_lower_kg_m2  = M_l,
        cum_upper_kg_m2          = res_u$value,
        cum_lower_kg_m2          = res_l$value,
        layer_stock_kg_m2        = res_l$value - res_u$value,
        extrapolated             = res_u$extrapolated | res_l$extrapolated,
        stringsAsFactors         = FALSE
      )
    })
    out_list[[i]] <- do.call(rbind, res)
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
