# =============================================================================
# esmc — Equivalent Soil Mass
# Unified comparison function: esm_compare()
#
# Version: 0.4.0
#
# Changes from v0.3.2:
#   - Default interpolation is now interp = "hyman" (cubic spline with Hyman
#     monotonicity filter), selected after a simulation benchmark
#   - interp = "hyman" uses the secant of the last interval as the slope at
#     the deepest node and extrapolates linearly with it
#   - Profiles are checked for contiguous layers from 0 cm (.esm_check_profiles)
#
# Changes from v0.3.1:
#   - New interpolation options: interp = "monotone" (Fritsch-Carlson monotone
#     cubic Hermite spline) and interp = "hyman" (cubic spline with Hyman
#     monotonicity filter, as in von Haden et al. 2020)
#   - All curve evaluations use a single helper, .esm_interp_fun()
#
# Changes from v0.3.0:
#   - Layers branch: profiles are evaluated on a non-overlapping partition of
#     the union of sampled boundaries within each stratum (no duplicated
#     boundaries when depth increments differ among profiles)
#   - Warning when reference masses are not monotone with depth
#   - Stocks computed as concentration x total soil mass (previously x mineral
#     mass); mineral mass is used only as the ESM coordinate
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
#' Stocks are computed as concentration times total soil mass, because
#' concentrations are expressed per unit of total soil mass. Mineral soil
#' mass (total mass minus SOM, when SOM is supplied or estimated) is used only
#' as the coordinate along which profiles are standardized.
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
#'   \code{c("site", "year")}). Must be a subset of \code{profile_by}. If
#'   NULL, a single global reference mass is used across all profiles.
#' @param ref_trt Value in column \code{trt} to use as the reference treatment.
#'   Required when \code{reference_mode = "ref_trt"}.
#' @param custom_mass Reference mineral soil mass(es) in kg/m2, required when
#'   \code{reference_mode = "custom"}. One of:
#'   \itemize{
#'     \item a single number: the cumulative reference mass. Defines no
#'       layers, so \code{output = "layers"} gives an error and
#'       \code{output = "both"} returns only the cumulative stock (with a
#'       warning).
#'     \item a named numeric vector, names = lower boundaries in cm (e.g.
#'       \code{c("30" = 380, "50" = 640, "100" = 1290)}): one target mass per
#'       boundary, the same for every stratum.
#'     \item a data.frame with columns \code{lower} and
#'       \code{target_mass_kg_m2}, plus the \code{reference_by} columns for
#'       targets that differ among strata.
#'   }
#'   With per-boundary targets, layers are defined by the supplied boundaries
#'   (0 to the first boundary, then between consecutive boundaries) and the
#'   cumulative stock is evaluated at the deepest supplied boundary.
#' @param interp Character. Interpolation method for the cumulative curve
#'   C(M) (cumulative stock vs cumulative mineral mass). Spline methods require
#'   at least 2 depth intervals (3 points including the (0, 0) anchor):
#'   \itemize{
#'     \item \code{"hyman"} (default): cubic spline with Hyman monotonicity
#'       filter (\code{stats::splinefun(method = "hyman")}), as in von Haden et
#'       al. (2020), except that the slope at the deepest node is the secant
#'       of the last interval (the concentration per unit mineral mass of the
#'       deepest sampled layer). Never decreasing, so layer stocks cannot be
#'       negative. Extrapolates linearly with that slope. With the end
#'       condition of \code{stats::splinefun} the last interval and the
#'       extrapolation were unstable with few depth intervals.
#'     \item \code{"spline"}: natural cubic spline (\code{method = "natural"}).
#'       Default up to v0.3.2. Extrapolates linearly.
#'     \item \code{"monotone"}: monotone cubic Hermite spline (Fritsch-Carlson;
#'       \code{method = "monoH.FC"}). Performs well with many depth intervals.
#'       Extrapolates linearly with the slope of the last segment, i.e. the
#'       concentration of the deepest sampled layer.
#'     \item \code{"linear"}: linear interpolation (recommended when only one
#'       depth interval is available). Extrapolates linearly.
#'     \item \code{"step"}: constant within interval.
#'   }
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
#'     C(M_lower) - C(M_upper). Layers are a non-overlapping partition of all
#'     sampled boundaries within each stratum (defined by \code{reference_by}),
#'     or the boundaries supplied in \code{custom_mass}, so every profile is
#'     evaluated on the same layers and layer stocks add up to the cumulative
#'     stock at the deepest boundary. This equals \code{$cumulative} when the
#'     cumulative reference mass is the target of the deepest boundary (always
#'     with per-boundary \code{custom_mass}; with \code{"min"} or
#'     \code{"ref_trt"}, when all profiles reach the deepest boundary). When depth increments differ among profiles (e.g.
#'     between years), layers not sampled in a given profile are estimated by
#'     interpolation along that profile's C(M) curve. Includes column
#'     \code{extrapolated}.}
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
#'   interp         = "hyman",
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
#'   interp         = "hyman",
#'   soil_prop_unit = "pct",
#'   som_from_soc   = TRUE
#' )
#' result2$cumulative
#' }
#'
#' @importFrom stats setNames
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
    interp         = c("hyman", "spline", "monotone", "linear", "step"),
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

  if (!is.null(reference_by)) {
    rb_out <- setdiff(reference_by, profile_by)
    if (length(rb_out) > 0)
      stop("reference_by must be a subset of profile_by (every profile must ",
           "belong to a single stratum). Not in profile_by: ",
           paste(rb_out, collapse = ", "), call. = FALSE)
  }

  custom_tbl <- NULL
  if (reference_mode == "custom") {
    custom_tbl <- .esm_custom_as_df(custom_mass, reference_by)
    if (is.null(custom_tbl)) {
      # single value: one cumulative reference mass, no per-boundary targets
      if (output == "layers")
        stop("A single custom_mass defines one cumulative reference mass; ",
             "layer stocks need one target mass per depth boundary. Supply ",
             "a named vector or a data.frame (see ?esm_compare), or use ",
             "output = 'cumulative'.", call. = FALSE)
      if (output == "both") {
        warning("A single custom_mass defines only the cumulative reference ",
                "mass: layer stocks are not computed (output = 'cumulative').",
                call. = FALSE)
        output <- "cumulative"
      }
    }
  }

  if (reference_mode == "ref_trt") {
    if (!ref_trt %in% data$trt)
      stop("ref_trt = '", ref_trt, "' not found in column 'trt'.", call. = FALSE)
    pid   <- interaction(data[profile_by], drop = TRUE, sep = " | ")
    n_trt <- tapply(as.character(data$trt), pid, function(x) length(unique(x)))
    if (any(n_trt > 1))
      stop("Column 'trt' varies within profile(s): ",
           paste(names(n_trt)[n_trt > 1], collapse = "; "),
           ". profile_by must identify single profiles.", call. = FALSE)
    if (!is.null(reference_by)) {
      skey    <- interaction(data[reference_by], drop = TRUE, sep = " | ")
      has_ref <- tapply(data$trt == ref_trt, skey, any)
      if (any(!has_ref))
        stop("ref_trt = '", ref_trt, "' is missing in stratum/strata: ",
             paste(names(has_ref)[!has_ref], collapse = "; "), call. = FALSE)
    }
  }

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

  # ---- Check profile geometry (contiguous layers from 0, no overlaps) --------
  df <- .esm_check_profiles(df, profile_by)

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
    layers_tbl  <- if (is.null(custom_tbl))
      .esm_layers_from_data(interval_full, reference_by) else
      .esm_layers_from_custom(custom_tbl, interval_full, reference_by)

    targets_tbl <- .esm_targets_by_boundary(
      interval_tbl   = interval_full,
      profile_by     = profile_by,
      reference_by   = reference_by,
      reference_mode = reference_mode,
      ref_trt        = ref_trt,
      custom_mass    = custom_mass
    )

    .esm_check_targets_monotone(targets_tbl, reference_by)

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
      version        = "0.4.0",
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
    # Stock = concentration (per unit total soil mass) x total soil mass.
    # Mineral mass is used only as the coordinate for mass standardization.
    g$stock_kg_m2     <- g$soil_prop_frac * g$soil_mass_kg_m2
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
  by_cols  <- unique(c(profile_by, reference_by,
                       if (reference_mode == "ref_trt") "trt"))
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
    ct <- .esm_custom_as_df(custom_mass, reference_by)
    if (is.null(ct)) {
      ref_tbl                <- unique(prof_max[strata_cols])
      ref_tbl$mass_ref_kg_m2 <- custom_mass
      return(ref_tbl)
    }
    # per-boundary targets: cumulative reference = deepest supplied boundary
    if (is.null(reference_by)) {
      deep <- ct[which.max(ct$lower), , drop = FALSE]
      return(data.frame(.__stratum__ = "ALL",
                        mass_ref_kg_m2 = deep$target_mass_kg_m2))
    }
    key  <- interaction(ct[reference_by], drop = TRUE)
    deep <- do.call(rbind, lapply(split(ct, key), function(g)
      g[which.max(g$lower), c(reference_by, "target_mass_kg_m2"), drop = FALSE]))
    names(deep)[names(deep) == "target_mass_kg_m2"] <- "mass_ref_kg_m2"
    rownames(deep) <- NULL
    return(deep)
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

    if (interp == "linear" && npts < 2) {
      note <- "INSUFFICIENT_POINTS_LINEAR"
    } else {
      f <- .esm_interp_fun(x, y, interp,
                           label = paste(pm[1, profile_by], collapse = " | "))
      esm_cum <- f(mass_ref)
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

  tidy <- function(x) {
    x <- safe_rbind(x)
    x$.__stratum__ <- NULL      # internal column (reference_by = NULL)
    rownames(x) <- NULL
    x
  }
  list(
    cumulative  = tidy(cum_list),
    diagnostics = tidy(diag_list)
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

# Interpolating function for the C(M) curve, shared by esm_compare() and
# esm_plot_spline(). Returns a function of mass (vectorized).
#   spline   : natural cubic spline, linear extrapolation
#   monotone : Fritsch-Carlson monotone cubic Hermite spline; linear
#              extrapolation with the end slopes (last-segment secant)
#   hyman    : cubic spline with Hyman monotonicity filter (von Haden et al.
#              2020) inside the observed range; linear extrapolation with the
#              end slopes of the spline (stats::spline would use a cubic)
#   linear   : piecewise linear, linear extrapolation with end-segment slopes
#   step     : constant within interval
.esm_interp_fun <- function(x, y, interp, label = NULL) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  if (interp == "step")
    return(function(xout) vapply(xout, function(m) .esm_eval_step(x, y, m), numeric(1)))
  if (interp == "linear")
    return(function(xout) vapply(xout, function(m) .esm_eval_linear(x, y, m), numeric(1)))
  if (length(unique(x)) < 3) {
    stop(
      "interp = '", interp, "' requires at least 2 depth intervals (3 C(M) ",
      "points including the (0, 0) anchor), but ",
      if (!is.null(label)) paste0("profile '", label, "'") else "a profile",
      " has only ", length(unique(x)) - 1, " interval(s). ",
      "Use interp = 'linear' instead."
    )
  }
  method <- switch(interp,
                   spline   = "natural",
                   monotone = "monoH.FC",
                   hyman    = "hyman",
                   stop("Unknown interp: ", interp))
  f <- stats::splinefun(x, y, method = method)
  if (interp != "hyman") return(f)

  # Hyman: node slopes from the Hyman-filtered spline, except at the last node,
  # where the slope is the secant of the last interval (the C/mass ratio of the
  # deepest sampled layer). The "fmm" end condition of stats::splinefun is
  # unreliable with few nodes and biased both the last interval and the linear
  # extrapolation. The neighbouring node is re-filtered so the Hermite cubic
  # stays monotone. Outside the observed range, splinefunH extrapolates
  # linearly with the end slopes.
  n   <- length(x)
  m   <- f(x, deriv = 1)
  sec <- diff(y) / diff(x)
  m[n] <- sec[n - 1]
  if (n >= 3)
    m[n - 1] <- max(0, min(m[n - 1], 3 * min(sec[n - 2], sec[n - 1])))
  stats::splinefunH(x, y, m)
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

# Layer definitions: non-overlapping partition of the union of sampled
# boundaries within each stratum. When depth increments differ among
# profiles (e.g. between years), every profile is evaluated on the same
# set of contiguous, non-overlapping layers, so each boundary appears once
# per profile and layer stocks add up to cumulative stocks.
.esm_layers_from_data <- function(interval_tbl, reference_by) {
  make_partition <- function(g) {
    b <- sort(unique(c(g$upper, g$lower)))
    data.frame(upper = utils::head(b, -1), lower = utils::tail(b, -1))
  }
  if (is.null(reference_by)) {
    out <- make_partition(interval_tbl)
    rownames(out) <- NULL
    return(out)
  }
  key    <- interaction(interval_tbl[reference_by], drop = TRUE)
  groups <- split(interval_tbl, key)
  out    <- do.call(rbind, lapply(groups, function(g) {
    part <- make_partition(g)
    cbind(g[rep(1, nrow(part)), reference_by, drop = FALSE], part)
  }))
  out <- out[do.call(order, out[, c(reference_by, "upper", "lower"),
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
    ct <- .esm_custom_as_df(custom_mass, reference_by)
    if (is.null(ct))
      stop("Layer targets need custom_mass with one value per boundary.")
    if (!is.null(reference_by) && !all(reference_by %in% names(ct))) {
      # same targets for every stratum
      strata <- unique(interval_tbl[, reference_by, drop = FALSE])
      ct <- merge(strata, ct, by = NULL)
    }
    ct <- ct[, c(reference_by, "lower", "target_mass_kg_m2"), drop = FALSE]
    ct <- ct[do.call(order, ct[, c(reference_by, "lower"), drop = FALSE]), ,
             drop = FALSE]
    rownames(ct) <- NULL
    return(ct)
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
      rownames(base) <- NULL
      cbind(base, tt, row.names = NULL)
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
    interp = c("hyman", "spline", "monotone", "linear", "step")
) {
  interp <- match.arg(interp)

  # Add M = 0 boundary
  if (is.null(reference_by)) {
    if (!any(targets_tbl$lower == 0))
      targets_tbl <- rbind(
        data.frame(lower = 0, target_mass_kg_m2 = 0), targets_tbl)
  } else {
    strata  <- unique(targets_tbl[, reference_by, drop = FALSE])
    rownames(strata) <- NULL
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

    value <- .esm_interp_fun(x, y, interp)(M)
    list(value = value, extrapolated = extrapolated)
  }

  gid       <- interaction(profile_mass[profile_by], drop = TRUE)
  pm_groups <- split(profile_mass, gid)

  out_list <- vector("list", length(pm_groups))
  names(out_list) <- names(pm_groups)

  for (i in seq_along(pm_groups)) {
    pm      <- pm_groups[[i]]
    prof_id <- pm[1, profile_by, drop = FALSE]
    rownames(prof_id) <- NULL

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


# Warn if reference masses do not increase with depth within a stratum
# (possible when boundaries come from different sets of profiles)
.esm_check_targets_monotone <- function(targets_tbl, reference_by) {
  key <- if (is.null(reference_by)) rep("ALL", nrow(targets_tbl)) else
    interaction(targets_tbl[reference_by], drop = TRUE)
  bad <- vapply(split(targets_tbl, key), function(g) {
    g <- g[order(g$lower), , drop = FALSE]
    any(diff(g$target_mass_kg_m2) <= 0)
  }, logical(1))
  if (any(bad))
    warning("Reference masses do not increase monotonically with depth in: ",
            paste(names(bad)[bad], collapse = ", "),
            ". Some layer stocks may be negative. This can happen when depth ",
            "increments differ among profiles; consider a finer reference_by.",
            call. = FALSE)
  invisible(NULL)
}


# Check that every profile is a contiguous sequence of layers starting at 0.
# The cumulative mass coordinate is only valid under that condition: a gap
# (e.g. a layer dropped because of a missing value) shifts every deeper node
# and, through reference_mode = "min", the targets of the whole stratum.
#   - overlapping or duplicated layers -> error
#   - profile not starting at 0        -> profile dropped, warning
#   - gap                              -> profile truncated at the gap, warning
.esm_check_profiles <- function(df, profile_by) {
  gid    <- interaction(df[profile_by], drop = TRUE, sep = " | ")
  groups <- split(df, gid)
  bad_overlap <- character(0); not_zero <- character(0); gapped <- character(0)
  out <- lapply(names(groups), function(nm) {
    g <- groups[[nm]]
    g <- g[order(g$upper, g$lower), , drop = FALSE]
    if (any(g$lower <= g$upper)) {
      bad_overlap <<- c(bad_overlap, nm); return(g)
    }
    if (nrow(g) > 1 && any(g$upper[-1] < g$lower[-nrow(g)])) {
      bad_overlap <<- c(bad_overlap, nm); return(g)
    }
    if (g$upper[1] != 0) { not_zero <<- c(not_zero, nm); return(NULL) }
    if (nrow(g) > 1) {
      gap <- which(g$upper[-1] > g$lower[-nrow(g)])
      if (length(gap) > 0) {
        gapped <<- c(gapped, nm)
        g <- g[seq_len(gap[1]), , drop = FALSE]
      }
    }
    g
  })
  if (length(bad_overlap) > 0)
    stop("Overlapping, duplicated or non-positive-thickness layers in profile(s): ",
         paste(bad_overlap, collapse = "; "),
         ". Check that profile_by uniquely identifies each profile.", call. = FALSE)
  if (length(not_zero) > 0)
    warning(length(not_zero), " profile(s) do not start at 0 cm (after removing ",
            "rows with missing values) and were dropped: ",
            paste(not_zero, collapse = "; "), call. = FALSE)
  if (length(gapped) > 0)
    warning(length(gapped), " profile(s) have a gap between layers (possibly ",
            "rows removed because of missing values) and were truncated at the ",
            "first gap: ", paste(gapped, collapse = "; "), call. = FALSE)
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0)
    stop("No valid profiles left after checking layer geometry.", call. = FALSE)
  rownames(out) <- NULL
  out
}


# Parse custom_mass. Returns NULL for a single value (cumulative reference
# only) or a data.frame with columns lower, target_mass_kg_m2 and, optionally,
# the reference_by columns.
.esm_custom_as_df <- function(custom_mass, reference_by) {
  err <- function(...) stop(..., call. = FALSE)
  if (is.null(custom_mass))
    err("reference_mode = 'custom' requires custom_mass (kg/m2).")
  if (is.data.frame(custom_mass)) {
    miss <- setdiff(c("lower", "target_mass_kg_m2"), names(custom_mass))
    if (length(miss) > 0)
      err("custom_mass data.frame is missing columns: ", paste(miss, collapse = ", "))
    extra_rb <- intersect(reference_by, names(custom_mass))
    if (length(extra_rb) > 0 && !setequal(extra_rb, reference_by))
      err("custom_mass data.frame must contain all reference_by columns or none.")
    if (length(extra_rb) == 0 &&
        any(!names(custom_mass) %in% c("lower", "target_mass_kg_m2")) &&
        is.null(reference_by))
      err("custom_mass data.frame has stratum columns but reference_by is NULL.")
    out <- custom_mass
    out$lower             <- as.numeric(out$lower)
    out$target_mass_kg_m2 <- as.numeric(out$target_mass_kg_m2)
  } else if (is.numeric(custom_mass) && length(custom_mass) == 1 &&
             is.null(names(custom_mass))) {
    if (!is.finite(custom_mass) || custom_mass <= 0)
      err("custom_mass must be a positive, finite number (kg/m2).")
    return(NULL)
  } else if (is.numeric(custom_mass) && !is.null(names(custom_mass))) {
    out <- data.frame(lower = suppressWarnings(as.numeric(names(custom_mass))),
                      target_mass_kg_m2 = as.numeric(custom_mass))
    if (anyNA(out$lower))
      err("Names of custom_mass must be lower boundaries in cm (e.g. c('30' = 380)).")
  } else {
    err("custom_mass must be one of:\n",
        "  - a single number (cumulative reference mass only)\n",
        "  - a named numeric vector (names = lower boundaries in cm)\n",
        "  - a data.frame with columns 'lower' and 'target_mass_kg_m2' ",
        "(plus reference_by columns, optionally)")
  }
  if (any(!is.finite(out$lower)) || any(!is.finite(out$target_mass_kg_m2)) ||
      any(out$lower <= 0) || any(out$target_mass_kg_m2 <= 0))
    err("custom_mass boundaries and target masses must be positive and finite.")
  key <- if (all(reference_by %in% names(out)) && !is.null(reference_by))
    interaction(out[reference_by], drop = TRUE) else rep("ALL", nrow(out))
  if (any(vapply(split(out$lower, key), anyDuplicated, integer(1)) > 0))
    err("custom_mass has duplicated lower boundaries.")
  out
}

# Layers defined by the boundaries supplied in custom_mass (0 to first
# boundary, then between consecutive boundaries), per stratum when needed
.esm_layers_from_custom <- function(ct, interval_tbl, reference_by) {
  part <- function(b) { b <- sort(unique(c(0, b)))
    data.frame(upper = utils::head(b, -1), lower = utils::tail(b, -1)) }
  if (is.null(reference_by)) return(part(ct$lower))
  strata <- unique(interval_tbl[, reference_by, drop = FALSE])
  out <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i) {
    st <- strata[i, , drop = FALSE]
    b  <- if (all(reference_by %in% names(ct))) {
      sel <- rep(TRUE, nrow(ct))
      for (cc in reference_by) sel <- sel & ct[[cc]] == st[[cc]]
      ct$lower[sel]
    } else ct$lower
    if (length(b) == 0)
      stop("custom_mass has no targets for stratum: ",
           paste(unlist(st), collapse = " | "), call. = FALSE)
    cbind(st[rep(1, length(unique(b))), , drop = FALSE], part(b))
  }))
  rownames(out) <- NULL
  out
}
