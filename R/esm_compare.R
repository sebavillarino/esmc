# =============================================================================
# esmc — Equivalent Soil Mass
# Core comparison: esm_compare()
#
# Returns ESM-adjusted cumulative stock at a common reference mineral mass (smm_ref).
# Uses a continuous C(M) curve per core (cumulative stock vs cumulative mineral mass),
# evaluated via step | linear | spline.
#
# Key design choices:
# - reference_by can be NULL (global reference mass across all strata)
# - outputs are "clean": interval output is slimmed (no *_g_cm2, no internal flags)
# - uses generic names: cum_stock_* (not cumC)
#
# NOTE:
# - interval_full is used for ALL internal computations.
# - interval_tbl is slimmed ONLY for returning to the user.
# =============================================================================


#' ESM comparison (Equivalent Soil Mass Comparison)
#'
#' Depth is sampling; mineral soil mass is comparison.
#'
#' @param data data.frame with required columns: upper, lower, bd, trt, plus rvar and core_by columns.
#' @param rvar string. Concentration column name (e.g., "soc").
#' @param core_by character vector defining a core/profile (e.g., c("Site","yr","trt","rep")).
#' @param reference_mode "ref_trt" | "min" | "custom"
#' @param reference_by character vector defining strata for reference mass, or NULL for global.
#' @param ref_trt value in `trt` used as reference treatment (required if reference_mode="ref_trt").
#' @param custom_smm numeric. Reference mineral soil mass (kg/m2) if reference_mode="custom".
#' @param interp "step" | "linear" | "spline"
#' @param rvar_unit "pct" | "g_kg"
#' @param som_col optional measured SOM column
#' @param som_unit "pct" | "g_kg" | "frac"
#' @param som_from_soc if TRUE, estimate SOM from SOC using som_factor
#' @param soc_col SOC column for SOM estimation (defaults to rvar)
#' @param soc_unit "pct" | "g_kg"
#' @param som_factor numeric > 0
#' @param extrapolation "none" | "to_depth" (reserved; v1 behaves as none)
#' @param min_core_length_cm numeric
#' @param out_unit "kg_m2" or "Mg_ha" (adds *_Mg_ha columns)
#'
#' @return list with interval, profile_mass, cumulative, diagnostics, meta
#' @export
esm_compare <- function(
    data,
    rvar,
    core_by,
    reference_mode = c("ref_trt", "min", "custom"),
    reference_by   = NULL,
    ref_trt        = NULL,
    custom_smm     = NULL,
    interp         = c("step", "linear", "spline"),
    rvar_unit      = c("pct", "g_kg"),

    som_col        = NULL,
    som_unit       = c("pct", "g_kg", "frac"),
    som_from_soc   = FALSE,
    soc_col        = NULL,
    soc_unit       = c("pct", "g_kg"),
    som_factor     = 2,

    extrapolation  = c("none", "to_depth"),
    min_core_length_cm = 0,
    out_unit       = c("kg_m2", "Mg_ha")
) {
  reference_mode <- match.arg(reference_mode)
  interp         <- match.arg(interp)
  rvar_unit      <- match.arg(rvar_unit)
  som_unit       <- match.arg(som_unit)
  soc_unit       <- match.arg(soc_unit)
  extrapolation  <- match.arg(extrapolation)
  out_unit       <- match.arg(out_unit)

  # ---- Validate required columns ----
  req_cols <- unique(c("upper", "lower", "bd", "trt", rvar, core_by, reference_by))
  missing <- setdiff(req_cols, names(data))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

  if (reference_mode == "ref_trt" && is.null(ref_trt)) {
    stop("reference_mode='ref_trt' requires ref_trt (a value present in column `trt`).")
  }
  if (reference_mode == "custom" && (is.null(custom_smm) || !is.finite(custom_smm))) {
    stop("reference_mode='custom' requires custom_smm (kg/m2).")
  }

  # ---- SOM config ----
  if (!is.null(som_col) && !som_col %in% names(data)) stop("som_col not found in data: ", som_col)

  if (isTRUE(som_from_soc)) {
    if (is.null(soc_col)) soc_col <- rvar
    if (!soc_col %in% names(data)) stop("soc_col not found in data: ", soc_col)
    if (!is.finite(som_factor) || som_factor <= 0) stop("som_factor must be > 0.")
  } else {
    soc_col <- NULL
  }

  # ---- Drop essential NAs (strict v1) ----
  df <- data[
    is.finite(data$upper) &
      is.finite(data$lower) &
      is.finite(data$bd) &
      is.finite(data[[rvar]]),
    ,
    drop = FALSE
  ]
  if (nrow(df) == 0) stop("No rows left after removing NAs in upper/lower/bd/rvar.")

  # ---- Interval computations (FULL for internal use) ----
  interval_full <- .esm_compute_intervals(
    df = df,
    rvar = rvar,
    core_by = core_by,
    rvar_unit = rvar_unit,
    min_core_length_cm = min_core_length_cm,
    som_col = som_col,
    som_unit = som_unit,
    som_from_soc = som_from_soc,
    soc_col = soc_col,
    soc_unit = if (isTRUE(som_from_soc)) (if (identical(soc_col, rvar)) rvar_unit else soc_unit) else soc_unit,
    som_factor = if (isTRUE(som_from_soc)) som_factor else NULL
  )

  # slim ONLY for output
  interval_tbl <- .esm_slim_interval(interval_full)

  # ---- Profile mass curves (use FULL; includes C(0)=0 anchor) ----
  profile_mass <- .esm_profile_mass(interval_full, core_by)

  # ---- Reference mass table (use FULL; smm_ref per stratum) ----
  ref_tbl <- .esm_define_reference_mass(
    interval_tbl    = interval_full,
    reference_mode  = reference_mode,
    reference_by    = reference_by,
    ref_trt         = ref_trt,
    custom_smm      = custom_smm,
    core_by         = core_by
  )

  # ---- Evaluate C(M_ref) per core ----
  fit_out <- .esm_fit_and_evaluate(
    profile_mass   = profile_mass,
    ref_tbl        = ref_tbl,
    core_by        = core_by,
    reference_by   = reference_by,
    interp         = interp,
    extrapolation  = extrapolation
  )

  # ---- Optional unit conversion ----
  if (out_unit == "Mg_ha") {
    interval_tbl <- .esm_add_Mg_ha(interval_tbl)
    profile_mass <- .esm_add_Mg_ha(profile_mass)
    if (is.data.frame(fit_out$cumulative))   fit_out$cumulative   <- .esm_add_Mg_ha(fit_out$cumulative)
    if (is.data.frame(fit_out$diagnostics))  fit_out$diagnostics  <- .esm_add_Mg_ha(fit_out$diagnostics)
  }

  list(
    interval     = interval_tbl,
    profile_mass = profile_mass,
    cumulative   = fit_out$cumulative,
    diagnostics  = fit_out$diagnostics,
    meta = list(
      package = "esmc",
      fn = "esm_compare",
      principle = "Depth is sampling; mineral soil mass is comparison.",
      reference_mode = reference_mode,
      reference_by = reference_by,
      ref_trt = ref_trt,
      custom_smm = custom_smm,
      interp = interp,
      rvar = rvar,
      rvar_unit = rvar_unit,
      mass_unit = out_unit,
      conversion = if (out_unit == "Mg_ha") "1 kg/m2 = 10 Mg/ha" else NULL,
      som = list(
        som_col = som_col,
        som_unit = som_unit,
        som_from_soc = som_from_soc,
        soc_col = soc_col,
        soc_unit = if (isTRUE(som_from_soc)) (if (identical(soc_col, rvar)) rvar_unit else soc_unit) else NULL,
        som_factor = if (isTRUE(som_from_soc)) som_factor else NULL
      ),
      note = "v1: spline requires >=3 points; extrapolation reserved."
    )
  )
}


# =============================================================================
# Internal helpers (shared)
# =============================================================================

.esm_add_Mg_ha <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  conv <- 10
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  kg_cols <- num_cols[grepl("_kg_m2$", num_cols)]
  if (length(kg_cols) == 0) return(df)

  for (cc in kg_cols) {
    new_name <- sub("_kg_m2$", "_Mg_ha", cc)
    if (!new_name %in% names(df)) df[[new_name]] <- df[[cc]] * conv
  }
  df
}

# Keep only user-facing columns (drop *_g_cm2, flags, etc.)
.esm_slim_interval <- function(interval_tbl) {
  keep <- c(
    "upper", "lower", "bd",
    "soil_mass_kg_m2",
    "smm_kg_m2", "cum_smm_kg_m2",
    "stock_kg_m2", "cum_stock_kg_m2"
  )
  keep <- intersect(keep, names(interval_tbl))

  # keep any id columns that are already there (core_by/reference_by/trt will be present in interval_full,
  # but this output is for the user only; we keep "whatever is not internal")
  id_cols <- setdiff(names(interval_tbl), c(
    "thick_cm","valid_thick","n","som_flag","rvar_frac","warning_code",
    "smm_g_cm2","cum_smm_g_cm2","soil_mass_g_cm2","som_frac","som_source","som_factor",
    keep
  ))

  interval_tbl[, c(id_cols, keep), drop = FALSE]
}

# Build clean C(M) points per core, with anchor C(0)=0
.esm_profile_mass <- function(interval_tbl, core_by) {
  gid <- interaction(interval_tbl[core_by], drop = TRUE)
  groups <- split(interval_tbl, gid)

  out <- lapply(groups, function(g) {
    ok <- is.finite(g$cum_smm_kg_m2) & is.finite(g$cum_stock_kg_m2)
    g2 <- g[ok, c(core_by, "cum_smm_kg_m2", "cum_stock_kg_m2"), drop = FALSE]
    g2 <- g2[order(g2$cum_smm_kg_m2), , drop = FALSE]

    if (nrow(g2) > 0 && g2$cum_smm_kg_m2[1] > 0) {
      z <- g2[1, , drop = FALSE]
      z$cum_smm_kg_m2   <- 0
      z$cum_stock_kg_m2 <- 0
      g2 <- rbind(z, g2)
    }
    g2
  })

  out <- out[!vapply(out, is.null, logical(1)) & vapply(out, nrow, integer(1)) > 0]
  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

# Interval computations (unchanged engine; produces extra fields, later slimmed)
.esm_compute_intervals <- function(
    df, rvar, core_by, rvar_unit, min_core_length_cm,
    som_col, som_unit, som_from_soc, soc_col, soc_unit, som_factor
) {
  gid <- interaction(df[core_by], drop = TRUE)
  groups <- split(df, gid)

  out <- lapply(groups, function(g) {
    g <- g[order(g$upper, g$lower), , drop = FALSE]

    g$thick_cm <- g$lower - g$upper
    g$valid_thick <- is.finite(g$thick_cm) & g$thick_cm > 0

    core_length_cm <- max(g$lower, na.rm = TRUE) - min(g$upper, na.rm = TRUE)
    meets_length   <- is.finite(core_length_cm) && core_length_cm >= min_core_length_cm

    # Total soil mass per interval (g/cm2)
    g$soil_mass_g_cm2 <- ifelse(g$valid_thick, g$bd * g$thick_cm, NA_real_)

    # SOM fraction (0-1)
    g$som_source <- "none"
    g$som_frac <- 0
    g$som_factor <- NA_real_

    if (!is.null(som_col)) {
      g$som_source <- "measured"
      som_vals <- g[[som_col]]
      g$som_frac <- switch(
        som_unit,
        pct  = som_vals / 100,
        g_kg = som_vals / 1000,
        frac = som_vals
      )
    } else if (isTRUE(som_from_soc)) {
      g$som_source <- paste0("estimated_from_", soc_col)
      g$som_factor <- som_factor

      soc_vals <- g[[soc_col]]
      soc_frac <- switch(
        soc_unit,
        pct  = soc_vals / 100,
        g_kg = soc_vals / 1000
      )
      g$som_frac <- soc_frac * som_factor
    }

    # Clamp SOM fraction to [0,1] for mineral mass calculation
    som_frac_clamped <- pmin(pmax(g$som_frac, 0), 1)

    # Mineral soil mass per interval (g/cm2)
    g$smm_g_cm2 <- g$soil_mass_g_cm2 * (1 - som_frac_clamped)
    g$cum_smm_g_cm2 <- cumsum(replace(g$smm_g_cm2, is.na(g$smm_g_cm2), 0))

    # Convert to kg/m2 (1 g/cm2 = 10 kg/m2)
    g$soil_mass_kg_m2 <- g$soil_mass_g_cm2 * 10
    g$smm_kg_m2       <- g$smm_g_cm2 * 10
    g$cum_smm_kg_m2   <- g$cum_smm_g_cm2 * 10

    # Response to fraction
    g$rvar_frac <- if (rvar_unit == "pct") g[[rvar]] / 100 else g[[rvar]] / 1000

    # Interval stock based on mineral mass
    g$stock_kg_m2 <- g$rvar_frac * g$smm_kg_m2
    g$cum_stock_kg_m2 <- cumsum(replace(g$stock_kg_m2, is.na(g$stock_kg_m2), 0))

    g$core_length_cm <- core_length_cm
    g$meets_min_core_length <- meets_length

    g
  })

  interval_tbl <- do.call(rbind, out)
  rownames(interval_tbl) <- NULL
  interval_tbl
}

# Reference mass per stratum (supports reference_by = NULL => global)
.esm_define_reference_mass <- function(interval_tbl, reference_mode, reference_by, ref_trt, custom_smm, core_by) {
  # max cumulative mineral mass per core
  by_cols <- unique(c(core_by, reference_by))
  by_cols <- by_cols[!is.null(by_cols)]
  prof_max <- stats::aggregate(
    interval_tbl$cum_smm_kg_m2,
    by = interval_tbl[by_cols],
    FUN = function(x) max(x, na.rm = TRUE)
  )
  names(prof_max)[ncol(prof_max)] <- "max_cum_smm_kg_m2"

  if (is.null(reference_by)) {
    prof_max$.__stratum__ <- "ALL"
    strata_cols <- ".__stratum__"
  } else {
    strata_cols <- reference_by
  }

  if (reference_mode == "custom") {
    ref_tbl <- unique(prof_max[strata_cols])
    ref_tbl$smm_ref_kg_m2 <- custom_smm
    return(ref_tbl)
  }

  if (reference_mode == "min") {
    ref_tbl <- stats::aggregate(
      prof_max$max_cum_smm_kg_m2,
      by = prof_max[strata_cols],
      FUN = function(x) min(x, na.rm = TRUE)
    )
    names(ref_tbl)[ncol(ref_tbl)] <- "smm_ref_kg_m2"
    return(ref_tbl)
  }

  if (reference_mode == "ref_trt") {
    if (!"trt" %in% names(prof_max)) {
      stop("reference_mode='ref_trt' requires column 'trt' to be present in core_by/data.")
    }
    prof_ref <- prof_max[prof_max$trt == ref_trt, , drop = FALSE]
    if (nrow(prof_ref) == 0) stop("No profiles found for ref_trt = ", ref_trt)

    # choose the *mean* max mass of ref_trt within stratum
    ref_tbl <- stats::aggregate(
      prof_ref$max_cum_smm_kg_m2,
      by = prof_ref[strata_cols],
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    names(ref_tbl)[ncol(ref_tbl)] <- "smm_ref_kg_m2"
    return(ref_tbl)
  }

  stop("Unknown reference_mode: ", reference_mode)
}

# Evaluate C(M_ref) per core using the profile_mass curves
.esm_fit_and_evaluate <- function(profile_mass, ref_tbl, core_by, reference_by, interp, extrapolation) {

  safe_rbind <- function(x) {
    x <- x[!vapply(x, is.null, logical(1))]
    if (length(x) == 0) return(data.frame())
    cols <- unique(unlist(lapply(x, names)))
    x2 <- lapply(x, function(d) {
      for (cc in cols) if (!cc %in% names(d)) d[[cc]] <- NA
      d[, cols, drop = FALSE]
    })
    do.call(rbind, x2)
  }

  # unify stratum handling
  if (is.null(reference_by)) {
    profile_mass$.__stratum__ <- "ALL"
    ref_tbl$.__stratum__ <- "ALL"
    strata_cols <- ".__stratum__"
  } else {
    strata_cols <- reference_by
  }

  gid <- interaction(profile_mass[core_by], drop = TRUE)
  groups <- split(profile_mass, gid)

  cum_list  <- list()
  diag_list <- list()

  for (nm in names(groups)) {
    pm <- groups[[nm]]
    pm <- pm[order(pm$cum_smm_kg_m2), , drop = FALSE]

    # stratum match -> smm_ref
    stratum_vals <- pm[1, strata_cols, drop = FALSE]
    idx <- .esm_match_row(ref_tbl[strata_cols], stratum_vals)
    if (is.na(idx)) stop("Could not match profile stratum to reference mass table.")
    smm_ref <- ref_tbl$smm_ref_kg_m2[idx]

    x <- pm$cum_smm_kg_m2
    y <- pm$cum_stock_kg_m2
    npts <- length(x)

    max_x <- if (npts > 0) max(x, na.rm = TRUE) else NA_real_
    reached <- is.finite(max_x) && is.finite(smm_ref) && max_x >= smm_ref

    esm_cum <- NA_real_
    note <- NULL

    if (!reached) {
      note <- if (extrapolation == "none") "TARGET_NOT_REACHED" else "TARGET_NOT_REACHED_EXTRAP_RESERVED"
    } else {
      if (interp == "step") {
        esm_cum <- .esm_eval_step(x, y, smm_ref)
      } else if (interp == "linear") {
        if (npts < 2) {
          note <- "INSUFFICIENT_POINTS_LINEAR"
        } else {
          esm_cum <- stats::approx(x, y, xout = smm_ref, method = "linear", rule = 1)$y
        }
      } else if (interp == "spline") {
        if (npts < 3) {
          note <- "INSUFFICIENT_POINTS_SPLINE"
        } else {
          f <- stats::splinefun(x, y, method = "natural")
          esm_cum <- f(smm_ref)
        }
      }
    }

    base_row <- pm[1, c(core_by, strata_cols), drop = FALSE]

    cum_row <- base_row
    cum_row$method <- "esm"
    cum_row$target_smm_kg_m2 <- smm_ref
    cum_row$cum_stock_kg_m2 <- esm_cum
    cum_list[[nm]] <- cum_row

    d <- base_row
    d$interp <- interp
    d$target_smm_kg_m2 <- smm_ref
    d$max_cum_smm_kg_m2 <- max_x
    d$reached_target <- reached
    d$note <- note
    diag_list[[nm]] <- d
  }

  list(
    cumulative  = safe_rbind(cum_list),
    diagnostics = safe_rbind(diag_list)
  )
}

.esm_eval_step <- function(x, y, xout) {
  ord <- order(x)
  x <- x[ord]; y <- y[ord]

  # guard: if xout is below the first knot, return first y
  if (xout <= x[1]) return(y[1])

  i <- max(which(x <= xout))
  y[i]
}

.esm_match_row <- function(tbl, rowvals) {
  if (nrow(tbl) == 0) return(NA_integer_)
  ok <- rep(TRUE, nrow(tbl))
  for (cn in names(rowvals)) ok <- ok & (tbl[[cn]] == rowvals[[cn]])
  idx <- which(ok)
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}
