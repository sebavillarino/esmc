# =============================================================================
# esmc — Equivalent Soil Mass
# Layer-wise comparison: esm_compare_layers()
#
# Principle:
#   Depth is sampling; mineral soil mass is comparison.
#
# Computes layer stocks as:
#   Stock_layer = C(M_lower) - C(M_upper)
# where C(M) is cumulative stock evaluated on a continuous C(M) curve.
#
# Design (aligned with esm_compare()):
# - interval_full used for ALL internal computations
# - interval output is slimmed for the user
# - reference_by can be NULL (global targets)
# =============================================================================


#' ESM comparison by soil layer (stocks)
#'
#' Computes layer stocks (e.g. SOC stocks) under an Equivalent Soil Mass (ESM)
#' framework using the soil layers defined in the input data.
#'
#' @param data data.frame with required columns: upper, lower, bd, trt, rvar,
#'   and grouping columns in core_by (and reference_by if not NULL).
#' @param rvar string. Concentration column name (e.g. "soc").
#' @param core_by character vector defining a soil core/profile
#'   (e.g. c("Site","yr","trt","rep")).
#' @param reference_by character vector defining strata within which ESM targets
#'   are computed (e.g. c("Site","yr")), or NULL for global targets.
#' @param reference_mode "min" | "ref_trt" | "custom".
#' @param ref_trt value in column `trt` used as reference treatment
#'   (required if reference_mode = "ref_trt").
#' @param custom_smm custom target mineral masses (required if reference_mode = "custom").
#'   Either:
#'     - named numeric vector: names are boundary depths (lower), values are target cum_smm_kg_m2
#'     - data.frame with columns: lower, target_smm_kg_m2 (and optionally reference_by cols)
#' @param interp interpolation method for C(M): "step", "linear", or "spline".
#' @param rvar_unit "pct" | "g_kg".
#' @param som_col optional measured SOM column.
#' @param som_unit "pct" | "g_kg" | "frac".
#' @param som_from_soc logical. If TRUE, estimate SOM from SOC using som_factor.
#' @param soc_col SOC column used for SOM estimation (defaults to rvar).
#' @param soc_unit "pct" | "g_kg".
#' @param som_factor numeric conversion factor (default 2).
#' @param min_core_length_cm numeric. Minimum profile length required.
#' @param out_unit "kg_m2" (default) or "Mg_ha" (adds *_Mg_ha columns).
#'
#' @return list with elements:
#' \itemize{
#'   \item interval: cleaned per-layer computations (slim output)
#'   \item profile_mass: cumulative stock vs cumulative mineral mass per core
#'   \item layers: layer definitions
#'   \item targets: target mineral mass per boundary
#'   \item layer: ESM stock per layer and core
#'   \item meta: metadata
#' }
#' @export
esm_compare_layers <- function(
    data,
    rvar,
    core_by,
    reference_by = NULL,

    reference_mode = c("min", "ref_trt", "custom"),
    ref_trt = NULL,
    custom_smm = NULL,

    interp = c("step", "linear", "spline"),
    rvar_unit = c("pct", "g_kg"),

    som_col = NULL,
    som_unit = c("pct", "g_kg", "frac"),
    som_from_soc = FALSE,
    soc_col = NULL,
    soc_unit = c("pct", "g_kg"),
    som_factor = 2,

    min_core_length_cm = 0,
    out_unit = c("kg_m2", "Mg_ha")
) {

  reference_mode <- match.arg(reference_mode)
  interp    <- match.arg(interp)
  rvar_unit <- match.arg(rvar_unit)
  som_unit  <- match.arg(som_unit)
  soc_unit  <- match.arg(soc_unit)
  out_unit  <- match.arg(out_unit)

  # ---------------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------------
  if (reference_mode == "ref_trt" && is.null(ref_trt)) {
    stop("reference_mode='ref_trt' requires ref_trt.")
  }
  if (reference_mode == "custom" && is.null(custom_smm)) {
    stop("reference_mode='custom' requires custom_smm.")
  }

  needed <- unique(c("upper", "lower", "bd", "trt", rvar, core_by, reference_by))
  missing <- setdiff(needed, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  if (isTRUE(som_from_soc)) {
    if (is.null(soc_col)) soc_col <- rvar
    if (!soc_col %in% names(data)) stop("soc_col not found: ", soc_col)
    if (!is.finite(som_factor) || som_factor <= 0) stop("som_factor must be > 0.")
  } else {
    soc_col <- NULL
  }

  # Drop essential NAs
  df <- data[
    is.finite(data$upper) &
      is.finite(data$lower) &
      is.finite(data$bd) &
      is.finite(data[[rvar]]),
    ,
    drop = FALSE
  ]
  if (nrow(df) == 0) stop("No valid rows after removing missing values in upper/lower/bd/rvar.")

  # ---------------------------------------------------------------------------
  # Interval computations (FULL for internal use)
  # ---------------------------------------------------------------------------
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

  # Slim only for output
  interval_tbl <- .esm_slim_interval(interval_full)

  # ---------------------------------------------------------------------------
  # Profile-level cumulative C(M) (use FULL)
  # ---------------------------------------------------------------------------
  profile_mass <- .esm_profile_mass(interval_full, core_by)

  # ---------------------------------------------------------------------------
  # Layer definitions (use FULL)
  # ---------------------------------------------------------------------------
  layers_tbl <- .esm_layers_from_data(interval_full, reference_by)

  # ---------------------------------------------------------------------------
  # Target mineral mass per boundary (use FULL)
  # ---------------------------------------------------------------------------
  targets_tbl <- .esm_targets_by_boundary(
    interval_tbl   = interval_full,
    core_by        = core_by,
    reference_by   = reference_by,
    reference_mode = reference_mode,
    ref_trt        = ref_trt,
    custom_smm     = custom_smm
  )

  # ---------------------------------------------------------------------------
  # Evaluate C(M) at boundaries and compute layer stocks
  # ---------------------------------------------------------------------------
  layer_tbl <- .esm_layer_stocks_from_targets(
    profile_mass = profile_mass,
    layers_tbl   = layers_tbl,
    targets_tbl  = targets_tbl,
    core_by      = core_by,
    reference_by = reference_by,
    interp       = interp
  )

  # ---------------------------------------------------------------------------
  # Optional unit conversion
  # ---------------------------------------------------------------------------
  if (out_unit == "Mg_ha") {
    interval_tbl  <- .esm_add_Mg_ha(interval_tbl)
    profile_mass  <- .esm_add_Mg_ha(profile_mass)
    targets_tbl   <- .esm_add_Mg_ha(targets_tbl)
    layer_tbl     <- .esm_add_Mg_ha(layer_tbl)
  }

  list(
    interval     = interval_tbl,
    profile_mass = profile_mass,
    layers       = layers_tbl,
    targets      = targets_tbl,
    layer        = layer_tbl,
    meta = list(
      package = "esmc",
      fn = "esm_compare_layers",
      principle = "Depth is sampling; mineral soil mass is comparison.",
      rvar = rvar,
      rvar_unit = rvar_unit,
      interp = interp,
      reference_by = reference_by,
      core_by = core_by,
      reference_mode = reference_mode,
      ref_trt = ref_trt,
      out_unit = out_unit
    )
  )
}


# =============================================================================
# Internal helpers (esm_compare_layers-specific)
# =============================================================================

# Layer definitions from the sampled layers in the data.
# If reference_by is NULL -> global layer set (upper/lower pairs).
# If reference_by provided -> layer set per stratum (reference_by + upper/lower).
.esm_layers_from_data <- function(interval_tbl, reference_by) {
  if (is.null(reference_by)) {
    out <- unique(interval_tbl[, c("upper", "lower"), drop = FALSE])
    out <- out[order(out$upper, out$lower), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }

  key <- interaction(interval_tbl[reference_by], drop = TRUE)
  groups <- split(interval_tbl, key)

  out <- lapply(groups, function(g) {
    unique(g[, c(reference_by, "upper", "lower"), drop = FALSE])
  })

  out <- do.call(rbind, out)
  out <- out[do.call(order, out[, c(reference_by, "upper", "lower"), drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Targets by boundary (lower) — supports:
# - reference_by = NULL => global targets across all cores
# - reference_mode:
#     "min"     => minimum cum_smm reached at each boundary
#     "ref_trt" => mean cum_smm of ref_trt at each boundary
#     "custom"  => user supplied targets
.esm_targets_by_boundary <- function(
    interval_tbl,
    core_by,
    reference_by = NULL,
    reference_mode = c("min", "ref_trt", "custom"),
    ref_trt = NULL,
    custom_smm = NULL
) {
  reference_mode <- match.arg(reference_mode)

  bnds <- sort(unique(interval_tbl$lower))

  if (reference_mode == "custom") {
    if (is.numeric(custom_smm) && !is.null(names(custom_smm))) {
      out <- data.frame(
        lower = as.numeric(names(custom_smm)),
        target_smm_kg_m2 = as.numeric(custom_smm),
        stringsAsFactors = FALSE
      )
      out <- out[order(out$lower), , drop = FALSE]
      rownames(out) <- NULL
      return(out)
    }

    if (is.data.frame(custom_smm)) {
      needed <- c("lower", "target_smm_kg_m2")
      miss <- setdiff(needed, names(custom_smm))
      if (length(miss) > 0) stop("custom_smm is missing: ", paste(miss, collapse = ", "))
      if (!is.null(reference_by)) {
        miss2 <- setdiff(reference_by, names(custom_smm))
        if (length(miss2) > 0) stop("custom_smm is missing reference_by cols: ", paste(miss2, collapse = ", "))
      }
      out <- custom_smm
      out$lower <- as.numeric(out$lower)
      out$target_smm_kg_m2 <- as.numeric(out$target_smm_kg_m2)
      out <- out[do.call(order, out[, c(reference_by, "lower"), drop = FALSE]), , drop = FALSE]
      rownames(out) <- NULL
      return(out)
    }

    stop("custom_smm must be a named numeric vector or a data.frame.")
  }

  # Helper builders
  build_global <- function(fun_boundary) {
    out <- lapply(bnds, function(b) {
      sub <- interval_tbl[interval_tbl$lower == b, , drop = FALSE]
      data.frame(lower = b, target_smm_kg_m2 = fun_boundary(sub), stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }

  build_by_group <- function(gcols, fun_boundary) {
    key <- interaction(interval_tbl[gcols], drop = TRUE)
    groups <- split(interval_tbl, key)

    out <- lapply(groups, function(g) {
      base <- unique(g[, gcols, drop = FALSE])[1, , drop = FALSE]
      tt <- lapply(bnds, function(b) {
        sub <- g[g$lower == b, , drop = FALSE]
        data.frame(lower = b, target_smm_kg_m2 = fun_boundary(sub), stringsAsFactors = FALSE)
      })
      cbind(base, do.call(rbind, tt))
    })
    do.call(rbind, out)
  }

  if (reference_mode == "min") {
    fun_boundary <- function(sub) min(sub$cum_smm_kg_m2, na.rm = TRUE)
  } else { # ref_trt
    if (is.null(ref_trt)) stop("reference_mode='ref_trt' requires ref_trt.")
    fun_boundary <- function(sub) {
      sub2 <- sub[sub$trt == ref_trt, , drop = FALSE]
      mean(sub2$cum_smm_kg_m2, na.rm = TRUE)
    }
  }

  if (is.null(reference_by)) {
    out <- build_global(fun_boundary)
  } else {
    out <- build_by_group(reference_by, fun_boundary)
  }

  out <- out[is.finite(out$target_smm_kg_m2) & out$target_smm_kg_m2 >= 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Compute layer stocks using targets + C(M) curves.
.esm_layer_stocks_from_targets <- function(
    profile_mass,
    layers_tbl,
    targets_tbl,
    core_by,
    reference_by = NULL,
    interp = c("step", "linear", "spline")
) {
  interp <- match.arg(interp)

  # add a 0-boundary target (M=0) per stratum
  targets0 <- targets_tbl
  if (is.null(reference_by)) {
    if (!any(targets0$lower == 0)) {
      targets0 <- rbind(data.frame(lower = 0, target_smm_kg_m2 = 0), targets0)
    }
  } else {
    strata <- unique(targets0[, reference_by, drop = FALSE])
    add0 <- cbind(strata, data.frame(lower = 0, target_smm_kg_m2 = 0))
    targets0 <- rbind(add0, targets0)
  }

  eval_C <- function(pm, M) {
    M <- as.numeric(M)
    if (!is.finite(M)) return(NA_real_)

    pm <- pm[order(pm$cum_smm_kg_m2), , drop = FALSE]
    x <- pm$cum_smm_kg_m2
    y <- pm$cum_stock_kg_m2

    if (length(x) == 0) return(NA_real_)
    if (M <= x[1]) return(y[1])
    if (M >= x[length(x)]) return(y[length(y)])

    if (interp == "step") {
      idx <- max(which(x <= M))
      return(y[idx])
    }
    if (interp == "linear") {
      return(stats::approx(x = x, y = y, xout = M, method = "linear", rule = 2)$y)
    }

    # spline
    if (length(unique(x)) < 3) {
      return(stats::approx(x = x, y = y, xout = M, method = "linear", rule = 2)$y)
    }
    stats::spline(x = x, y = y, xout = M, method = "natural")$y
  }

  gid <- interaction(profile_mass[core_by], drop = TRUE)
  pm_groups <- split(profile_mass, gid)

  out_list <- vector("list", length(pm_groups))
  names(out_list) <- names(pm_groups)

  for (i in seq_along(pm_groups)) {
    pm <- pm_groups[[i]]
    core_id <- pm[1, core_by, drop = FALSE]

    # pick the relevant stratum targets + layers
    if (is.null(reference_by)) {
      tsub <- targets0
      lyr  <- layers_tbl
    } else {
      key <- core_id[, reference_by, drop = FALSE]

      tsub <- targets0
      for (cc in reference_by) tsub <- tsub[tsub[[cc]] == key[[cc]], , drop = FALSE]

      lyr <- layers_tbl
      for (cc in reference_by) lyr <- lyr[lyr[[cc]] == key[[cc]], , drop = FALSE]
    }

    # boundary -> target mass map
    m_map <- setNames(tsub$target_smm_kg_m2, as.character(tsub$lower))

    res <- lapply(seq_len(nrow(lyr)), function(j) {
      u <- lyr$upper[j]
      l <- lyr$lower[j]

      M_u <- if (as.character(u) %in% names(m_map)) m_map[[as.character(u)]] else NA_real_
      M_l <- if (as.character(l) %in% names(m_map)) m_map[[as.character(l)]] else NA_real_

      C_u <- eval_C(pm, M_u)
      C_l <- eval_C(pm, M_l)

      data.frame(
        core_id,
        upper = u,
        lower = l,
        target_smm_upper_kg_m2 = M_u,
        target_smm_lower_kg_m2 = M_l,
        cum_upper_kg_m2 = C_u,
        cum_lower_kg_m2 = C_l,
        layer_stock_kg_m2 = C_l - C_u,
        stringsAsFactors = FALSE
      )
    })

    out_list[[i]] <- do.call(rbind, res)
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
