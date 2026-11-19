#' Uncertainty of cumulative fixed-depth stocks with depth
#'
#' Estimates how the uncertainty of cumulative fixed-depth (FD) stocks grows
#' with depth, using replicate cores within each stratum (defined by
#' `core_by`), and suggests a stopping depth based on a signal-to-noise ratio
#' (SNR), confidence-interval overlap, or both.
#'
#' Stocks are computed as concentration times total soil mass (bulk density
#' times thickness), as in [esm_compare()]. Uncertainty of the cumulative
#' stock is computed from the cumulative stock of each core, so the positive
#' covariance among layers of the same core is taken into account (summing
#' layer variances as if layers were independent underestimates it).
#'
#' Cores must be contiguous sequences of layers starting at 0 cm, and all
#' cores of a stratum must share the same depth increments (a core may end
#' shallower than the others). Rows with missing values are removed first;
#' a core with a gap is truncated at the gap (with a warning).
#'
#' @param data A data.frame with soil layers and replicates. Must include
#'   `upper`, `lower` (cm), `bd` (g/cm3), the replicate column `rep`, the
#'   concentration column `rvar`, and all columns listed in `core_by`.
#' @param rvar Character. Name of the concentration column (e.g. SOC),
#'   interpreted using `rvar_unit`.
#' @param core_by Character vector. Columns defining the reporting stratum
#'   (e.g. `c("trt", "site", "time")`). One stopping depth is estimated per
#'   stratum.
#' @param rep Character. Column identifying replicate cores within a stratum.
#' @param rvar_unit Character. Unit of `rvar`: `"pct"` or `"g_kg"`.
#' @param method Deprecated. Uncertainty is now always computed from the
#'   per-core cumulative stocks; `"mc"` gives the same result with a warning.
#' @param n_mc,seed Deprecated and ignored (kept for compatibility).
#' @param uncertainty Character. Measure used for the SNR and the intervals:
#'   `"sd"` (standard deviation among cores, default; describes the spread of
#'   single cores) or `"se"` (standard error of the mean, SD / sqrt(n); describes
#'   the precision of the stratum mean). Both are returned.
#' @param k Numeric. Multiplier for the intervals (mean +/- k * uncertainty).
#' @param criterion Character. Stopping criterion: `"snr"`, `"ci_overlap"`, or
#'   `"both"` (the shallowest of the two).
#' @param snr_threshold Numeric. Stop at the first layer where SNR < threshold.
#' @param overlap_threshold Numeric. Stop at the first layer where the overlap
#'   length exceeds this threshold (0 flags any overlap).
#' @param min_reps Integer. Minimum number of cores per layer. Layers with
#'   fewer cores are flagged (`rep_ok = FALSE`) and are not used by the
#'   stopping rule; with fewer than 2 cores the uncertainty is NA.
#' @param som_col,som_unit,som_from_soc,soc_col,soc_unit,som_factor
#'   Deprecated and ignored: FD stocks use total soil mass, so SOM is not
#'   needed. Kept for compatibility.
#' @param out_unit Character. `"kg_m2"` (default) or `"Mg_ha"` (adds columns
#'   with suffix `_Mg_ha`; 1 kg/m2 = 10 Mg/ha).
#'
#' @return A list with `depth` (depth-wise diagnostics per stratum), `summary`
#'   (one row per stratum with the suggested `stop_depth_cm`) and `meta`.
#' @export
esm_uncert_depth <- function(
    data,
    rvar,
    core_by,
    rep = "rep",
    rvar_unit = c("pct", "g_kg"),
    method = c("analytic", "mc"),
    n_mc = NULL,
    uncertainty = c("sd", "se"),
    k = 1.96,
    criterion = c("snr", "ci_overlap", "both"),
    snr_threshold = 1,
    overlap_threshold = 0,
    seed = NULL,
    min_reps = 2,
    som_col = NULL,
    som_unit = c("pct", "g_kg", "frac"),
    som_from_soc = FALSE,
    soc_col = NULL,
    soc_unit = c("pct", "g_kg"),
    som_factor = 2,
    out_unit = c("kg_m2", "Mg_ha")
) {
  rvar_unit   <- match.arg(rvar_unit)
  method      <- match.arg(method)
  uncertainty <- match.arg(uncertainty)
  criterion   <- match.arg(criterion)
  out_unit    <- match.arg(out_unit)

  if (method == "mc")
    warning("method = 'mc' is deprecated: uncertainty is computed from the ",
            "per-core cumulative stocks (no simulation needed).", call. = FALSE)
  if (!is.null(som_col) || isTRUE(som_from_soc))
    warning("SOM arguments are ignored: fixed-depth stocks use total soil ",
            "mass, so SOM is not needed.", call. = FALSE)

  req  <- c("upper", "lower", "bd", rep, rvar, core_by)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  df <- data[is.finite(data$upper) & is.finite(data$lower) &
             is.finite(data$bd) & is.finite(data[[rvar]]), , drop = FALSE]
  if (!nrow(df)) stop("No rows left after removing NAs in upper/lower/bd/rvar.")

  # Same geometry checks as esm_compare(): contiguous from 0, no overlaps
  df <- .esm_check_profiles(df, unique(c(core_by, rep)))

  conv <- if (rvar_unit == "pct") 1 / 100 else 1 / 1000
  # stock (kg/m2) = fraction x BD (g/cm3) x thickness (cm) x 10
  df$stock_layer_kg_m2 <- df[[rvar]] * conv * df$bd * (df$lower - df$upper) * 10

  gid    <- interaction(df[core_by], drop = TRUE)
  groups <- split(df, gid)
  depth_list <- list(); summ_list <- list()

  for (nm in names(groups)) {
    g    <- groups[[nm]]
    cid  <- as.character(g[[rep]])
    lay  <- unique(g[, c("upper", "lower")])
    lay  <- lay[order(lay$upper), , drop = FALSE]
    # cores must share increments (a shallower core is a prefix of the others)
    if (any(duplicated(lay$upper)) ||
        (nrow(lay) > 1 && any(lay$upper[-1] != lay$lower[-nrow(lay)])))
      stop("Cores in stratum '", nm, "' do not share the same depth ",
           "increments.", call. = FALSE)
    L <- nrow(lay)
    if (L < 2) next

    # core x layer matrix of layer stocks; cumulative per core
    cores <- unique(cid)
    S <- matrix(NA_real_, nrow = length(cores), ncol = L)
    S[cbind(match(cid, cores), match(g$upper, lay$upper))] <- g$stock_layer_kg_m2
    C <- t(apply(S, 1, cumsum))

    n_rep  <- colSums(is.finite(S))
    rep_ok <- n_rep >= min_reps
    sd_fin <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else stats::sd(x) }
    inc_mu <- colMeans(S, na.rm = TRUE)
    inc_sd <- apply(S, 2, sd_fin)
    cum_mu <- colMeans(C, na.rm = TRUE)
    cum_sd <- apply(C, 2, sd_fin)
    inc_se <- inc_sd / sqrt(n_rep)
    cum_se <- cum_sd / sqrt(colSums(is.finite(C)))

    inc_u <- if (uncertainty == "sd") inc_sd else inc_se
    cum_u <- if (uncertainty == "sd") cum_sd else cum_se
    cum_prev_mu <- c(NA, cum_mu[-L])
    cum_prev_u  <- c(NA, cum_u[-L])

    snr <- inc_mu / sqrt(cum_prev_u^2 + inc_u^2)
    snr[1] <- NA
    snr[!rep_ok] <- NA

    overlap_len <- pmax(0, pmin(cum_prev_mu + k * cum_prev_u, inc_mu + k * inc_u) -
                          pmax(cum_prev_mu - k * cum_prev_u, inc_mu - k * inc_u))
    overlap_any <- overlap_len > overlap_threshold
    overlap_any[1] <- NA
    overlap_any[!rep_ok] <- NA

    stop_idx <- NA_integer_
    if (criterion %in% c("snr", "both")) {
      idx <- which(!is.na(snr) & snr < snr_threshold)
      if (length(idx)) stop_idx <- idx[1]
    }
    if (criterion %in% c("ci_overlap", "both")) {
      idx <- which(!is.na(overlap_any) & overlap_any)
      if (length(idx)) stop_idx <- if (is.na(stop_idx)) idx[1] else min(stop_idx, idx[1])
    }
    stop_depth_cm <- if (!is.na(stop_idx)) lay$lower[stop_idx] else NA_real_

    d0 <- g[1, core_by, drop = FALSE]
    depth_tbl <- cbind(
      d0[rep(1, L), , drop = FALSE], lay,
      n_rep = n_rep, rep_ok = rep_ok,
      stock_layer_mu_kg_m2 = inc_mu, stock_layer_sd_kg_m2 = inc_sd,
      stock_layer_se_kg_m2 = inc_se,
      cum_stock_mu_kg_m2   = cum_mu, cum_stock_sd_kg_m2   = cum_sd,
      cum_stock_se_kg_m2   = cum_se,
      snr = snr, ci_overlap_len = overlap_len, ci_overlap_any = overlap_any,
      stop_depth_cm = stop_depth_cm, stop_idx = stop_idx
    )
    rownames(depth_tbl) <- NULL
    if (out_unit == "Mg_ha") {
      for (cc in grep("_kg_m2$", names(depth_tbl), value = TRUE))
        depth_tbl[[sub("_kg_m2$", "_Mg_ha", cc)]] <- depth_tbl[[cc]] * 10
    }
    depth_list[[nm]] <- depth_tbl

    srow <- g[1, core_by, drop = FALSE]
    srow$stop_depth_cm <- stop_depth_cm
    srow$stop_idx <- stop_idx
    srow$any_rep_insufficient <- any(!rep_ok)
    summ_list[[nm]] <- srow
  }

  depth_out <- if (length(depth_list)) do.call(rbind, depth_list) else data.frame()
  summ_out  <- if (length(summ_list))  do.call(rbind, summ_list)  else data.frame()
  rownames(depth_out) <- NULL; rownames(summ_out) <- NULL

  list(
    depth = depth_out,
    summary = summ_out,
    meta = list(
      mode = "fd", rvar = rvar, rvar_unit = rvar_unit,
      uncertainty = uncertainty, k = k, criterion = criterion,
      snr_threshold = snr_threshold, overlap_threshold = overlap_threshold,
      rep = rep, min_reps = min_reps, out_unit = out_unit,
      stock_basis = "total soil mass",
      conversion = "1 kg/m2 = 10 Mg/ha"
    )
  )
}
