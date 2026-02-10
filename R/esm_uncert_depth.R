#' Uncertainty propagation with depth (FD) using replicates
#'
#' Estimates uncertainty propagation of cumulative stocks with depth using
#' replicate cores within each stratum (defined by `core_by`). Returns depth-wise
#' diagnostics and a suggested stopping depth based on signal-to-noise ratio (SNR),
#' confidence-interval overlap, or both.
#'
#' By default, this function uses fixed-depth (FD) intervals defined by `upper` and
#' `lower`. Stocks are computed on a mineral-soil-mass basis by subtracting SOM mass
#' from total soil mass, using either measured SOM (`som_col`) or optional estimation
#' from SOC (`som_from_soc`).
#'
#' @param data A data.frame containing soil layers and replicates. Must include:
#'   `upper`, `lower` (cm), `bd` (g/cm3), replicate column `rep`, concentration column
#'   `rvar`, and all columns listed in `core_by`. If SOM correction is used, must also
#'   include `som_col` (measured) or `soc_col` (for estimation).
#' @param rvar Character string. Name of the concentration column (e.g., SOC). Interpreted
#'   using `rvar_unit`.
#' @param core_by Character vector. Columns defining the reporting stratum (e.g.,
#'   `c("trt","site","time")`). A single stopping depth is estimated per unique stratum.
#' @param rep Character string. Column name identifying replicate cores within each stratum.
#'   Default `"rep"`.
#' @param rvar_unit Character. Unit of `rvar`: `"pct"` for percent (0–100) or `"g_kg"` for g/kg.
#' @param method Character. Uncertainty propagation method: `"mc"` (Monte Carlo) or `"analytic"`.
#' @param n_mc Integer. Number of Monte Carlo draws per stratum (only used when `method="mc"`).
#' @param k Numeric. Confidence multiplier (e.g., 1.96 for ~95% intervals).
#' @param criterion Character. Stopping-depth criterion: `"snr"`, `"ci_overlap"`, or `"both"`.
#' @param snr_threshold Numeric. Stopping rule for SNR: stop at first layer where SNR < threshold.
#' @param overlap_threshold Numeric. Stopping rule for CI overlap: stop at first layer where the
#'   overlap length exceeds this threshold. Use 0 to flag any overlap.
#' @param seed Integer. Random seed for Monte Carlo draws.
#'
#' @param min_reps Integer. Minimum number of unique replicates required per depth interval.
#'   If fewer are available, the layer is flagged and its SD is set to 0 (uncertainty likely
#'   underestimated).
#'
#' @param som_col Optional character string. Name of measured SOM column (mass basis).
#'   If provided, SOM is converted using `som_unit` and subtracted from total soil mass.
#' @param som_unit Character. Unit of `som_col`: `"pct"` (percent), `"g_kg"` (g/kg), or `"frac"` (0–1).
#' @param som_from_soc Logical. If TRUE and `som_col` is NULL, estimate SOM from SOC using
#'   `SOM = SOC * som_factor`. If TRUE and `soc_col` is NULL, `soc_col` defaults to `rvar`.
#' @param soc_col Optional character string. Name of SOC column used for SOM estimation when
#'   `som_from_soc=TRUE`. Defaults to `rvar` if NULL.
#' @param soc_unit Character. Unit of `soc_col` when `soc_col` differs from `rvar`: `"pct"` or `"g_kg"`.
#' @param som_factor Numeric > 0. Conversion factor used when estimating SOM from SOC
#'   (default 2).
#'
#' @param out_unit Character. Output unit for convenience columns: `"kg_m2"` (default) or `"Mg_ha"`.
#'   When `"Mg_ha"`, additional columns with suffix `_Mg_ha` are returned (1 kg/m2 = 10 Mg/ha).
#'
#' @return A list with:
#'   \itemize{
#'     \item `depth`: depth-wise diagnostics per stratum (means, SDs, SNR, CI overlap, replicate counts).
#'     \item `summary`: one row per stratum with suggested `stop_depth_cm` and flags.
#'     \item `meta`: settings used.
#'   }
#' @export
esm_uncert_depth <- function(
    data,
    rvar,
    core_by,
    rep = "rep",
    rvar_unit = c("pct","g_kg"),
    method = c("mc","analytic"),
    n_mc = 1000,
    k = 1.96,
    criterion = c("snr","ci_overlap","both"),
    snr_threshold = 1,
    overlap_threshold = 0,
    seed = 1,

    # replicate QC
    min_reps = 2,

    # SOM options (aligned with esm_compare)
    som_col        = NULL,
    som_unit       = c("pct","g_kg","frac"),
    som_from_soc   = FALSE,
    soc_col        = NULL,
    soc_unit       = c("pct","g_kg"),
    som_factor     = 2,

    # output units
    out_unit       = c("kg_m2","Mg_ha")
) {
  rvar_unit <- match.arg(rvar_unit)
  method <- match.arg(method)
  criterion <- match.arg(criterion)
  som_unit <- match.arg(som_unit)
  soc_unit <- match.arg(soc_unit)
  out_unit <- match.arg(out_unit)

  req <- c("upper","lower","bd", rep, rvar, core_by)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))

  if (!is.null(som_col) && !som_col %in% names(data)) stop("som_col not found in data: ", som_col)

  if (isTRUE(som_from_soc)) {
    if (is.null(soc_col)) soc_col <- rvar
    if (!soc_col %in% names(data)) stop("soc_col not found in data: ", soc_col)
    if (!is.finite(som_factor) || som_factor <= 0) stop("som_factor must be > 0.")
  } else {
    soc_col <- NULL
  }

  df <- data[!is.na(data$upper) & !is.na(data$lower) & !is.na(data$bd) & !is.na(data[[rvar]]), , drop=FALSE]
  if (!nrow(df)) stop("No rows left after removing NAs in upper/lower/bd/rvar.")

  df$thick_cm <- df$lower - df$upper
  if (any(df$thick_cm <= 0, na.rm=TRUE)) stop("Found non-positive thickness (lower - upper).")

  # convert concentration to fraction
  df$rvar_frac <- if (rvar_unit == "pct") df[[rvar]]/100 else df[[rvar]]/1000

  # FD soil mass (kg/m2): bd(g/cm3)*thick(cm)=g/cm2; *10 => kg/m2
  df$soil_mass_kg_m2 <- df$bd * df$thick_cm * 10

  # SOM fraction (0-1)
  df$som_source <- "none"
  df$som_frac <- 0
  df$som_factor_used <- NA_real_

  if (!is.null(som_col)) {
    df$som_source <- "measured"
    som_vals <- df[[som_col]]
    df$som_frac <- switch(
      som_unit,
      pct  = som_vals/100,
      g_kg = som_vals/1000,
      frac = som_vals
    )
  } else if (isTRUE(som_from_soc)) {
    df$som_source <- paste0("estimated_from_", soc_col)
    df$som_factor_used <- som_factor

    soc_vals <- df[[soc_col]]
    # soc_unit inherits rvar_unit if soc_col == rvar
    soc_unit_eff <- if (identical(soc_col, rvar)) rvar_unit else soc_unit
    soc_frac <- switch(
      soc_unit_eff,
      pct  = soc_vals/100,
      g_kg = soc_vals/1000
    )
    df$som_frac <- soc_frac * som_factor
  }

  # flags + clamp
  df$som_flag <- NA_character_
  flags <- character(0)
  if (any(!is.finite(df$som_frac), na.rm=TRUE)) flags <- c(flags, "SOM_NONFINITE")
  if (any(df$som_frac < 0, na.rm=TRUE)) flags <- c(flags, "SOM_NEGATIVE")
  if (any(df$som_frac > 1, na.rm=TRUE)) flags <- c(flags, "SOM_GT1")
  if (any(df$som_frac > 0.3, na.rm=TRUE)) flags <- c(flags, "SOM_GT0.3")
  if (length(flags)) df$som_flag <- paste(unique(flags), collapse=";")

  som_frac_clamped <- pmin(pmax(df$som_frac, 0), 1)

  # mineral mass + stock per replicate
  df$mineral_mass_kg_m2 <- df$soil_mass_kg_m2 * (1 - som_frac_clamped)
  df$stock_layer_kg_m2 <- df$rvar_frac * df$mineral_mass_kg_m2

  # group by stratum core_by
  gid <- interaction(df[core_by], drop = TRUE)
  groups <- split(df, gid)

  depth_list <- list()
  summ_list  <- list()

  for (nm in names(groups)) {
    g <- groups[[nm]]

    layers <- unique(g[, c("upper","lower"), drop=FALSE])
    layers <- layers[order(layers$upper, layers$lower), , drop=FALSE]
    L <- nrow(layers)
    if (L < 2) next

    mu_S <- sd_S <- numeric(L)
    n_rep <- integer(L)
    rep_ok <- logical(L)

    for (i in seq_len(L)) {
      ui <- layers$upper[i]; li <- layers$lower[i]
      gi <- g[g$upper == ui & g$lower == li, , drop=FALSE]

      # replicate stocks (drop NA)
      s <- gi$stock_layer_kg_m2
      reps <- gi[[rep]]
      ok <- is.finite(s) & !is.na(reps)

      n_rep[i] <- length(unique(reps[ok]))
      rep_ok[i] <- n_rep[i] >= min_reps

      mu_S[i] <- mean(s[ok], na.rm=TRUE)
      sd_S[i] <- stats::sd(s[ok], na.rm=TRUE)
      if (!is.finite(sd_S[i])) sd_S[i] <- 0
      if (!rep_ok[i]) {
        # keep running but flag: uncertainty likely underestimated
        sd_S[i] <- 0
      }
    }

    if (method == "analytic") {
      cum_mu <- cumsum(mu_S)
      cum_sd <- sqrt(cumsum(sd_S^2))
      inc_mu <- mu_S
      inc_sd <- sd_S
    } else {
      set.seed(seed)
      S_draw <- matrix(0, nrow=n_mc, ncol=L)
      for (i in seq_len(L)) S_draw[, i] <- stats::rnorm(n_mc, mean = mu_S[i], sd = sd_S[i])
      S_draw[S_draw < 0] <- 0

      cum_draw <- t(apply(S_draw, 1, cumsum))
      cum_mu <- apply(cum_draw, 2, mean)
      cum_sd <- apply(cum_draw, 2, stats::sd)

      inc_mu <- apply(S_draw, 2, mean)
      inc_sd <- apply(S_draw, 2, stats::sd)
    }

    cum_prev_mu <- c(NA, cum_mu[-L])
    cum_prev_sd <- c(NA, cum_sd[-L])

    snr <- inc_mu / sqrt(pmax(cum_prev_sd^2 + inc_sd^2, 0))
    snr[1] <- NA

    A_lo <- cum_prev_mu - k*cum_prev_sd
    A_hi <- cum_prev_mu + k*cum_prev_sd
    B_lo <- inc_mu - k*inc_sd
    B_hi <- inc_mu + k*inc_sd

    overlap_len <- pmax(0, pmin(A_hi, B_hi) - pmax(A_lo, B_lo))
    overlap_any <- overlap_len > overlap_threshold
    overlap_any[1] <- NA

    stop_idx <- NA_integer_
    if (criterion %in% c("snr","both")) {
      idx <- which(!is.na(snr) & snr < snr_threshold)
      if (length(idx)) stop_idx <- idx[1]
    }
    if (criterion %in% c("ci_overlap","both")) {
      idx <- which(!is.na(overlap_any) & overlap_any)
      if (length(idx)) stop_idx <- if (is.na(stop_idx)) idx[1] else min(stop_idx, idx[1])
    }
    stop_depth_cm <- if (!is.na(stop_idx)) layers$lower[stop_idx] else NA_real_

    d0 <- g[1, core_by, drop=FALSE]
    depth_tbl <- cbind(
      d0[rep(1, L), , drop=FALSE],
      layers,
      n_rep = n_rep,
      rep_ok = rep_ok,
      stock_layer_mu_kg_m2 = inc_mu,
      stock_layer_sd_kg_m2 = inc_sd,
      cum_stock_mu_kg_m2   = cum_mu,
      cum_stock_sd_kg_m2   = cum_sd,
      snr                  = snr,
      ci_overlap_len       = overlap_len,
      ci_overlap_any       = overlap_any,
      stop_depth_cm        = stop_depth_cm,
      stop_idx             = stop_idx
    )
    rownames(depth_tbl) <- NULL

    # Optional unit conversion columns
    if (out_unit == "Mg_ha") {
      conv <- 10
      depth_tbl$stock_layer_mu_Mg_ha <- depth_tbl$stock_layer_mu_kg_m2 * conv
      depth_tbl$stock_layer_sd_Mg_ha <- depth_tbl$stock_layer_sd_kg_m2 * conv
      depth_tbl$cum_stock_mu_Mg_ha   <- depth_tbl$cum_stock_mu_kg_m2   * conv
      depth_tbl$cum_stock_sd_Mg_ha   <- depth_tbl$cum_stock_sd_kg_m2   * conv
    }

    depth_list[[nm]] <- depth_tbl

    srow <- g[1, core_by, drop=FALSE]
    srow$stop_depth_cm <- stop_depth_cm
    srow$stop_idx <- stop_idx
    srow$any_rep_insufficient <- any(!rep_ok, na.rm=TRUE)
    srow$som_source <- {
      tab <- sort(table(g$som_source), decreasing = TRUE)
      names(tab)[1]
    }
    srow$som_factor <- unique(g$som_factor_used)[1]
    summ_list[[nm]] <- srow
  }

  depth_out <- if (length(depth_list)) do.call(rbind, depth_list) else data.frame()
  summ_out  <- if (length(summ_list))  do.call(rbind, summ_list)  else data.frame()

  list(
    depth = depth_out,
    summary = summ_out,
    meta = list(
      mode = "fd",
      rvar = rvar,
      rvar_unit = rvar_unit,
      method = method,
      n_mc = if (method == "mc") n_mc else NA_integer_,
      k = k,
      criterion = criterion,
      snr_threshold = snr_threshold,
      overlap_threshold = overlap_threshold,
      rep = rep,
      min_reps = min_reps,
      out_unit = out_unit,
      conversion = "1 kg/m2 = 10 Mg/ha",
      som = list(
        som_col = som_col,
        som_unit = som_unit,
        som_from_soc = som_from_soc,
        soc_col = soc_col,
        soc_unit = if (isTRUE(som_from_soc)) (if (identical(soc_col, rvar)) rvar_unit else soc_unit) else NULL,
        som_factor = if (isTRUE(som_from_soc)) som_factor else NULL
      )
    )
  )
}
