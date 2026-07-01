# =============================================================================
# esmc — Equivalent Soil Mass
# Visualization: esm_plot_spline()
#
# Plots the cubic spline (or linear) interpolation curve fitted to each
# core's C(M) data — cumulative SOC stock as a function of cumulative
# mineral soil mass — alongside the observed data points.
#
# When the reference mass (target) requires extrapolation beyond the
# observed data range, the curve is extended to the target and the
# extrapolated segment is drawn with a distinct linetype, color, and
# a shaded background region — making the extrapolation visually obvious.
# =============================================================================

utils::globalVariables(c(
  "xmin", "xmax",
  "cum_mineral_mass_kg_m2", "cum_stock_kg_m2",
  "target", "target_type"
))

#' Plot spline (or linear) interpolation fit for ESM profiles
#'
#' For each soil core, plots cumulative SOC stock vs cumulative mineral soil
#' mass: observed nodes as points and the fitted interpolation curve as a
#' line. When the reference (target) mineral soil mass requires extrapolation
#' beyond the observed data range, the curve is extended to the target and the
#' extrapolated segment is drawn with a distinct style (dashed line +
#' shaded background), making it immediately clear why ESM estimates can be
#' unreliable at deep cumulative depths (e.g. 0–100 cm).
#'
#' @param esm_out List. Output from \code{\link{esm_compare}}, which must
#'   contain elements \code{$profile_mass} and \code{$diagnostics}.
#'   This requires \code{esm_compare()} to be called with
#'   \code{output = "both"} (default) or \code{output = "cumulative"}.
#'   It will not work with \code{output = "layers"} because the plot
#'   requires the cumulative C(M) values stored in \code{$diagnostics}.
#' @param profile_by Character vector. Same \code{profile_by} used in
#'   \code{esm_compare()} — identifies individual cores/profiles.
#' @param interp Character. Interpolation method: \code{"spline"} (default)
#'   or \code{"linear"}. Should match the method used in \code{esm_compare()}.
#' @param n_grid Integer. Number of points used to draw each curve segment.
#'   Default 200.
#' @param x_lab Character. X-axis label.
#' @param y_lab Character. Y-axis label.
#' @param facet_ncol Integer or NULL. Number of columns in the facet grid.
#'   NULL lets ggplot2 choose automatically.
#' @param point_size Numeric. Size of observed data points. Default 2.5.
#' @param line_width Numeric. Width of the fitted curve. Default 0.8.
#' @param show_target Logical. If TRUE (default), draw a vertical line at
#'   the reference mineral soil mass.
#' @param color_interp Character. Color for interpolation region elements
#'   (curve, target line). Default \code{"#2166ac"} (blue).
#' @param color_extrap Character. Color for extrapolation region elements
#'   (extended curve, target line, shading). Default \code{"#d73027"} (red).
#' @param extrap_alpha Numeric in the range 0 to 1. Transparency of the
#' extrapolation shaded region. Default 0.10.
#' @param free_scales Logical. If TRUE, each panel gets its own axis scales.
#'   Default FALSE (shared scales for comparability).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' result <- esm_compare(
#'   data = my_data, rvar = "soc",
#'   profile_by = c("site", "yr", "trt", "rep"),
#'   output = "both",
#'   reference_mode = "min", interp = "spline"
#' )
#' esm_plot_spline(result, profile_by = c("site", "yr", "trt", "rep"))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_vline geom_rect
#'   facet_wrap labs theme_bw theme scale_color_manual scale_linetype_manual
#' @export

esm_plot_spline <- function(
    esm_out,
    profile_by,
    interp        = c("spline", "linear"),
    n_grid        = 200,
    x_lab         = "Cumulative mineral soil mass (kg m\u207b\u00b2)",
    y_lab         = "Cumulative SOC stock (kg C m\u207b\u00b2)",
    facet_ncol    = NULL,
    point_size    = 2.5,
    line_width    = 0.8,
    show_target   = TRUE,
    color_interp  = "#2166ac",
    color_extrap  = "#d73027",
    extrap_alpha  = 0.10,
    free_scales   = FALSE
) {
  interp <- match.arg(interp)

  # ---- Validate input -------------------------------------------------------
  if (!is.list(esm_out) || !"profile_mass" %in% names(esm_out)) {
    stop(
      "esm_out must be the list returned by esm_compare(), ",
      "containing at least $profile_mass."
    )
  }
  if (!"diagnostics" %in% names(esm_out)) {
    stop(
      "esm_out does not contain $diagnostics, which is required for esm_plot_spline().\n",
      "Re-run esm_compare() with output = \"both\" or output = \"cumulative\" to include diagnostics."
    )
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install it with: install.packages('ggplot2')")
  }

  pm   <- esm_out$profile_mass
  diag <- esm_out$diagnostics

  miss_pm <- setdiff(profile_by, names(pm))
  if (length(miss_pm) > 0) {
    stop("profile_by columns not found in $profile_mass: ", paste(miss_pm, collapse = ", "))
  }

  # ---- Core labels (facet strips) -------------------------------------------
  pm$.__core_label__   <- .esm_core_label(pm,   profile_by)
  diag$.__core_label__ <- .esm_core_label(diag, profile_by)

  # ---- Build target info per core -------------------------------------------
  # We need target + max_obs to know how far to extend the curve
  has_target_cols <- all(c("target_mass_kg_m2", "max_cum_mass_kg_m2") %in% names(diag))

  target_df <- NULL
  if (show_target && nrow(diag) > 0 && has_target_cols) {
    target_df <- data.frame(
      .__core_label__ = diag$.__core_label__,
      target          = diag$target_mass_kg_m2,
      max_obs         = diag$max_cum_mass_kg_m2,
      stringsAsFactors = FALSE
    )
    target_df <- target_df[is.finite(target_df$target), , drop = FALSE]
    target_df$extrap      <- target_df$target > target_df$max_obs
    target_df$target_type <- ifelse(target_df$extrap, "Extrapolation", "Interpolation")
  }

  # ---- Build fitted curves (with optional extrapolation extension) -----------
  core_labels <- unique(pm$.__core_label__)

  # Two data frames: one for the interpolation segment, one for extrapolation
  curve_interp_list <- list()
  curve_extrap_list <- list()

  for (lbl in core_labels) {
    sub <- pm[pm$.__core_label__ == lbl, , drop = FALSE]
    sub <- sub[order(sub$cum_mineral_mass_kg_m2), , drop = FALSE]

    x <- sub$cum_mineral_mass_kg_m2
    y <- sub$cum_stock_kg_m2
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    npts <- length(x)
    if (npts < 2) next

    x_min <- min(x)
    x_max <- max(x)

    # How far to extend? Up to the target, if there is one for this core
    x_end <- x_max
    if (!is.null(target_df)) {
      td <- target_df[target_df$.__core_label__ == lbl & target_df$extrap, , drop = FALSE]
      if (nrow(td) > 0) x_end <- max(td$target, na.rm = TRUE)
    }

    # Build the spline/linear function over the full observed range
    fit_fn <- tryCatch({
      if (interp == "spline") {
        if (npts < 3) {
          warning("Core '", lbl, "' has < 3 points; using linear interpolation.")
          function(xout) stats::approx(x, y, xout = xout, method = "linear", rule = 2)$y
        } else {
          f <- stats::splinefun(x, y, method = "natural")
          f
        }
      } else {
        function(xout) stats::approx(x, y, xout = xout, method = "linear", rule = 2)$y
      }
    }, error = function(e) {
      warning("Could not fit curve for core '", lbl, "': ", conditionMessage(e))
      NULL
    })
    if (is.null(fit_fn)) next

    # --- Interpolation segment: x_min to x_max --------------------------------
    x_seq_i <- seq(x_min, x_max, length.out = n_grid)
    curve_interp_list[[lbl]] <- data.frame(
      .__core_label__ = lbl,
      x = x_seq_i,
      y = fit_fn(x_seq_i),
      stringsAsFactors = FALSE
    )

    # --- Extrapolation segment: x_max to x_end (only if needed) --------------
    if (x_end > x_max) {
      x_seq_e <- seq(x_max, x_end, length.out = max(10, round(n_grid * (x_end - x_max) / (x_max - x_min + 1e-9))))
      curve_extrap_list[[lbl]] <- data.frame(
        .__core_label__ = lbl,
        x = x_seq_e,
        y = fit_fn(x_seq_e),
        stringsAsFactors = FALSE
      )
    }
  }

  curve_interp_df <- do.call(rbind, curve_interp_list[lengths(curve_interp_list) > 0])
  curve_extrap_df <- do.call(rbind, curve_extrap_list[lengths(curve_extrap_list) > 0])

  # ---- Observed points (exclude C(0)=0 anchor) ------------------------------
  obs_df <- pm[is.finite(pm$cum_mineral_mass_kg_m2) & is.finite(pm$cum_stock_kg_m2) &
                 pm$cum_mineral_mass_kg_m2 > 0, , drop = FALSE]

  # ---- Shading rectangle for extrapolation region per panel -----------------
  shade_df <- NULL
  if (!is.null(target_df)) {
    extrap_cores <- target_df[target_df$extrap, , drop = FALSE]
    if (nrow(extrap_cores) > 0) {
      # xmin of shade = max observed mass for that core
      shade_df <- data.frame(
        .__core_label__ = extrap_cores$.__core_label__,
        xmin            = extrap_cores$max_obs,
        xmax            = extrap_cores$target,
        stringsAsFactors = FALSE
      )
    }
  }

  # ---- Assemble ggplot -------------------------------------------------------
  scales_arg <- if (free_scales) "free" else "fixed"

  p <- ggplot2::ggplot()

  # Extrapolation shading (drawn first, behind everything)
  if (!is.null(shade_df) && nrow(shade_df) > 0) {
    p <- p + ggplot2::geom_rect(
      data = shade_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill  = color_extrap,
      alpha = extrap_alpha,
      inherit.aes = FALSE
    )
  }

  # Interpolation curve segment (solid)
  if (!is.null(curve_interp_df) && nrow(curve_interp_df) > 0) {
    p <- p + ggplot2::geom_line(
      data = curve_interp_df,
      ggplot2::aes(x = x, y = y),
      color     = "grey30",
      linewidth = line_width,
      linetype  = "solid"
    )
  }

  # Extrapolation curve segment (dashed, red)
  if (!is.null(curve_extrap_df) && nrow(curve_extrap_df) > 0) {
    p <- p + ggplot2::geom_line(
      data = curve_extrap_df,
      ggplot2::aes(x = x, y = y),
      color     = color_extrap,
      linewidth = line_width,
      linetype  = "dashed"
    )
  }

  # Observed nodes
  p <- p + ggplot2::geom_point(
    data = obs_df,
    ggplot2::aes(x = cum_mineral_mass_kg_m2, y = cum_stock_kg_m2),
    size   = point_size,
    color  = "grey20",
    shape  = 21,
    fill   = "white",
    stroke = 1
  )

  # Target vertical lines (colored by interp / extrap)
  if (!is.null(target_df) && nrow(target_df) > 0) {
    p <- p +
      ggplot2::geom_vline(
        data = target_df,
        ggplot2::aes(xintercept = target, color = target_type),
        linetype  = "dashed",
        linewidth = 0.7
      ) +
      ggplot2::scale_color_manual(
        name   = "Target mass",
        values = c("Interpolation" = color_interp, "Extrapolation" = color_extrap),
        labels = c(
          "Interpolation" = "Interpolation\n(target within data range)",
          "Extrapolation" = "Extrapolation\n(target beyond data range)"
        )
      )
  }

  p <- p +
    ggplot2::facet_wrap(
      ~ .__core_label__,
      ncol   = facet_ncol,
      scales = scales_arg
    ) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey60"),
      strip.text       = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(size = 9),
      legend.text      = ggplot2::element_text(size = 8)
    )

  p
}


# =============================================================================
# Internal helper: build a readable core label from profile_by columns
# =============================================================================
.esm_core_label <- function(df, profile_by) {
  cols_present <- intersect(profile_by, names(df))
  if (length(cols_present) == 0) return(rep("core", nrow(df)))
  apply(df[cols_present], 1, paste, collapse = " | ")
}
