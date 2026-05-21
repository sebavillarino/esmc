# esmc — Equivalent Soil Mass Comparisons

<!-- badges: start -->
<!-- badges: end -->

## The problem

Comparing soil organic carbon (SOC) stocks — or any other soil property that varies with depth — between treatments or over time is a common goal in soil science. The standard approach is to sample soil at fixed depths (e.g., 0–10, 10–30, 30–50 cm) and compare stocks across those intervals. However, the mass of soil within a fixed-depth interval depends on bulk density (BD), which varies with management, texture, compaction, and organic matter content. When BD differs between treatments or changes over time, samples collected at the same depth contain different amounts of soil mass. Comparing SOC stocks across unequal soil masses can produce artifacts — creating apparent differences where none exist, or masking real ones.

The **Equivalent Soil Mass (ESM) framework** corrects for this by standardizing SOC stocks to a common mineral soil mass reference rather than a common depth. In practice, ESM is applied as a post-hoc correction of fixed-depth samples: a continuous curve of cumulative SOC stock as a function of cumulative mineral soil mass — C(M) — is fitted for each soil profile, and stocks are evaluated at a common reference mass via interpolation.

Despite its wide use, key methodological choices in ESM implementation remain unresolved: which interpolation method to use, how to select the reference mass, and how sampling design affects accuracy. `esmc` provides a unified, reproducible framework for applying ESM based on the recommendations of Villarino et al. (2026).

## Installation

```r
# install.packages("remotes")
remotes::install_github("sebavillarino/esmc")
```

## Main functions

| Function | Description |
|---|---|
| `esm_compare()` | Core function — computes ESM-adjusted stocks by layer or cumulatively |
| `esm_plot_spline()` | Diagnostic plot of the C(M) spline fit per profile |
| `run_esmc_app()` | Launch the interactive Shiny simulator |

## Basic usage

```r
library(esmc)

# Compute ESM-adjusted SOC stocks by depth layer
result <- esm_compare(
  data           = my_data,
  soil_prop      = "soc",
  profile_by     = c("site", "year", "trt", "rep"),
  output         = "layers",
  reference_mode = "min",
  interp         = "spline",
  soil_prop_unit = "pct",
  som_from_soc   = TRUE
)

# ESM stocks per layer
result$layer

# Diagnostic plot: C(M) spline fit per profile
esm_plot_spline(result, profile_by = c("site", "year", "trt", "rep"))
```

### Input data format

`esm_compare()` expects a data frame with one row per depth interval per profile, with at least these columns:

| Column | Description |
|---|---|
| `upper` | Upper boundary of depth interval (cm) |
| `lower` | Lower boundary of depth interval (cm) |
| `bd` | Bulk density (g cm⁻³) |
| `trt` | Treatment identifier |
| `soc` | SOC concentration (or any soil property) |
| + grouping columns | e.g. `site`, `year`, `rep` |

### Key arguments

**`output`**: what to return.
- `"layers"` — ESM stock per depth layer (default, most common)
- `"cumulative"` — cumulative ESM stock at a single reference mass
- `"both"` — both of the above

**`reference_mode`**: how to define the reference mineral soil mass.
- `"min"` — minimum observed mass across profiles (recommended; avoids extrapolation)
- `"ref_trt"` — mean mass of a reference treatment
- `"custom"` — user-supplied value

**`interp`**: interpolation method for the C(M) curve.
- `"spline"` — cubic spline (recommended when ≥ 3 depth intervals)
- `"linear"` — recommended when only 2 depth intervals are available
- `"step"` — constant within interval

### Interactive simulator

```r
run_esmc_app()
```

Launches a Shiny app for exploring how bulk density changes affect Type I error, bias, and statistical power when comparing SOC stocks using fixed-depth vs ESM approaches. Useful for study design and for understanding when ESM corrections matter most.

> Requires: `shiny`, `ggplot2`, `future`, `future.apply`

## Citation

If you use `esmc` in your work, please cite:

> Villarino, S.H., Castellano, M., and Miguez, F.E. (2026). Towards a Unified Equivalent Soil Mass Framework. *In review*.

## License

MIT
