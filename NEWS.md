# esmc 0.4.0

* **New default interpolation:** `esm_compare()` and `esm_plot_spline()` now use
  `interp = "hyman"` (cubic spline with Hyman monotonicity filter) by default.
  The choice is based on a simulation benchmark of linear, natural spline,
  Hyman and Fritsch-Carlson monotone interpolation across fine and coarse
  sampling designs and contrasting SOC stratification. With the end slope
  described below, Hyman had the lowest error in the upper layers of coarse
  designs and was comparable to the best method elsewhere. Results obtained
  with earlier versions can be reproduced with `interp = "spline"`.
* `interp = "hyman"`: the slope at the deepest node is the secant of the last
  interval (the concentration per unit mineral mass of the deepest sampled
  layer) instead of the end condition of `stats::splinefun()`, and the curve is
  extrapolated linearly with that slope. With few depth intervals the end
  condition of `stats::splinefun()` biased the last interval and, above all,
  the extrapolation (simulation benchmark: RMSE of the 0-100 cm stock change
  0.51 vs 0.20 kg C m-2 with coarse sampling and a 20% BD decrease). Inside all
  intervals but the last one, results equal `stats::spline(method = "hyman")`
  (used by von Haden et al. 2020, which extrapolates with the cubic of the end
  interval).
* `esm_compare()` now checks that each profile is a contiguous sequence of
  layers starting at 0 cm. Overlapping or duplicated layers give an error;
  profiles not starting at 0 cm are dropped, and profiles with a gap (e.g.
  rows removed because of missing values) are truncated at the first gap,
  both with a warning. Previously these cases silently produced wrong
  cumulative masses.
* `esm_compare()`, `reference_mode = "custom"`: `custom_mass` now works as a
  named vector (names = lower boundaries) or a data.frame (`lower`,
  `target_mass_kg_m2` and, optionally, the `reference_by` columns); both
  previously failed in argument validation. With per-boundary targets, layers
  are defined by the supplied boundaries and the cumulative stock is evaluated
  at the deepest one. A single `custom_mass` defines only the cumulative
  reference mass: `output = "layers"` now gives an error and `output = "both"`
  returns only the cumulative stock with a warning (previously every layer
  below the first one had zero stock).
* `esm_compare()`: `reference_by` must be a subset of `profile_by`
  (informative error instead of "undefined columns selected");
  `reference_mode = "ref_trt"` no longer requires `trt` in `profile_by`, and
  checks that `ref_trt` exists in every stratum.
* `esm_uncert_depth()`:
  * the uncertainty of cumulative stocks is computed from the per-core
    cumulative stocks, so the covariance among layers of the same core is
    included (previously layer variances were summed as if independent, which
    underestimated it);
  * new argument `uncertainty = c("sd", "se")`; both are returned;
  * `method = "mc"`, `n_mc` and `seed` are deprecated (the global random seed
    is no longer modified); SOM arguments are ignored with a warning, since
    fixed-depth stocks use total soil mass, and the `som_source`,
    `som_factor` and SOM flag columns are no longer returned;
  * cores are checked as in `esm_compare()`, and cores of a stratum must share
    the same depth increments.
* `esm_compare()`: `$cumulative` and `$diagnostics` no longer include the
  internal `.__stratum__` column, and their row names are reset.
* `esm_plot_spline()`: corrected an outdated argument name in the example.
* Vignette and README updated: per-boundary `custom_mass`, profile geometry
  checks, the Hyman end slope, and `esm_uncert_depth()`.
* Shiny app: Hyman interpolation by default (all methods selectable); the
  simulated BD change now declines exponentially with depth (e-folding depth
  29 cm, calibrated with GSOCS-LULCC) and can be negative; the lower bound on
  observed BD was reduced from 0.6 to 0.05 g cm-3.

# esmc 0.3.2

* `esm_compare()`: two new interpolation methods for the cumulative curve C(M):
  * `interp = "monotone"`: monotone cubic Hermite spline (Fritsch-Carlson,
    `stats::splinefun(method = "monoH.FC")`). Locally defined and never
    decreasing, so layer stocks cannot be negative. Extrapolates linearly with
    the slope of the last segment (concentration of the deepest sampled layer).
  * `interp = "hyman"`: cubic spline with Hyman monotonicity filter, the method
    used in the ESM script of von Haden et al. (2020). Extrapolates with the
    cubic polynomial of the end interval.
* The default (`interp = "spline"`, natural cubic spline) is unchanged.
* `esm_plot_spline()`: supports the new methods and now uses the same
  interpolating function as `esm_compare()`. With `interp = "linear"` the
  extrapolated segment now follows the slope of the last segment (as in
  `esm_compare()`) instead of a constant value.

# esmc 0.3.1

* `esm_compare()`: SOC (or other property) stocks are now computed as
  concentration times **total** soil mass. Previously they were computed with
  mineral soil mass, which underestimated stocks by roughly the SOM fraction.
  Mineral soil mass is still used as the coordinate for mass standardization.
* `esm_compare()`: in the layers branch, profiles are now evaluated on a
  non-overlapping partition of all sampled boundaries within each stratum.
  Previously, when depth increments differed among profiles (e.g. between
  years), overlapping layers were returned and some boundaries appeared more
  than once per profile, which could duplicate observations in downstream
  tests.
* `esm_compare()`: new warning when reference masses do not increase
  monotonically with depth within a stratum.
* `esm_compare()`: removed spurious "row names were found from a short
  variable" warnings.
* `esm_uncert_depth()`: stocks are now computed with total soil mass.
* Shiny app: fixed a lazy-evaluation bug that made all simulated cores share
  the bulk density multiplier of the last core; default to unpaired cores and
  a between-core BD CV of 0.10.

# esmc 0.3.0

* Unified `esm_compare()` API replacing `esm_compare_layers()`.
* Extrapolation beyond the observed mass range is performed and flagged
  (`extrapolated = TRUE`), with a warning.
