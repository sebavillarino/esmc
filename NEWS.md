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
