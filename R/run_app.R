# =============================================================================
# esmc — Equivalent Soil Mass
# run_esmc_app(): launch the SOC Sampling Design Simulator
# =============================================================================


#' Launch the SOC Sampling Design Simulator
#'
#' Opens an interactive Shiny application for exploring how bulk density
#' changes affect Type I error rates, bias, and statistical power when
#' comparing soil organic carbon (SOC) stocks between treatments using
#' fixed-depth (FD) and equivalent soil mass (ESM) approaches.
#'
#' The app allows users to:
#' \itemize{
#'   \item Set up H0 (no true SOC change) or H1 (true SOC change) scenarios
#'   \item Control bulk density change, SOC stratification, and noise parameters
#'   \item Compare FD vs ESM performance via Monte Carlo simulation
#'   \item Compute power curves as a function of sample size
#' }
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}
#'   (e.g. \code{port}, \code{launch.browser}).
#'
#' @return No return value. Called for its side effect of launching the app.
#'
#' @examples
#' \dontrun{
#' run_esmc_app()
#' }
#'
#' @export
run_esmc_app <- function(...) {

  # Check required packages
  required_pkgs <- c("shiny", "ggplot2", "future", "future.apply")
  missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                          logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required to run the app but are not installed:\n  ",
      paste(missing_pkgs, collapse = ", "), "\n",
      "Install them with:\n  ",
      'install.packages(c("', paste(missing_pkgs, collapse = '", "'), '"))',
      call. = FALSE
    )
  }

  app_dir <- system.file("app", package = "esmc")

  if (app_dir == "") {
    stop(
      "Could not find the app directory in the esmc package.\n",
      "Try re-installing the package with: remotes::install_github('sebavillarino/esmc')",
      call. = FALSE
    )
  }

  shiny::runApp(app_dir, ...)
}
