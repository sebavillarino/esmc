# Two groups (years) with different depth increments
make_data <- function() {
  lay <- list(
    "2002" = data.frame(upper = c(0, 15, 30, 50), lower = c(15, 30, 50, 100)),
    "2014" = data.frame(upper = c(0, 15, 30, 60), lower = c(15, 30, 60, 100))
  )
  out <- list()
  for (yr in names(lay)) for (trt in c("a", "b")) for (rep in 1:2) {
    l  <- lay[[yr]]
    bd <- 1.1 + 0.003 * l$lower + 0.02 * rep + ifelse(trt == "b", 0.05, 0)
    soc <- 2.5 * exp(-l$lower / 60) + 0.1 * rep
    out[[length(out) + 1]] <- data.frame(yr = yr, trt = trt, rep = rep, l,
                                         bd = bd, soc = soc)
  }
  do.call(rbind, out)
}

test_that("each boundary appears once per profile when increments differ", {
  res <- suppressWarnings(esm_compare(
    make_data(), soil_prop = "soc", profile_by = c("yr", "trt", "rep"),
    output = "layers", reference_mode = "min", interp = "spline",
    som_from_soc = TRUE))
  L <- res$layer
  n <- table(interaction(L$yr, L$trt, L$rep, L$lower, drop = TRUE))
  expect_true(all(n == 1))
  expect_setequal(unique(L$lower), c(15, 30, 50, 60, 100))
})

test_that("layer stocks add up to the cumulative stock", {
  res <- suppressWarnings(esm_compare(
    make_data(), soil_prop = "soc", profile_by = c("yr", "trt", "rep"),
    output = "layers", reference_mode = "min", interp = "spline",
    som_from_soc = TRUE))
  L <- res$layer
  for (p in split(L, interaction(L$yr, L$trt, L$rep, drop = TRUE))) {
    expect_equal(sum(p$layer_stock_kg_m2),
                 p$cum_lower_kg_m2[p$lower == max(p$lower)])
  }
})

test_that("stocks use total soil mass, not mineral mass", {
  d   <- make_data()
  res <- esm_compare(d, soil_prop = "soc", profile_by = c("yr", "trt", "rep"),
                     output = "cumulative", reference_mode = "min",
                     som_from_soc = TRUE)
  iv  <- res$interval
  expect_equal(iv$stock_kg_m2,
               d$soc[match(paste(iv$yr, iv$trt, iv$rep, iv$upper),
                           paste(d$yr, d$trt, d$rep, d$upper))] / 100 *
                 iv$soil_mass_kg_m2)
})

test_that("results are unchanged when all profiles share increments", {
  d <- make_data()
  d <- d[d$yr == "2002", ]
  res <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                     output = "layers", reference_mode = "min",
                     som_from_soc = TRUE)
  expect_setequal(unique(res$layer$lower), c(15, 30, 50, 100))
})
