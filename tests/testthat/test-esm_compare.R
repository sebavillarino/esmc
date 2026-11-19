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

# --- Interpolation methods (v0.3.2) -------------------------------------------
one_profile_data <- function() {
  data.frame(trt = c("a", "a", "a", "a", "b", "b", "b", "b"), rep = 1,
             upper = c(0, 10, 30, 60, 0, 10, 30, 60),
             lower = c(10, 30, 60, 100, 10, 30, 60, 100),
             bd  = c(1.0, 1.2, 1.3, 1.35, 1.2, 1.3, 1.35, 1.4),
             soc = c(3.0, 1.5, 0.8, 0.4, 3.0, 1.5, 0.8, 0.4))
}

test_that("all interpolation methods reproduce the observed nodes", {
  d <- one_profile_data()
  for (m in c("spline", "monotone", "hyman", "linear")) {
    r  <- esm_compare(d[d$trt == "a", ], soil_prop = "soc",
                      profile_by = c("trt", "rep"), output = "layers",
                      reference_mode = "min", interp = m)
    iv <- r$interval
    expect_equal(r$layer$cum_lower_kg_m2, iv$cum_stock_kg_m2, tolerance = 1e-10)
  }
})

test_that("interp = 'hyman' matches stats::spline(method = 'hyman')", {
  d  <- one_profile_data()
  r  <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                    output = "layers", reference_mode = "min", interp = "hyman")
  pm <- r$profile_mass[r$profile_mass$trt == "b", ]
  tg <- r$targets
  ref <- stats::spline(pm$cum_mineral_mass_kg_m2, pm$cum_stock_kg_m2,
                       xout = tg$target_mass_kg_m2, method = "hyman")$y
  got <- r$layer$cum_lower_kg_m2[r$layer$trt == "b"]
  expect_equal(got, ref[match(r$layer$lower[r$layer$trt == "b"], tg$lower)],
               tolerance = 1e-10)
})

test_that("interp = 'monotone' gives non-negative layer stocks", {
  d <- one_profile_data()
  r <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                   output = "layers", reference_mode = "min", interp = "monotone")
  expect_true(all(r$layer$layer_stock_kg_m2 >= 0))
})

test_that("interp = 'monotone' extrapolates with the last-layer concentration", {
  d <- one_profile_data()[1:4, ]
  expect_warning(
    r <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                     output = "cumulative", reference_mode = "custom",
                     custom_mass = 2000, interp = "monotone",
                     som_from_soc = FALSE),
    "extrapolation")
  iv   <- r$interval
  last <- iv[which.max(iv$lower), ]
  conc <- last$stock_kg_m2 / last$mineral_mass_kg_m2
  expected <- last$cum_stock_kg_m2 + conc * (2000 - last$cum_mineral_mass_kg_m2)
  expect_true(r$cumulative$extrapolated)
  expect_equal(r$cumulative$cum_stock_kg_m2, expected, tolerance = 1e-10)
})

# --- Default interpolation and Hyman extrapolation (v0.4.0) --------------------
test_that("the default interpolation is 'hyman'", {
  expect_identical(eval(formals(esm_compare)$interp)[1], "hyman")
})

test_that("interp = 'hyman' extrapolates linearly with the end slope", {
  d <- one_profile_data()[1:4, ]
  r <- suppressWarnings(
    esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                output = "cumulative", reference_mode = "custom",
                custom_mass = 2000, interp = "hyman", som_from_soc = FALSE))
  pm <- r$profile_mass
  f  <- stats::splinefun(pm$cum_mineral_mass_kg_m2, pm$cum_stock_kg_m2,
                         method = "hyman")
  xn <- max(pm$cum_mineral_mass_kg_m2)
  expected <- f(xn) + f(xn, deriv = 1) * (2000 - xn)
  expect_equal(r$cumulative$cum_stock_kg_m2, expected, tolerance = 1e-10)
})
