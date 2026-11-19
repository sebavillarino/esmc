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

test_that("interp = 'hyman' matches stats::spline(method = 'hyman') except in the last interval", {
  d  <- one_profile_data()
  r  <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                    output = "layers", reference_mode = "min", interp = "hyman")
  pm <- r$profile_mass[r$profile_mass$trt == "b", ]
  tg <- r$targets
  # targets inside all intervals but the last one
  inner <- tg$target_mass_kg_m2 < sort(pm$cum_mineral_mass_kg_m2, decreasing = TRUE)[2]
  ref <- stats::spline(pm$cum_mineral_mass_kg_m2, pm$cum_stock_kg_m2,
                       xout = tg$target_mass_kg_m2[inner], method = "hyman")$y
  got <- r$layer$cum_lower_kg_m2[r$layer$trt == "b"][match(tg$lower[inner],
                                   r$layer$lower[r$layer$trt == "b"])]
  expect_equal(got, ref, tolerance = 1e-10)
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

test_that("interp = 'hyman' extrapolates with the last-layer concentration", {
  d <- one_profile_data()[1:4, ]
  r <- suppressWarnings(
    esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                output = "cumulative", reference_mode = "custom",
                custom_mass = 2000, interp = "hyman", som_from_soc = FALSE))
  iv   <- r$interval
  last <- iv[which.max(iv$lower), ]
  conc <- last$stock_kg_m2 / last$mineral_mass_kg_m2
  expected <- last$cum_stock_kg_m2 + conc * (2000 - last$cum_mineral_mass_kg_m2)
  expect_equal(r$cumulative$cum_stock_kg_m2, expected, tolerance = 1e-10)
})

test_that("interp = 'hyman' gives non-negative layer stocks", {
  d <- one_profile_data()
  r <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                   output = "layers", reference_mode = "min", interp = "hyman")
  expect_true(all(r$layer$layer_stock_kg_m2 >= 0))
})

# --- Profile geometry (v0.4.0) ------------------------------------------------
test_that("a gap from a missing value truncates the profile with a warning", {
  d <- one_profile_data()
  d$soc[d$trt == "b" & d$upper == 10] <- NA
  expect_warning(
    r <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                     output = "cumulative", reference_mode = "min",
                     interp = "linear"),
    "gap")
  expect_equal(max(r$interval$lower[r$interval$trt == "b"]), 10)
})

test_that("profiles not starting at 0 are dropped with a warning", {
  d <- one_profile_data()
  d <- d[!(d$trt == "b" & d$upper == 0), ]
  expect_warning(
    r <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                     output = "cumulative", reference_mode = "min",
                     interp = "linear"),
    "do not start at 0")
  expect_false("b" %in% r$cumulative$trt)
})

test_that("overlapping layers give an error", {
  d <- rbind(one_profile_data(),
             data.frame(trt = "b", rep = 1, upper = 0, lower = 30, bd = 1.2, soc = 2))
  expect_error(esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                           interp = "linear"), "Overlapping")
})

# --- custom_mass, reference_by and ref_trt validation (v0.4.0) ---------------
test_that("custom_mass as a named vector defines layers and the cumulative", {
  d  <- one_profile_data()
  cm <- c("30" = 300, "60" = 650)
  r  <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                    reference_mode = "custom", custom_mass = cm,
                    interp = "linear")
  expect_setequal(unique(r$layer$lower), c(30, 60))
  expect_equal(unique(r$layer$target_mass_lower_kg_m2[r$layer$lower == 60]), 650)
  expect_equal(unique(r$cumulative$target_mass_kg_m2), 650)
  for (t in c("a", "b"))
    expect_equal(sum(r$layer$layer_stock_kg_m2[r$layer$trt == t]),
                 r$cumulative$cum_stock_kg_m2[r$cumulative$trt == t])
})

test_that("custom_mass as a data.frame by stratum", {
  d <- rbind(cbind(site = "s1", one_profile_data()),
             cbind(site = "s2", one_profile_data()))
  cm <- data.frame(site = c("s1", "s1", "s2", "s2"), lower = c(30, 60, 30, 60),
                   target_mass_kg_m2 = c(300, 650, 280, 600))
  r <- esm_compare(d, soil_prop = "soc", profile_by = c("site", "trt", "rep"),
                   reference_by = "site", reference_mode = "custom",
                   custom_mass = cm, interp = "linear")
  expect_equal(sort(unique(r$cumulative$target_mass_kg_m2)), c(600, 650))
  expect_equal(nrow(r$layer), 2 * 2 * 2)
})

test_that("a single custom_mass gives only the cumulative stock", {
  d <- one_profile_data()
  expect_error(esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                           output = "layers", reference_mode = "custom",
                           custom_mass = 500), "one target mass per")
  expect_warning(r <- esm_compare(d, soil_prop = "soc",
                                  profile_by = c("trt", "rep"),
                                  reference_mode = "custom", custom_mass = 500,
                                  interp = "linear"),
                 "layer stocks are not computed")
  expect_null(r$layer)
  expect_equal(unique(r$cumulative$target_mass_kg_m2), 500)
})

test_that("invalid custom_mass gives informative errors", {
  d <- one_profile_data()
  f <- function(cm) esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                                reference_mode = "custom", custom_mass = cm,
                                output = "cumulative", interp = "linear")
  expect_error(f(c(300, 650)), "custom_mass must be one of")
  expect_error(f(c(a = 300)), "lower boundaries")
  expect_error(f(-5), "positive")
  expect_error(f(data.frame(lower = 30)), "missing columns")
})

test_that("reference_by must be a subset of profile_by", {
  d <- cbind(site = "s1", one_profile_data())
  d$core <- paste(d$trt, d$rep)
  expect_error(esm_compare(d, soil_prop = "soc", profile_by = "core",
                           reference_by = "site", interp = "linear"),
               "subset of profile_by")
})

test_that("ref_trt works when trt is not in profile_by", {
  d <- one_profile_data(); d$core <- paste(d$trt, d$rep)
  r1 <- esm_compare(d, soil_prop = "soc", profile_by = "core",
                    reference_mode = "ref_trt", ref_trt = "a", interp = "linear")
  r2 <- esm_compare(d, soil_prop = "soc", profile_by = c("trt", "rep"),
                    reference_mode = "ref_trt", ref_trt = "a", interp = "linear")
  expect_equal(r1$cumulative$cum_stock_kg_m2, r2$cumulative$cum_stock_kg_m2)
})

test_that("ref_trt validation", {
  d <- rbind(cbind(site = "s1", one_profile_data()),
             cbind(site = "s2", one_profile_data()[5:8, ]))
  expect_error(esm_compare(d, soil_prop = "soc", profile_by = c("site", "trt", "rep"),
                           reference_by = "site", reference_mode = "ref_trt",
                           ref_trt = "a", interp = "linear"), "missing in stratum")
  expect_error(esm_compare(d, soil_prop = "soc", profile_by = c("site", "trt", "rep"),
                           reference_mode = "ref_trt", ref_trt = "zz",
                           interp = "linear"), "not found")
})
