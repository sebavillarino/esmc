make_cores <- function() {
  set.seed(42)
  do.call(rbind, lapply(c("a", "b"), function(t) do.call(rbind, lapply(1:4, function(r) {
    f <- exp(rnorm(1, 0, 0.15))          # core effect: correlated among layers
    data.frame(trt = t, rep = r, upper = c(0, 10, 30, 50), lower = c(10, 30, 50, 100),
               bd = c(1.1, 1.25, 1.35, 1.4), soc = c(2.5, 1.4, 0.8, 0.5) * f)
  }))))
}

test_that("stocks use total soil mass", {
  d <- make_cores()
  u <- esm_uncert_depth(d, rvar = "soc", core_by = "trt")
  a <- d[d$trt == "a" & d$upper == 0, ]
  expect_equal(u$depth$stock_layer_mu_kg_m2[u$depth$trt == "a" & u$depth$upper == 0],
               mean(a$soc / 100 * a$bd * 10 * 10))
})

test_that("cumulative SD is the SD of per-core cumulative stocks", {
  d <- make_cores()
  u <- esm_uncert_depth(d, rvar = "soc", core_by = "trt")
  a <- d[d$trt == "a", ]
  a$st <- a$soc / 100 * a$bd * (a$lower - a$upper) * 10
  cum100 <- tapply(a$st, a$rep, sum)
  got <- u$depth$cum_stock_sd_kg_m2[u$depth$trt == "a" & u$depth$lower == 100]
  expect_equal(got, sd(cum100))
  # with correlated cores it exceeds the independent-layer sum
  indep <- sqrt(sum(tapply(a$st, a$upper, var)))
  expect_gt(got, indep)
})

test_that("se = sd / sqrt(n)", {
  u <- esm_uncert_depth(make_cores(), rvar = "soc", core_by = "trt")
  expect_equal(u$depth$cum_stock_se_kg_m2,
               u$depth$cum_stock_sd_kg_m2 / sqrt(u$depth$n_rep))
})

test_that("the global RNG state is not modified", {
  d <- make_cores()
  set.seed(1); x1 <- runif(1)
  set.seed(1); invisible(suppressWarnings(
    esm_uncert_depth(d, rvar = "soc", core_by = "trt", method = "mc", seed = 99)))
  x2 <- runif(1)
  expect_equal(x1, x2)
})

test_that("a core truncated by a missing value is excluded from deeper layers", {
  d <- make_cores()
  d$soc[d$trt == "a" & d$rep == 1 & d$upper == 30] <- NA
  expect_warning(u <- esm_uncert_depth(d, rvar = "soc", core_by = "trt"), "gap")
  a <- u$depth[u$depth$trt == "a", ]
  expect_equal(a$n_rep, c(4, 4, 3, 3))
})

test_that("cores with different increments give an error", {
  d <- make_cores()
  d2 <- d[!(d$trt == "a" & d$rep == 1), ]
  d2 <- rbind(d2, data.frame(trt = "a", rep = 1, upper = c(0, 15, 30, 50),
                             lower = c(15, 30, 50, 100), bd = 1.2, soc = 1))
  expect_error(esm_uncert_depth(d2, rvar = "soc", core_by = "trt"),
               "same depth")
})

test_that("SOM arguments warn and do not change results", {
  d  <- make_cores()
  u1 <- esm_uncert_depth(d, rvar = "soc", core_by = "trt")
  expect_warning(u2 <- esm_uncert_depth(d, rvar = "soc", core_by = "trt",
                                        som_from_soc = TRUE), "ignored")
  expect_equal(u1$depth, u2$depth)
})
