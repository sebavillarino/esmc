# =============================================================================
# esmc Shiny App — SOC Sampling Design Simulator
# Version: 0.3.0 (updated to esm_compare() v0.3.0 API)
#
# Dependencies (all in Suggests in DESCRIPTION):
#   shiny, ggplot2, future, future.apply
# =============================================================================

library(shiny)
library(ggplot2)
library(future)
library(future.apply)
library(esmc)

future::plan(multisession, workers = max(1, future::availableCores() - 1))

# -----------------------------------------------------------------------------
# Treatment labels
# -----------------------------------------------------------------------------
recode_trt <- function(x) {
  ifelse(x == "ref", "Reference", ifelse(x == "t1", "Treatment", as.character(x)))
}

fix_trt <- function(x) {
  factor(recode_trt(x), levels = c("Reference", "Treatment"))
}

# -----------------------------------------------------------------------------
# Robust treatment test (base R, no broom/dplyr)
# -----------------------------------------------------------------------------
safe_lm_trt <- function(dat, y = "layer_stock_kg_m2") {
  dat$trt <- fix_trt(dat$trt)
  tab <- table(dat$trt)

  empty_row <- data.frame(
    term      = "trtTreatment",
    estimate  = NA_real_, std.error = NA_real_,
    statistic = NA_real_, p.value   = NA_real_,
    conf.low  = NA_real_, conf.high = NA_real_,
    stringsAsFactors = FALSE
  )

  if (length(tab) < 2 || any(tab == 0)) return(empty_row)

  m           <- lm(stats::as.formula(paste(y, "~ trt")), data = dat)
  coef_summ   <- summary(m)$coefficients
  term_target <- paste0("trt", levels(dat$trt)[2])

  if (!term_target %in% rownames(coef_summ)) return(empty_row)

  ci <- tryCatch(confint(m), error = function(e) NULL)
  conf.low  <- if (!is.null(ci) && term_target %in% rownames(ci)) ci[term_target, 1] else NA_real_
  conf.high <- if (!is.null(ci) && term_target %in% rownames(ci)) ci[term_target, 2] else NA_real_

  data.frame(
    term      = term_target,
    estimate  = coef_summ[term_target, "Estimate"],
    std.error = coef_summ[term_target, "Std. Error"],
    statistic = coef_summ[term_target, "t value"],
    p.value   = coef_summ[term_target, "Pr(>|t|)"],
    conf.low  = conf.low,
    conf.high = conf.high,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Tooltip helper
# -----------------------------------------------------------------------------
with_info <- function(input_ui, content) {
  tagList(
    tags$div(
      style = "display:flex; align-items:flex-end; gap:8px; width:100%;",
      tags$div(style = "flex: 1 1 auto; min-width: 0;", input_ui),
      tags$span(
        class = "tip-icon",
        `data-toggle` = "tooltip",
        `data-bs-toggle` = "tooltip",
        title = content,
        style = "cursor: help; color:#2C7FB8; font-weight:700; user-select:none; margin-bottom:10px;",
        "\u2139"
      )
    )
  )
}

# -----------------------------------------------------------------------------
# Freeze inputs into a config object
# -----------------------------------------------------------------------------
make_cfg <- function(input) {
  list(
    mode          = input$mode,
    n_rep         = input$n_rep,
    alpha         = input$alpha,
    n_cores       = input$n_cores,
    paired        = input$paired,
    design        = input$design,
    interval      = as.numeric(input$interval),
    soc_prof      = input$soc_prof,
    bd_pct        = input$bd_pct,
    soc_eff       = if (input$mode == "null") 0 else input$soc_eff,
    soc_eff_depth = if (input$mode == "null") "full" else input$soc_eff_depth,
    bd_core_cv    = input$bd_core_cv,
    bd_lab_sd     = input$bd_lab_sd,
    c_lab_cv      = input$c_lab_cv,
    esm_interp    = input$esm_interp,
    cores_min     = if (!is.null(input$cores_min)) input$cores_min else NA_integer_,
    cores_max     = if (!is.null(input$cores_max)) input$cores_max else NA_integer_,
    cores_step    = if (!is.null(input$cores_step)) input$cores_step else NA_integer_
  )
}

# -----------------------------------------------------------------------------
# Simulation core (base R — no dplyr)
# -----------------------------------------------------------------------------
simulate_once <- function(
    bd_increase_pct  = 0,
    soc_effect_pct   = 0,
    soc_effect_depth = c("full", "surface"),
    seed      = NULL,
    n_cores   = 3,
    paired    = TRUE,
    layers    = list(c(0,15), c(15,30), c(30,50), c(50,70), c(70,100)),
    soc_prof  = c("high", "mid", "low"),
    c0 = NULL, c_end = NULL, m0 = NULL,
    z0_cm      = 0.5,
    bd_core_cv = 0.07,
    bd_lab_sd  = 0.03,
    c_lab_cv   = 0.05,
    dz_mass    = 0.1,
    z_effect_max = 50,
    w_transition = 6
) {
  if (!is.null(seed)) set.seed(seed)
  soc_prof         <- match.arg(soc_prof)
  soc_effect_depth <- match.arg(soc_effect_depth)

  bd_wiese <- function(z_cm, z0_cm = 0.5) {
    z_m <- (z_cm + z0_cm) / 100
    0.1140 * log(z_m) + 1.23
  }

  r_lognorm_mean1 <- function(cv) {
    sigma <- sqrt(log(1 + cv^2)); mu <- -0.5 * sigma^2
    exp(rnorm(1, mean = mu, sd = sigma))
  }

  add_lab_noise <- function(x, cv) {
    sigma <- sqrt(log(1 + cv^2)); mu <- -0.5 * sigma^2
    x * exp(rnorm(length(x), mu, sigma))
  }

  mk_mprof <- function(bd_fun, z_max = 100, dz = 0.1) {
    z <- seq(0, z_max, by = dz)
    bd <- bd_fun(z); dm <- bd * dz
    m <- c(0, cumsum(dm[-1]))
    data.frame(z = z, bd = bd, dm = dm, m = m)
  }

  c_mass_fun <- function(m, c0 = 3, c_end = 0.5, m_end, m0 = 1) {
    s <- log(c_end / c0) / log((m_end + m0) / m0)
    c0 * ((m + m0) / m0)^s
  }

  depth_multiplier_logistic <- function(z, bd_increase_pct, z_effect_max = 50, w = 6) {
    if (bd_increase_pct <= 0) return(rep(1, length(z)))
    inc <- bd_increase_pct / 100; z0 <- z_effect_max
    g <- plogis((z0 - z) / w); g0 <- plogis(0); gS <- plogis(z0 / w)
    f <- pmin(pmax((g - g0) / (gS - g0), 0), 1)
    1 + inc * f
  }

  make_bd_fun <- function(trt_name, core_mult, bd_increase_pct,
                           z0_cm = 0.5, z_effect_max = 50, w_transition = 6) {
    function(z) {
      bd0 <- bd_wiese(z, z0_cm = z0_cm) * core_mult
      if (trt_name == "t1")
        bd0 <- bd0 * depth_multiplier_logistic(z, bd_increase_pct,
                                                z_effect_max = z_effect_max,
                                                w = w_transition)
      bd0
    }
  }

  mk_bd_true_layers_from_fun <- function(bd_fun, layers, trt_name, core_id, core_mult) {
    res <- do.call(rbind, lapply(layers, function(z) {
      d <- seq(z[1], z[2], by = 0.5)
      data.frame(trt = trt_name, core = core_id,
                 upper = z[1], lower = z[2], bd_true = mean(bd_fun(d)))
    }))
    res$mult_core_bd <- core_mult
    res
  }

  soc_params <- list(
    high = list(c0 = 3.5, c_end = 0.3, m0 = 1.0),
    mid  = list(c0 = 2.5, c_end = 0.6, m0 = 3.0),
    low  = list(c0 = 2.0, c_end = 0.9, m0 = 12.0)
  )
  if (is.null(c0))    c0    <- soc_params[[soc_prof]]$c0
  if (is.null(c_end)) c_end <- soc_params[[soc_prof]]$c_end
  if (is.null(m0))    m0    <- soc_params[[soc_prof]]$m0

  bd_ref_mult <- replicate(n_cores, r_lognorm_mean1(bd_core_cv))
  bd_t1_mult  <- if (paired) bd_ref_mult else replicate(n_cores, r_lognorm_mean1(bd_core_cv))

  bd_fun_ref_list <- vector("list", n_cores)
  bd_fun_t1_list  <- vector("list", n_cores)
  for (i in seq_len(n_cores)) {
    bd_fun_ref_list[[i]] <- make_bd_fun("ref", bd_ref_mult[i], bd_increase_pct,
                                        z0_cm, z_effect_max, w_transition)
    bd_fun_t1_list[[i]]  <- make_bd_fun("t1",  bd_t1_mult[i],  bd_increase_pct,
                                        z0_cm, z_effect_max, w_transition)
  }

  ref_layers <- do.call(rbind, lapply(seq_len(n_cores), function(i)
    mk_bd_true_layers_from_fun(bd_fun_ref_list[[i]], layers, "ref", i, bd_ref_mult[i])))
  t1_layers  <- do.call(rbind, lapply(seq_len(n_cores), function(i)
    mk_bd_true_layers_from_fun(bd_fun_t1_list[[i]],  layers, "t1",  i, bd_t1_mult[i])))

  out <- rbind(ref_layers, t1_layers)
  out <- out[order(out$trt, out$core, out$upper), ]

  mprof_ref1  <- mk_mprof(bd_fun_ref_list[[1]], z_max = 100, dz = dz_mass)
  m_end       <- max(mprof_ref1$m)
  eff_mult    <- 1 + soc_effect_pct / 100

  c_of_m_ref <- function(m) c_mass_fun(m, c0 = c0, c_end = c_end, m_end = m_end, m0 = m0)

  out$c_true_layer     <- NA_real_
  out$mU_ref           <- NA_real_
  out$mL_ref           <- NA_real_
  out$c_true_esm_stock <- NA_real_

  ref_mass_bounds <- do.call(rbind, lapply(seq_len(n_cores), function(i) {
    mprof_ref <- mk_mprof(bd_fun_ref_list[[i]], z_max = 100, dz = dz_mass)
    data.frame(
      core  = i,
      upper = vapply(layers, function(z) z[1], numeric(1)),
      lower = vapply(layers, function(z) z[2], numeric(1)),
      mU_ref = approx(mprof_ref$z, mprof_ref$m,
                      xout = vapply(layers, function(z) z[1], numeric(1)), rule = 2)$y,
      mL_ref = approx(mprof_ref$z, mprof_ref$m,
                      xout = vapply(layers, function(z) z[2], numeric(1)), rule = 2)$y
    )
  }))

  for (i in seq_len(n_cores)) {
    bnd_i <- ref_mass_bounds[ref_mass_bounds$core == i, ]
    bnd_i <- bnd_i[order(bnd_i$upper), ]

    for (tt in c("ref", "t1")) {
      bd_fun <- if (tt == "ref") bd_fun_ref_list[[i]] else bd_fun_t1_list[[i]]
      mprof  <- mk_mprof(bd_fun, z_max = 100, dz = dz_mass)

      c_layer <- vapply(seq_len(nrow(bnd_i)), function(k) {
        idx   <- mprof$z >= bnd_i$upper[k] & mprof$z < bnd_i$lower[k]
        dm    <- mprof$dm[idx]; m <- mprof$m[idx]; z <- mprof$z[idx]
        cvals <- c_of_m_ref(m)
        if (tt == "t1") {
          if (soc_effect_depth == "full") {
            cvals <- cvals * eff_mult
          } else if (soc_effect_depth == "surface") {
            cvals <- ifelse(z <= 50, cvals * eff_mult, cvals)
          }
        }
        sum(cvals * dm) / sum(dm)
      }, numeric(1))

      c_esm_stock <- vapply(seq_len(nrow(bnd_i)), function(k) {
        m_seq <- seq(bnd_i$mU_ref[k], bnd_i$mL_ref[k], by = dz_mass)
        if (length(m_seq) < 2) m_seq <- c(bnd_i$mU_ref[k], bnd_i$mL_ref[k])
        cvals <- c_of_m_ref(m_seq)
        if (tt == "t1") {
          if (soc_effect_depth == "full") {
            cvals <- cvals * eff_mult
          } else if (soc_effect_depth == "surface") {
            z_seq <- approx(mprof$m, mprof$z, xout = m_seq, rule = 2)$y
            cvals <- ifelse(z_seq <= 50, cvals * eff_mult, cvals)
          }
        }
        if (length(m_seq) < 2) mean(cvals) * (bnd_i$mL_ref[k] - bnd_i$mU_ref[k])
        else sum(cvals * dz_mass)
      }, numeric(1))

      sel <- out$core == i & out$trt == tt
      out$c_true_layer[sel]     <- c_layer
      out$mU_ref[sel]           <- bnd_i$mU_ref
      out$mL_ref[sel]           <- bnd_i$mL_ref
      out$c_true_esm_stock[sel] <- c_esm_stock
    }
  }

  out$bd_obs       <- pmax(out$bd_true + rnorm(nrow(out), 0, bd_lab_sd), 0.6)
  out$c_obs        <- add_lab_noise(out$c_true_layer, c_lab_cv)
  out$tk           <- out$lower - out$upper
  out$c_obs_zstock <- out$c_obs * out$bd_obs * out$tk * 0.1
  out$soc_prof     <- soc_prof
  out
}

# -----------------------------------------------------------------------------
# Long formats (base R — updated to new esm_compare() API)
# -----------------------------------------------------------------------------
get_fd_long <- function(df_sim) {
  data.frame(
    method            = "FD",
    trt               = recode_trt(df_sim$trt),
    core              = df_sim$core,
    upper             = df_sim$upper,
    lower             = df_sim$lower,
    layer_stock_kg_m2 = df_sim$c_obs_zstock,
    stringsAsFactors  = FALSE
  )
}

get_esm_long <- function(df_sim, interp = "spline") {
  df_sim$bd        <- df_sim$bd_obs
  df_sim$soc       <- df_sim$c_obs

  res <- esm_compare(
    data           = df_sim,
    soil_prop      = "soc",
    profile_by     = c("trt", "core"),
    output         = "layers",
    reference_mode = "min",
    interp         = interp,
    soil_prop_unit = "pct",
    som_from_soc   = TRUE,
    soc_unit       = "pct"
  )

  lyr <- res$layer
  data.frame(
    method            = paste0("ESM_", interp),
    trt               = recode_trt(lyr$trt),
    core              = lyr$core,
    upper             = lyr$upper,
    lower             = lyr$lower,
    layer_stock_kg_m2 = lyr$layer_stock_kg_m2,
    stringsAsFactors  = FALSE
  )
}

integrate_to <- function(long_df, max_lower = 100) {
  sub <- long_df[long_df$lower <= max_lower, , drop = FALSE]
  # aggregate sum of layer_stock by method + trt + core
  agg <- aggregate(
    layer_stock_kg_m2 ~ method + trt + core,
    data = sub,
    FUN  = function(x) sum(x, na.rm = TRUE)
  )
  agg
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$script(HTML("
      function initTooltips(){
        var plc = (window.innerWidth < 768) ? 'bottom' : 'left';
        if (window.bootstrap && bootstrap.Tooltip) {
          var els = [].slice.call(document.querySelectorAll('[data-bs-toggle=\"tooltip\"]'));
          els.forEach(function(el){
            var inst = bootstrap.Tooltip.getInstance(el);
            if (inst) inst.dispose();
            new bootstrap.Tooltip(el, { placement: plc, trigger: 'hover focus', container: 'body' });
          });
          return;
        }
        if (window.jQuery && jQuery.fn && jQuery.fn.tooltip) {
          var $els = jQuery('[data-toggle=\"tooltip\"]');
          $els.tooltip('destroy');
          $els.tooltip({ placement: plc, trigger: 'hover focus', container: 'body' });
        }
      }
      $(document).on('shiny:connected', function(){ initTooltips(); });
      $(document).on('shiny:value',     function(){ initTooltips(); });
      $(window).on('resize',            function(){ initTooltips(); });
    "))
  ),

  titlePanel("SOC Sampling Design Simulator"),

  sidebarLayout(
    sidebarPanel(
      h4("Step 0: Choose scenario"),
      with_info(
        radioButtons(
          "mode",
          "Is there a true SOC change between treatments?",
          choices  = c("No (H0) — Type I error / bias" = "null",
                       "Yes (H1) — Power / design"     = "effect"),
          selected = "null"
        ),
        "H0: no true SOC change (Type I error focus). H1: true SOC change (power/design focus)."
      ),

      tags$hr(),
      h4("Monte Carlo"),
      with_info(
        sliderInput("n_rep", "Monte Carlo replicates",
                    min = 10, max = 500, value = 50, step = 10),
        "Number of simulated experiments. Larger values reduce Monte Carlo noise but increase runtime."
      ),
      with_info(
        sliderInput("alpha", "alpha",
                    min = 0.001, max = 0.20, value = 0.05, step = 0.005),
        "Significance threshold used to compute Type I error (H0) or power (H1)."
      ),

      tags$hr(),
      h4("Sampling design"),
      with_info(
        sliderInput("n_cores", "Cores per treatment (single run)",
                    min = 2, max = 10, value = 3, step = 1),
        "Cores per treatment used in the main simulation."
      ),
      with_info(
        checkboxInput("paired", "Paired cores (same BD core multiplier)", value = TRUE),
        "If TRUE, each Treatment core shares the same BD multiplier as its paired Reference core."
      ),
      with_info(
        selectInput("design", "Layering design",
                    choices = c("fine", "coarse"), selected = "fine"),
        "Fine = more layers (0-10, 10-30, 30-50, 50-70, 70-100 cm). Coarse = fewer layers (0-30, 30-50, 50-100 cm)."
      ),
      with_info(
        selectInput("interval", "Integrated interval",
                    choices  = c("0-30" = 30, "0-50" = 50, "0-100" = 100),
                    selected = 30),
        "Stocks are summed up to this depth before testing."
      ),

      tags$hr(),
      h4("Scenario"),
      with_info(
        selectInput("soc_prof", "SOC stratification",
                    choices = c("low", "mid", "high"), selected = "mid"),
        "'high' = strong stratification (SOC decreases sharply with depth). 'low' = weak stratification."
      ),
      with_info(
        sliderInput("bd_pct", "BD increase in treatment (%)",
                    min = 0, max = 30, value = 15, step = 1),
        "Bulk density increase applied to Treatment near the surface (logistic decay with depth)."
      ),

      conditionalPanel(
        condition = "input.mode == 'effect'",
        with_info(
          sliderInput("soc_eff", "True SOC effect in treatment (%)",
                      min = -30, max = 30, value = 5, step = 1),
          "Multiplicative effect on Treatment SOC profile."
        ),
        with_info(
          selectInput(
            "soc_eff_depth", "SOC effect depth",
            choices  = c("Whole profile (0-100 cm)" = "full",
                         "Surface only (0-50 cm)"   = "surface"),
            selected = "full", width = "100%"
          ),
          "Whether the true SOC effect applies to the whole profile or only to the upper 0-50 cm."
        )
      ),
      conditionalPanel(
        condition = "input.mode == 'null'",
        helpText("In H0 mode, the true SOC effect is fixed at 0%.")
      ),

      tags$hr(),
      h4("Noise model"),
      with_info(
        sliderInput("bd_core_cv", "BD core CV",
                    min = 0, max = 0.30, value = 0.07, step = 0.01),
        "Between-core BD variability (lognormal multiplier with mean 1)."
      ),
      with_info(
        sliderInput("bd_lab_sd", "BD lab SD (additive)",
                    min = 0, max = 0.20, value = 0.03, step = 0.005),
        "BD measurement error: BD_obs = BD_true + N(0, SD)."
      ),
      with_info(
        sliderInput("c_lab_cv", "C lab CV",
                    min = 0, max = 0.30, value = 0.05, step = 0.01),
        "SOC analytical error (multiplicative lognormal, mean-preserving)."
      ),

      tags$hr(),
      h4("ESM settings"),
      with_info(
        selectInput("esm_interp", "Interpolation",
                    choices = c("linear", "spline"), selected = "spline"),
        "Interpolation method used for the cumulative SOC vs cumulative mineral mass curve."
      ),

      tags$hr(),
      h4("Power curve (H1 mode only)"),
      conditionalPanel(
        condition = "input.mode == 'effect'",
        with_info(
          sliderInput("cores_min", "Min cores", min = 2, max = 10, value = 2, step = 1),
          "Minimum cores per treatment in the power curve."
        ),
        with_info(
          sliderInput("cores_max", "Max cores", min = 3, max = 10, value = 5, step = 1),
          "Maximum cores per treatment in the power curve."
        ),
        with_info(
          sliderInput("cores_step", "Step", min = 1, max = 5, value = 1, step = 1),
          "Grid step for n_cores in the power curve."
        )
      ),

      tags$hr(),
      uiOutput("run_buttons")
    ),

    mainPanel(
      uiOutput("results_ui"),
      verbatimTextOutput("txt")
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  layers_fine   <- list(c(0,10), c(10,30), c(30,50), c(50,70), c(70,100))
  layers_coarse <- list(c(0,30), c(30,50), c(50,100))

  main_results  <- reactiveVal(NULL)
  power_results <- reactiveVal(NULL)

  # Dynamic buttons
  output$run_buttons <- renderUI({
    if (input$mode == "effect") {
      tagList(
        actionButton("run", "Run simulation"),
        tags$div(style = "margin-top:8px;"),
        actionButton("run_power", "Run power curve"),
        tags$div(style = "margin-top:8px;"),
        actionButton("clear", "Clear results")
      )
    } else {
      tagList(
        actionButton("run", "Run simulation"),
        tags$div(style = "margin-top:8px;"),
        actionButton("clear", "Clear results")
      )
    }
  })

  # Dynamic results UI
  output$results_ui <- renderUI({
    if (input$mode == "null") {
      fluidRow(
        column(6, plotOutput("p_pvals")),
        column(6, plotOutput("p_est"))
      )
    } else {
      tagList(
        fluidRow(
          column(6, plotOutput("p_est_true")),
          column(6, plotOutput("p_power"))
        )
      )
    }
  })

  # Clear results
  observeEvent(input$clear, {
    main_results(NULL); power_results(NULL)
  }, ignoreInit = TRUE)

  # Ground-truth delta for H1 (base R)
  compute_true_delta_single <- function(cfg) {
    req(cfg$mode == "effect")
    layers    <- if (cfg$design == "fine") layers_fine else layers_coarse
    max_lower <- cfg$interval

    df_sim <- simulate_once(
      layers           = layers,
      bd_increase_pct  = cfg$bd_pct,
      soc_effect_pct   = cfg$soc_eff,
      soc_effect_depth = cfg$soc_eff_depth,
      seed       = 123, n_cores = 1, paired = TRUE,
      soc_prof   = cfg$soc_prof,
      bd_core_cv = 0, bd_lab_sd = 0, c_lab_cv = 0
    )

    df_true      <- df_sim
    df_true$bd   <- df_true$bd_true
    df_true$soc  <- df_true$c_true_layer

    res_true <- esm_compare(
      data           = df_true,
      soil_prop      = "soc",
      profile_by     = c("trt", "core"),
      output         = "layers",
      reference_mode = "min",
      interp         = cfg$esm_interp,
      soil_prop_unit = "pct",
      som_from_soc   = TRUE,
      soc_unit       = "pct"
    )

    lyr <- res_true$layer
    lyr$trt <- recode_trt(lyr$trt)
    lyr <- lyr[lyr$lower <= max_lower, , drop = FALSE]

    # stock per trt x core
    stock_by_core <- aggregate(
      layer_stock_kg_m2 ~ trt + core, data = lyr,
      FUN = function(x) sum(x, na.rm = TRUE)
    )

    mean_trt <- mean(stock_by_core$layer_stock_kg_m2[stock_by_core$trt == "Treatment"])
    mean_ref <- mean(stock_by_core$layer_stock_kg_m2[stock_by_core$trt == "Reference"])
    out      <- mean_trt - mean_ref

    if (abs(cfg$soc_eff) < 1e-12) out <- 0
    out
  }

  # Core engine for one Monte Carlo replicate
  run_one_rep <- function(rep_id, n_cores_use, cfg) {
    layers    <- if (cfg$design == "fine") layers_fine else layers_coarse
    max_lower <- cfg$interval

    df_sim <- simulate_once(
      layers           = layers,
      bd_increase_pct  = cfg$bd_pct,
      soc_effect_pct   = cfg$soc_eff,
      soc_effect_depth = cfg$soc_eff_depth,
      seed       = 1000 + rep_id,
      n_cores    = n_cores_use,
      paired     = cfg$paired,
      soc_prof   = cfg$soc_prof,
      bd_core_cv = cfg$bd_core_cv,
      bd_lab_sd  = cfg$bd_lab_sd,
      c_lab_cv   = cfg$c_lab_cv,
      w_transition = 6
    )

    fd_long  <- get_fd_long(df_sim)
    esm_long <- suppressWarnings(get_esm_long(df_sim, interp = cfg$esm_interp))

    int_df <- integrate_to(rbind(fd_long, esm_long), max_lower = max_lower)

    # test by method
    methods <- unique(int_df$method)
    res_list <- lapply(methods, function(m) {
      sub  <- int_df[int_df$method == m, , drop = FALSE]
      test <- safe_lm_trt(sub, y = "layer_stock_kg_m2")
      data.frame(
        rep      = rep_id,
        method   = m,
        estimate = test$estimate,
        p_value  = test$p.value,
        conf.low = test$conf.low,
        conf.high = test$conf.high,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, res_list)
  }

  # Main simulation
  observeEvent(input$run, {
    cfg <- make_cfg(input)
    n   <- cfg$n_rep
    nc  <- cfg$n_cores

    out_list <- future_lapply(
      X          = seq_len(n),
      FUN        = function(i) run_one_rep(i, nc, cfg),
      future.seed = TRUE
    )
    main_results(do.call(rbind, out_list))
  }, ignoreInit = TRUE)

  # Power curve
  observeEvent(input$run_power, {
    cfg <- make_cfg(input)
    req(cfg$mode == "effect")

    cores_seq <- seq(cfg$cores_min, cfg$cores_max, by = cfg$cores_step)
    cores_seq <- cores_seq[cores_seq >= 2]
    if (length(cores_seq) < 1) { power_results(NULL); return() }

    alpha <- cfg$alpha
    n     <- cfg$n_rep

    out_pc <- do.call(rbind, lapply(cores_seq, function(nc2) {
      df_nc <- do.call(rbind, future_lapply(
        X          = seq_len(n),
        FUN        = function(i) run_one_rep(i, nc2, cfg),
        future.seed = TRUE
      ))
      # power by method
      methods <- unique(df_nc$method)
      do.call(rbind, lapply(methods, function(m) {
        sub <- df_nc[df_nc$method == m, , drop = FALSE]
        data.frame(
          method  = m,
          n_cores = nc2,
          power   = mean(sub$p_value < alpha, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }))
    }))

    power_results(out_pc)
  }, ignoreInit = TRUE)

  # H0 plots
  output$p_pvals <- renderPlot({
    req(main_results())
    df <- main_results()
    ggplot(df, aes(x = p_value)) +
      geom_histogram(bins = 25) +
      geom_vline(xintercept = input$alpha, linetype = "dashed") +
      coord_cartesian(xlim = c(0, 1)) +
      facet_wrap(~ method, ncol = 1) +
      labs(x = "p-value", y = "Count", title = "p-value distribution") +
      theme_bw()
  })

  output$p_est <- renderPlot({
    req(main_results())
    df <- main_results()
    ggplot(df, aes(x = estimate)) +
      geom_histogram(bins = 25) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      facet_wrap(~ method, ncol = 1) +
      labs(x = expression(hat(Delta) ~ Stock ~ (kg ~ m^{-2}) ~ "(Treatment - Reference)"),
           y = "Count", title = "Estimated treatment effect distribution (H0)") +
      theme_bw()
  })

  # H1 plots
  output$p_power <- renderPlot({
    req(power_results())
    df <- power_results()
    ggplot(df, aes(x = n_cores, y = power, color = method, linetype = method)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      ylim(0, 1) +
      labs(x = "Cores per treatment (n)", y = "Empirical power  P(p < alpha)",
           title = "Power curve (by method)", color = "Method") +
      theme_bw()
  })

  output$p_est_true <- renderPlot({
    req(main_results())
    df         <- main_results()
    true_delta <- compute_true_delta_single(make_cfg(input))
    ggplot(df, aes(x = estimate)) +
      geom_histogram(bins = 25) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = true_delta, color = "red", linewidth = 1) +
      facet_wrap(~ method, ncol = 1, scales = "fixed") +
      labs(x = expression(hat(Delta) ~ Stock ~ (kg ~ m^{-2}) ~ "(Treatment - Reference)"),
           y = "Count", title = "Estimated effect distribution (H1)",
           subtitle = "Red line = ground-truth \u0394 stock (ESM-truth)") +
      theme_bw()
  })

  # Text summary
  output$txt <- renderText({
    req(main_results())
    cfg   <- make_cfg(input)
    df    <- main_results()
    alpha <- cfg$alpha

    scenario <- paste0(
      "Scenario\n",
      "- Mode: ", ifelse(cfg$mode == "null", "H0 (no true SOC change)", "H1 (true SOC change)"), "\n",
      "- Design: ", cfg$design, " | interval: 0-", cfg$interval, " cm",
      " | n_cores: ", cfg$n_cores, " | paired: ", cfg$paired, "\n",
      "- SOC profile: ", cfg$soc_prof, " | BD increase (%): ", cfg$bd_pct, "\n",
      "- True SOC effect (%): ", cfg$soc_eff,
      " | depth: ", cfg$soc_eff_depth, "\n",
      "- Noise: bd_core_cv = ", cfg$bd_core_cv,
      ", bd_lab_sd = ", cfg$bd_lab_sd,
      ", c_lab_cv = ", cfg$c_lab_cv, "\n",
      "- Alpha: ", alpha, " | Monte Carlo reps: ", cfg$n_rep, "\n"
    )

    methods <- unique(df$method)

    if (cfg$mode == "null") {
      lines <- paste(sapply(methods, function(m) {
        sub <- df[df$method == m, , drop = FALSE]
        paste0("* ", m,
               ": Type I error = ", round(mean(sub$p_value < alpha, na.rm = TRUE), 2),
               " | Bias = ",        round(mean(sub$estimate, na.rm = TRUE), 2),
               " | SD(estimate) = ", round(sd(sub$estimate, na.rm = TRUE), 2))
      }), collapse = "\n")

      paste0(scenario, "\n", lines, "\n\nNotes\n",
             "- In H0 mode, Type I error should be close to alpha if well-calibrated.\n")

    } else {
      true_delta <- compute_true_delta_single(cfg)

      lines <- paste(sapply(methods, function(m) {
        sub <- df[df$method == m, , drop = FALSE]
        paste0("* ", m,
               ": Power = ",      round(mean(sub$p_value < alpha, na.rm = TRUE), 2),
               " | Bias = ",      round(mean(sub$estimate - true_delta, na.rm = TRUE), 2),
               " | RMSE = ",      round(sqrt(mean((sub$estimate - true_delta)^2, na.rm = TRUE)), 2),
               " | Coverage95 = ",round(mean(sub$conf.low <= true_delta &
                                               sub$conf.high >= true_delta, na.rm = TRUE), 2))
      }), collapse = "\n")

      pc_note <- "\n\nPower curve has not been run yet."
      if (!is.null(power_results())) {
        pc_note <- paste0("\n\nPower curve computed for n_cores = ",
                          cfg$cores_min, " to ", cfg$cores_max,
                          " (step ", cfg$cores_step, ").")
      }

      paste0(scenario,
             "\nGround-truth \u0394 stock (ESM-truth) = ", round(true_delta, 3), " kg m-2\n\n",
             lines, pc_note, "\n\nNotes\n",
             "- Power is P(p < alpha).\n",
             "- Bias/RMSE/Coverage computed against ESM ground truth.\n",
             "- Coverage close to 0.95 indicates well-calibrated confidence intervals.\n")
    }
  })
}

shinyApp(ui, server)
