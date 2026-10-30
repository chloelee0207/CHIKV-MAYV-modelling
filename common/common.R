# ============================================================
# common.R
# ------------------------------------------------------------
# Single source of truth for code shared across the CHIKV engine (CHIKV_engine.R) and MAYV_engine.R. 
# Source it once near the top of each.
# ============================================================
library(readxl)
suppressMessages({library(dplyr); library(tidyr)})
dir.create("outputs", showWarnings = FALSE)

# The 2025 epi calendar has a Semana 53, so 2025 contributes 30 weeks (W24..W53):
#   2025-W24 = index 1  ; 2025-W53 = index 30   (2025 weeks: index = week - 23)
#   2026-W01 = index 31                          (2026 weeks: index = 30 + week)
week_to_index <- function(year, week) ifelse(year == 2025, week - 23L, 30L + week)

# Summarise a Monte Carlo draw vector
fmtq <- function(v, d = 0) {
  q <- quantile(v, c(.5, .025, .975), na.rm = TRUE)
  f <- function(x) formatC(round(x, d), big.mark = ",", format = "f", digits = d)
  sprintf("%s (%s - %s)", f(q[1]), f(q[2]), f(q[3]))
}

qs <- function(v) as.numeric(quantile(v, c(.5, .025, .975), na.rm = TRUE))

if (!exists("age_weight")) age_weight <- 1

# Folder structure
DRAWS_DIR <- "outputs/lhs_draws"

# Estimate burden
# hr: hosp rate
# cv: CFR-by-age
# w re-weights the age distribution for deaths only, redistributing symptomatic cases across ages 
burden <- function(out, hr = hosp_rate, cv = cfr_vec, w = age_weight) {
  symp_age <- rowSums(out$new_symptomatic)
  symp     <- sum(symp_age)
  symp_w   <- symp_age * w
  sw       <- sum(symp_w)
  symp_dw  <- if (sw > 0) symp_w * (symp / sw) else symp_age   # renormalise: keep total
  c(infections       = sum(out$new_infections),
    symptomatic      = symp,
    hospitalisations = symp * hr,
    deaths           = sum(symp_dw * cv))
}

# Read the disease-progression Beta(alpha, beta) (Kang et al. 2025 Table S4, disease_progression.xlsx) 
load_burden_params <- function(A,
                               dp_path  = "data/disease_progression.xlsx",
                               dp_sheet = "disease_progression") {
  dp <- read_excel(dp_path, sheet = dp_sheet)
  names(dp)[1:10] <- c("parameter", "group", "median", "ui_lo", "ui_hi",
                       "dist", "p1", "alpha", "p2", "beta")
  dp <- dp[dp$dist == "Beta", ]

  get_beta_ab <- function(param, group_regex = NULL) {
    d <- dp[dp$parameter == param, ]
    if (!is.null(group_regex)) d <- d[grepl(group_regex, d$group), ]
    lo <- suppressWarnings(as.numeric(sub(".*\\[\\s*([0-9]+).*", "\\1", d$group)))
    if (!all(is.na(lo))) d <- d[order(lo), ]
    list(a = d$alpha, b = d$beta)
  }

  ps <- get_beta_ab("Probability of symptomatic cases among infections", "^Overall$")
  ps_a <- ps$a; ps_b <- ps$b

  hp <- get_beta_ab("Probability of hospitalisation among symptomatic cases")
  hosp_a <- hp$a; hosp_b <- hp$b
  hosp_rate <- hosp_a / (hosp_a + hosp_b)

  ch <- get_beta_ab("Probability of death among hospitalised cases")
  cn <- get_beta_ab("Probability of death among non-hospitalised cases")
  cfr_hosp_a <- ch$a; cfr_hosp_b <- ch$b
  cfr_nonh_a <- cn$a; cfr_nonh_b <- cn$b
  stopifnot(length(ps_a) == 1, length(hosp_a) == 1,
            length(cfr_hosp_a) == 9, length(cfr_nonh_a) == 9)

  cfr_h_mean <- cfr_hosp_a / (cfr_hosp_a + cfr_hosp_b)
  cfr_n_mean <- cfr_nonh_a / (cfr_nonh_a + cfr_nonh_b)
  cfr_band   <- hosp_rate * cfr_h_mean + (1 - hosp_rate) * cfr_n_mean  # death per symptomatic, by band

  age_to_band <- c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9)   # 12 model groups -> 9 decadal bands
  stopifnot(length(age_to_band) == A)
  cfr_vec <- cfr_band[age_to_band]

  list(ps_a = ps_a, ps_b = ps_b,
       hosp_a = hosp_a, hosp_b = hosp_b, hosp_rate = hosp_rate,
       cfr_hosp_a = cfr_hosp_a, cfr_hosp_b = cfr_hosp_b,
       cfr_nonh_a = cfr_nonh_a, cfr_nonh_b = cfr_nonh_b,
       cfr_h_mean = cfr_h_mean, cfr_n_mean = cfr_n_mean, cfr_band = cfr_band,
       age_to_band = age_to_band, cfr_vec = cfr_vec)
}

# ------------------------------------------------------------
# Caldas Novas age-stratified cases (the "ca_combined" SINAN sheet).
# Single source of truth for cases within the study period
# ------------------------------------------------------------
load_caldas_age_cases <- function(path = "data/weekly_case.xlsx", sheet = "ca_combined") {
  age_levels <- c("<1 Ano", "1-4", "5-9", "10-14", "15-19",
                  "20-39", "40-59", "60-64", "65-69", "70-79", "80 e +")
  ca_long <- read_excel(path, sheet = sheet) |>
    rename(Year = Ano, semana = Semana) |>
    mutate(Year = as.integer(Year), week = as.integer(sub("Semana ", "", semana))) |>
    filter((Year == 2025 & week >= 24) | (Year == 2026 & week <= 22)) |>
    pivot_longer(all_of(age_levels), names_to = "age_group", values_to = "cases") |>
    mutate(cases = ifelse(is.na(cases), 0, cases))

  week_grid <- bind_rows(tibble(Year = 2025L, week = 24:53),
                         tibble(Year = 2026L, week = 1:22)) |>
    arrange(Year, week) |>
    mutate(week_index = row_number(),
           week_label = paste0(Year, "-W", sprintf("%02d", week)))

  caldas_obs <- week_grid |>
    left_join(ca_long |> group_by(Year, week) |>
                summarise(cases = sum(cases), .groups = "drop"),
              by = c("Year", "week")) |>
    mutate(cases = ifelse(is.na(cases), 0, cases))

  ca_age <- tidyr::expand_grid(week_grid, age_group = factor(age_levels, age_levels)) |>
    left_join(ca_long |> mutate(age_group = factor(age_group, age_levels)),
              by = c("Year", "week", "age_group")) |>
    mutate(cases = ifelse(is.na(cases), 0, cases)) |>
    dplyr::select(week_index, week_label, Year, week, age_group, cases)

  age_totals <- ca_age |> group_by(age_group) |>
    summarise(cases = sum(cases), .groups = "drop")

  group_to_bands <- list("<1 Ano" = 1, "1-4" = 1, "5-9" = 1, "10-14" = 2, "15-19" = 2,
                         "20-39" = c(3, 4), "40-59" = c(5, 6),
                         "60-64" = 7, "65-69" = 7, "70-79" = 8, "80 e +" = 9)
  obs_band <- numeric(9)
  for (g in age_levels) {
    b <- group_to_bands[[g]]
    n <- age_totals$cases[age_totals$age_group == g]
    obs_band[b] <- obs_band[b] + n / length(b)
  }

  list(observed_cases = caldas_obs$cases, caldas_obs = caldas_obs,
       ca_age = ca_age, age_totals = age_totals,
       obs_band_prop = obs_band / sum(obs_band), age_levels = age_levels)
}

# Age re-weighting vector w_a (length A) that maps the model's baseline infection age split to the observed case age split, per decadal band.
compute_age_weight <- function(inf_age_model, obs_band_prop, age_to_band) {
  mod_band <- as.numeric(tapply(inf_age_model, age_to_band, sum)[as.character(seq_along(obs_band_prop))])
  mod_band[is.na(mod_band)] <- 0
  mod_prop <- mod_band / sum(mod_band)
  w_band   <- ifelse(mod_prop > 0, obs_band_prop / mod_prop, 1)
  w_band[age_to_band]
}

# ============================================================
# Shared SEIRV engine + DALY parameter
# ============================================================

# VE_inf moves S->immune (infection blocking)
# VE_block scales symptomatic among the covered (disease blocking)
# 1 channel per arm (VE_block = 0 when VE_inf > 0
# symptomatic formula applies VE_block to all infections
seirv_vaccinated <- function(
    T_weeks, 
    A, 
    N, 
    R_init_prop, 
    I0, 
    E0, 
    base_beta, 
    sigma, 
    gamma, 
    rho,
    target_age, 
    total_coverage, 
    weekly_delivery_speed, 
    delay,
    VE_inf = 0,
    VE_block = 0,
    immun_delay = 2, 
    prop_symp = 0.5242478, 
    sub_steps = 7) {
  pmax0 <- function(x) pmax(0, x); N_total <- sum(N); dt <- 1/sub_steps
  S <- E <- I <- R <- V <- matrix(0, A, T_weeks)
  V_covered <- vacc_delayed <- coverage_frac <- matrix(0, A, T_weeks)
  new_infections <- new_symptomatic <- matrix(0, A, T_weeks)
  target_idx <- which(target_age == 1); target_pop <- sum(N[target_idx])
  total_supply <- target_pop * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  total_avail_age <- rep(0, A); total_avail_age[target_idx] <- total_supply * (N[target_idx]/target_pop)
  total_used_age <- rep(0, A); unvaccinated <- N
  S_now <- pmax0(N - I0 - E0 - R_init_prop*N); E_now <- E0; I_now <- I0
  R_now <- R_init_prop*N; V_now <- rep(0, A)
  for (t in 1:T_weeks) {
    prev_V_covered <- if (t > 1) V_covered[, t-1] else rep(0, A)
    if (t - immun_delay >= 1) {
      effective_dose <- vacc_delayed[, t-immun_delay]
      immunized <- round(VE_inf * effective_dose)
      V_covered[, t] <- prev_V_covered + effective_dose
    } else { immunized <- rep(0, A); V_covered[, t] <- prev_V_covered }
    S_now <- pmax0(S_now - immunized); V_now <- V_now + immunized
    coverage_frac[, t] <- V_covered[, t] / N
    if (t >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(ceiling(weekly_dose_total*(N[a]/target_pop)), rem,
                     unvaccinated[a], total_avail_age[a]-total_used_age[a])
        if (alloc > 0) {
          prop_S <- if (N[a] > 0) S_now[a]/N[a] else 0
          vacc_to_S <- round(alloc*prop_S)
          vacc_delayed[a, t] <- vacc_to_S
          total_used_age[a] <- total_used_age[a] + alloc
          unvaccinated[a] <- unvaccinated[a] - alloc; rem <- rem - alloc
        }
      }
    }
    new_I_week <- rep(0, A); beta_t <- base_beta[t]
    for (k in 1:sub_steps) {
      foi <- beta_t * sum(I_now)/N_total
      new_E <- foi*S_now*dt; 
      new_I <- sigma*E_now*dt; 
      new_R <- gamma*I_now*dt
      S_now <- pmax0(S_now-new_E); 
      E_now <- pmax0(E_now+new_E-new_I)
      I_now <- pmax0(I_now+new_I-new_R); 
      R_now <- pmax0(R_now+new_R)
      new_I_week <- new_I_week + new_I
    }
    S[,t]<-S_now; E[,t]<-E_now; I[,t]<-I_now; R[,t]<-R_now; V[,t]<-V_now
    new_infections[,t] <- new_I_week
    new_symptomatic[,t] <- prop_symp*new_I_week*(1 - VE_block*coverage_frac[,t])
  }
  list(new_infections=new_infections, 
       new_symptomatic=new_symptomatic,
       new_reported=rho*new_symptomatic, 
       V_covered=V_covered, 
       total_used_age=total_used_age)
}

# Read disease-progression distribution hyperparameters from disease_progression.xlsx
load_daly_params <- function(dp_path = "data/disease_progression.xlsx",
                             dp_sheet = "disease_progression") {
  dp <- read_excel(dp_path, sheet = dp_sheet)
  names(dp)[1:10] <- c("parameter","group","median","ui_lo","ui_hi","dist","p1","v1","p2","v2")
  row <- function(param, grp = NULL) {
    d <- dp[dp$parameter == param, ]
    if (!is.null(grp)) d <- d[grepl(grp, d$group), ]
    d
  }
  beta_ab  <- function(param, grp = NULL) { d <- row(param, grp); list(a = d$v1, b = d$v2) }
  lnorm_ms <- function(param, grp = NULL) { d <- row(param, grp); list(m = d$v1, s = d$v2, med = d$median) }
  dw_mm  <- beta_ab("Disability weight for mild and moderate chikungunya")
  dw_sev <- beta_ab("Disability weight for severe chikungunya")
  dw_chr <- beta_ab("Disability weight for chronic chikungunya")
  du_mm  <- lnorm_ms("Duration of illness for mild and moderate chikungunya (years)")
  du_sev <- lnorm_ms("Duration of illness for severe chikungunya (years)")
  du_sub <- lnorm_ms("Duration of illness for sub-acute chikungunya (years)")
  du_chr <- lnorm_ms("Duration of illness for chronic chikungunya (years)")
  # No separate disability weight is published for the sub-acute phase; assume chronic dw, sampled from the same draw
  le <- row("Remaining life-years")
  lo <- as.numeric(sub(".*\\[\\s*([0-9]+).*", "\\1", le$group)); le <- le[order(lo), ]
  le_ms <- list(m = le$v1, s = le$v2, med = le$median)
  stopifnot(length(le_ms$m) == 9)
  p14_y <- beta_ab("Probability of recovery within 14 days after onset of symptoms", "< 40")
  p14_o <- beta_ab("Probability of recovery within 14 days after onset of symptoms", "> 40")
  p90_y <- beta_ab("Probability of recovery within 90 days after acute period", "< 40")
  p90_o <- beta_ab("Probability of recovery within 90 days after acute period", "> 40")
  # The five recovery rows are marginal proportions of one cohort, not conditional
  # probabilities: they sum to 1 in each age group. So the chronic proportion is the sum of the 6m/12m/30m rows
  p6_y  <- beta_ab("Probability of recovery within 6 months after sub-acute period", "< 40")
  p6_o  <- beta_ab("Probability of recovery within 6 months after sub-acute period", "> 40")
  p12_y <- beta_ab("Probability of recovery within 12 months after 6 months of chronicity", "< 40")
  p12_o <- beta_ab("Probability of recovery within 12 months after 6 months of chronicity", "> 40")
  p30_y <- beta_ab("Probability of recovery within 30 months after 12 months of chronicity", "< 40")
  p30_o <- beta_ab("Probability of recovery within 30 months after 12 months of chronicity", "> 40")
  list(dw_mm=dw_mm, dw_sev=dw_sev, dw_chr=dw_chr,
       du_mm=du_mm, du_sev=du_sev, du_sub=du_sub, du_chr=du_chr, le=le_ms,
       p14_y=p14_y, p14_o=p14_o, p90_y=p90_y, p90_o=p90_o,
       p6_y=p6_y, p6_o=p6_o, p12_y=p12_y, p12_o=p12_o, p30_y=p30_y, p30_o=p30_o)
}
