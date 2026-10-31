# ============================================================
# MAYV_threshold_sensitivity.R
# Run after: MAYV_lhs.R (writes MAYV_lhs_ensemble_high.rds)
# ============================================================
DEFS_ONLY <- TRUE
suppressMessages({library(writexl); source("common/common.R"); source("MAYV/MAYV_engine.R")})

E <- readRDS("outputs/MAYV_lhs_ensemble_high.rds")
N <- E$N; A <- E$A; T_weeks <- E$T_weeks; ND <- length(E$gamma)
season <- { v <- as.numeric(readRDS("outputs/caldas_hybrid_season.rds")); v / max(v) } 
target_age <- rep(0, A); target_age[4:8] <- 1
start_pre <- 17L; immun_delay <- 2                        # 2025-W40

# Vaccine parameters are read from the engine, not redrawn
Mres <- readRDS("outputs/MAYV_engine_results_high.rds")
stopifnot(!is.null(Mres$del_d))          # re-run MAYV_engine.R if this fails
cov_d <- Mres$cov_d; veb_d <- Mres$veb_d; del_d <- Mres$del_d; dly_d <- Mres$delay_d
stopifnot(length(cov_d) == ND, length(veb_d) == ND, length(del_d) == ND, length(dly_d) == ND)

wk_lab <- function(i) { o <- 23 + i
  if (o <= 53) sprintf("2025-W%02d", o) else sprintf("2026-W%02d", o - 53) }
q3  <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
fmt <- function(x, dp = 0) { q <- q3(x); sprintf("%s (%s-%s)",
        formatC(round(q[1], dp), big.mark=",", format="f", digits=dp),
        formatC(round(q[2], dp), big.mark=",", format="f", digits=dp),
        formatC(round(q[3], dp), big.mark=",", format="f", digits=dp)) }

# baseline and vaccinated symptomatic totals
one <- function(i, R0, seed_week) {
  imm <- E$immune_frac[i]; Rimm <- rep(imm, A); sus <- N * (1 - imm)
  bt  <- R0 * E$gamma[i] * season
  sim <- function(cov, vb) seirv_vaccinated_MAYV(
    T_weeks, A, N, Rimm, E$I0_total * sus/sum(sus), bt, E$sigma[i], E$gamma[i], E$rho[i],
    target_age, cov, del_d[i], start_pre + dly_d[i], VE_inf = 0, VE_block = vb,
    immun_delay = immun_delay, prop_symp = E$prop_symp[i], E0 = E$E0, seed_week = seed_week)
  b <- sim(0, 0); v <- sim(cov_d[i], veb_d[i]); psi <- E$prop_symp[i]
  c(base  = psi * sum(b$new_infections),
    vacc  = sum(colSums(psi * v$new_infections * (1 - veb_d[i] * v$coverage_frac))),
    doses = sum(v$total_used_age))
}
sweep_cfg <- function(R0, seed_week) {
  m <- vapply(seq_len(ND), function(i) one(i, R0, seed_week), numeric(3))
  av <- m["base", ] - m["vacc", ]
  list(base = m["base", ], averted = av, doses = m["doses", ],
       pct = 100 * av / m["base", ],
       # per 100,000 doses
       per100k = 1e5 * av / m["doses", ])
}

# 1. transmission potential
R0_MAIN <- 2.47                      # median of the sampled high range
SCEN <- list(low  = "1.1-1.3 (incidental transmission, outside Amazon basin)",
             high = "2.1-2.9 (urban transmission, Amazon-basin intensity)")

scen_sampled <- function(sc) {
  Es <- readRDS(sprintf("outputs/MAYV_lhs_ensemble_%s.rds", sc))
  Ms <- readRDS(sprintf("outputs/MAYV_engine_results_%s.rds", sc))
  stopifnot(identical(Es$R0_scenario, Ms$R0_scenario), isTRUE(Es$R0_sampled))
  nd <- length(Es$gamma)
  cv <- as.numeric(Ms$cov_d); ve <- as.numeric(Ms$veb_d)
  dl <- as.numeric(Ms$del_d); dy <- as.numeric(Ms$delay_d)
  base <- av <- dose <- numeric(nd)
  mR <- pR <- wk <- numeric(nd)
  for (i in seq_len(nd)) {
    imm <- Es$immune_frac[i]; SNi <- 1 - imm
    R0i <- Es$R0[i]
    mR[i] <- R0i * mean(season) * SNi           # mean R_eff over the 52 weeks
    pR[i] <- R0i * max(season)  * SNi           # peak R_eff (season is peak-normalised)
    wk[i] <- sum(R0i * season * SNi > 1)        # weeks with R_eff > 1
    Rimm <- rep(imm, A); sus <- N * (1 - imm)
    sim <- function(cov, vb) seirv_vaccinated_MAYV(
      T_weeks, A, N, Rimm, Es$I0_total * sus/sum(sus), R0i * Es$gamma[i] * season,
      Es$sigma[i], Es$gamma[i], Es$rho[i], target_age, cov, dl[i],
      start_pre + dy[i], VE_inf = 0, VE_block = vb, immun_delay = immun_delay,
      prop_symp = Es$prop_symp[i], E0 = Es$E0, seed_week = Es$seed_week)
    r <- sim(cv[i], ve[i]); ni <- r$new_infections; psi <- Es$prop_symp[i]
    base[i] <- sum(colSums(psi * ni))
    av[i]   <- base[i] - sum(colSums(psi * ni * (1 - ve[i] * r$coverage_frac)))
    dose[i] <- Ms$target_pop_elig * cv[i]
  }
  cat(sprintf("  %s scenario done\n", sc))
  data.frame(scenario = SCEN[[sc]],
             sampled_R0 = reff(Es$R0),
             mean_R_eff = reff(mR),
             peak_R_eff = reff(pR),
             weeks_R_eff_above_1 = sprintf("%d (%d-%d)", median(wk),
                                           quantile(wk, .025), quantile(wk, .975)),
             baseline_symptomatic = fmt(base),
             averted_symptomatic  = fmt(av),
             averted_per_100k_doses = fmt(1e5 * av / dose),
             pct_averted = fmt(100 * av / base, 1),
             stringsAsFactors = FALSE)
}
SN_d <- 1 - E$immune_frac; SN <- median(SN_d)
reff <- function(x, d = 3) sprintf("%.*f (%.*f-%.*f)", d, median(x), d,
                                   quantile(x, .025), d, quantile(x, .975))
cat(sprintf("Propagating %d draws over the two sampled scenarios...\n", ND))
tp <- do.call(rbind, lapply(names(SCEN), scen_sampled))

# 2. seeding week 
SEED_GRID <- c(1, 9, 17, 25)
cat("Propagating seeding weeks...\n")
sd <- do.call(rbind, lapply(SEED_GRID, function(sw) {
  r <- sweep_cfg(R0_MAIN, sw); cat(sprintf("  seed week %d done\n", sw))
  data.frame(seed_week_index = sw, seed_week = wk_lab(sw),
             season_at_seeding = round(season[sw], 3),
             R_eff_at_seeding  = round(R0_MAIN * season[sw] * SN, 3),
             baseline_symptomatic = fmt(r$base),
             averted_symptomatic  = fmt(r$averted),
             averted_per_100k_doses = fmt(r$per100k),
             pct_averted          = fmt(r$pct, 2),
             stringsAsFactors = FALSE) }))

# 3. seasonal envelope 
v_mean1 <- as.numeric(readRDS("outputs/caldas_hybrid_season.rds"))
SPLICE  <- 40L                                   # 2026-W10; weeks 1-40 CHIKV beta, 41-52 rainfall
env <- data.frame(
  week_index      = seq_along(v_mean1),
  epi_week        = vapply(seq_along(v_mean1), wk_lab, character(1)),
  source          = ifelse(seq_along(v_mean1) <= SPLICE,
                           "fitted CHIKV beta_t", "CHIRPS rainfall climatology"),
  envelope_mean1  = round(v_mean1, 4),
  season_peaknorm = round(v_mean1 / max(v_mean1), 4),
  R0_t_at_scenario    = round(R0_MAIN * v_mean1 / max(v_mean1), 3),
  R_eff_t_at_scenario = round(R0_MAIN * v_mean1 / max(v_mean1) * SN, 3),
  stringsAsFactors = FALSE)
env_summary <- data.frame(quantity = c(
  "mean of the stored mean-1 envelope", "max of the stored mean-1 envelope",
  "mean(season) after peak-normalisation", "derivation",
  "min(season)", "weeks with season > 0.5", "peak week",
  "source split", "why the tail is rainfall"),
  value = c(
  sprintf("%.6f (1 by construction)", mean(v_mean1)),
  sprintf("%.6f", max(v_mean1)),
  sprintf("%.4f", mean(v_mean1)/max(v_mean1)),
  sprintf("mean(season) = mean(v)/max(v) = 1/%.4f = %.4f -- a consequence of peak-normalising, not a separate estimate", max(v_mean1), 1/max(v_mean1)),
  sprintf("%.4f", min(v_mean1)/max(v_mean1)),
  sprintf("%d of %d", sum(v_mean1/max(v_mean1) > 0.5), length(v_mean1)),
  sprintf("week %d (%s)", which.max(v_mean1), wk_lab(which.max(v_mean1))),
  sprintf("weeks 1-%d (%.0f%%) fitted CHIKV beta_t; weeks %d-%d (%.0f%%) rainfall climatology, rescaled to join continuously at 2026-W10",
          SPLICE, 100*SPLICE/52, SPLICE+1, length(v_mean1), 100*(52-SPLICE)/52),
  "The fitted CHIKV beta_t tail is confounded with susceptible depletion, not short of data: 39.5% of observed CHIKV cases fall at or after the join. Depletion alone explains the case decline, so beta_t stays flat (0.638 -> 0.630) and carries no seasonal signal. Rainfall is exogenous to the epidemic."),
  stringsAsFactors = FALSE)

write_xlsx(list(transmission_potential = tp, seeding_week = sd,
                envelope = env, envelope_summary = env_summary),
           "outputs/MAYV_threshold_sensitivity.xlsx")
cat("\nWrote MAYV_threshold_sensitivity.xlsx\n")
print(tp[, c("scenario","sampled_R0","mean_R_eff","peak_R_eff","weeks_R_eff_above_1","baseline_symptomatic","averted_symptomatic")], row.names = FALSE)
cat("\n"); print(sd[, c("seed_week","R_eff_at_seeding","baseline_symptomatic","averted_symptomatic","pct_averted")], row.names = FALSE)
