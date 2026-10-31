# ============================================================
# MAYV_owsa.R
# Run order:  MAYV_lhs.R -> MAYV_engine.R (once) -> this file
# ============================================================
suppressMessages({library(dplyr); library(ggplot2); library(writexl)})

DEFS_ONLY <- TRUE
source("MAYV/MAYV_engine.R")        # seirv_vaccinated_MAYV, outcome_one, severity/DALY parameters
E <- readRDS("outputs/MAYV_lhs_ensemble.rds")
M <- readRDS("outputs/MAYV_engine_results.rds")

if (!identical(E$R0_scenario, M$R0_scenario))
  stop("MAYV_lhs_ensemble.rds is scenario '", E$R0_scenario,
       "' but MAYV_engine_results.rds is '", M$R0_scenario,
       "'. Re-run MAYV_lhs.R then MAYV_engine.R for one scenario before this script.")

GAMMA <- E$base_gamma; SIGMA <- E$base_sigma
PSYMP <- E$base_prop_symp        # MAYV-specific median (~0.787)
EVAL_WIN <- 1:T_weeks

# ------------------------------------------------------------
# 1. Central values and one-way bounds
# ------------------------------------------------------------
ob <- M$outbreak           
R0_BASE <- unname(M$R0_median)    
stopifnot(is.finite(R0_BASE))

R0_MARKS <- c(low = 1.20, high = 2.47)             # scenario medians
R0_SWEEP <- sort(unique(round(c(seq(1.1, 3.6, by = 0.1), unname(R0_MARKS)), 4)))

# Prior immunity is 0
IMM_BASE <- 0
BASE <- list(R0 = R0_BASE, imm = IMM_BASE, rho = 0.25, ve = 4/12,
             cov = 0.30, deliv = 0.10, delay = 2, immun = 2, env = "hybrid",
             psymp = PSYMP,  
             latent = 1/SIGMA)  

BOUNDS <- list(
  R0     = unname(c(E$R0_lo, E$R0_hi)),                # the scenario's sampled range
  ve     = c(0.109, 0.610),                            # Kostecki 2026, Beta(4,8) 95% UI
  cov   = c(0.20, 0.40),
  deliv = c(0.09, 0.11),
  delay = c(1, 3),
  immun = c(1, 3))

PAR_LAB <- c(R0 = "R0",
             rho = "Reporting rate", ve = "Vaccine efficacy",
             cov = "Vaccine coverage", deliv = "Weekly delivery speed",
             delay = "Delay in deployment", immun = "Time to immunity",
             env = "Seasonal envelope (inherited)")

# Seasonal envelopes
env_of <- function(tag) {
  v <- switch(tag,
    hybrid = readRDS("outputs/caldas_hybrid_season.rds"),
    beta   = readRDS("outputs/caldas_beta_season.rds"),
    rain   = { r <- readRDS("outputs/caldas_rain_season.rds"); r / mean(r) })
  v <- as.numeric(v); stopifnot(length(v) == T_weeks)
  v / max(v) 
}

# ------------------------------------------------------------
# 2. baseline + disease-blocking arm, deterministic
# ------------------------------------------------------------
run_scenario <- function(p) {
  season <- env_of(p$env)
  Rimm   <- rep(p$imm, A)  
  sus    <- N * (1 - p$imm)
  I0i    <- E$I0_total * sus / sum(sus)
  beta_t <- p$R0 * GAMMA * season
  psi    <- if (!is.null(p$psymp)) p$psymp else PSYMP   # symptomatic fraction 
  sg     <- if (!is.null(p$latent)) 1/p$latent else SIGMA   # latent period (week) -> rate
  sim <- function(cov, vb)
    seirv_vaccinated_MAYV(T_weeks, A, N, Rimm, I0i, beta_t, sg, GAMMA, p$rho,
      target_age, cov, p$deliv, start_pre + p$delay, VE_inf = 0, VE_block = vb,
      immun_delay = p$immun, prop_symp = psi, E0 = E$E0, seed_week = E$seed_week)
  b <- sim(0, 0); v <- sim(p$cov, p$ve)
  symp_b <- sum(colSums(psi * b$new_infections)[EVAL_WIN])
  symp_v <- sum(colSums(psi * v$new_infections * (1 - p$ve * v$coverage_frac))[EVAL_WIN])
  inf    <- sum(b$new_infections[, EVAL_WIN])
  pool   <- sum(N * (1 - p$imm))
  # Doses = the whole eligible 18-59 cohort times coverage
  doses  <- target_pop_elig * p$cov
  list(symptomatic_base = symp_b, infections_base = inf, averted = symp_b - symp_v,
       # reported = rho x true
       reported_symptomatic = p$rho * symp_b,
       reported_averted     = p$rho * (symp_b - symp_v),
       attack = 100 * inf / pool,
       hosp_averted = (symp_b - symp_v) * hosp_rate,
       doses = doses, averted_per_100k = 1e5 * (symp_b - symp_v) / doses)
}

# ------------------------------------------------------------
# 3. Base case + each parameter at both bounds
# ------------------------------------------------------------
cat(sprintf("MAYV OWSA: base peak R0 = %.2f (FIXED, '%s' scenario; %d draws, no conditioning)\n",
            R0_BASE, M$R0_scenario, length(ob)))
base_run <- run_scenario(BASE)
cat(sprintf("  base: attack %.1f%% | symptomatic %.0f | averted %.0f\n",
            base_run$attack, base_run$symptomatic_base, base_run$averted))

rows <- list()
for (nm in names(BOUNDS)) {
  for (side in c("lower", "upper")) {
    p <- BASE; p[[nm]] <- BOUNDS[[nm]][if (side == "lower") 1 else 2]
    r <- run_scenario(p)
    rows[[length(rows)+1]] <- data.frame(
      parameter = unname(PAR_LAB[nm]), par_key = nm, bound = side,
      value = as.character(p[[nm]]),
      symptomatic_base = r$symptomatic_base, averted = r$averted,
      hosp_averted = r$hosp_averted, attack_pct = r$attack,
      pct_reduction = 100 * r$averted / r$symptomatic_base, row.names = NULL)
  }
  cat("  ", nm, " done\n", sep = "")
}
owsa <- do.call(rbind, rows)

# ---- structural sensitivities ----
srow <- function(input, setting, p) {
  r <- run_scenario(p)
  data.frame(input = input, setting = setting,
             symptomatic_base = r$symptomatic_base, averted = r$averted,
             averted_per_100k_doses = r$averted_per_100k,
             reported_symptomatic = r$reported_symptomatic,
             reported_averted = r$reported_averted,
             attack_pct = r$attack, pct_reduction = 100*r$averted/r$symptomatic_base,
             row.names = NULL)
}
struct <- rbind(
  do.call(rbind, lapply(unname(R0_MARKS), function(x)
    srow("Peak R0 (scenario median)", sprintf("%.2f", x), modifyList(BASE, list(R0 = x))))),
  do.call(rbind, lapply(c("beta", "hybrid", "rain"), function(x)
    srow("Seasonal envelope (structure)",
         c(beta = "pure CHIKV beta", hybrid = "hybrid (base case)",
           rain = "pure rainfall")[[x]], modifyList(BASE, list(env = x))))),
  do.call(rbind, lapply(list(c(PSYMP,     sprintf("%.1f%% (MAYV-specific, base case)", 100*PSYMP)),
                             c(0.5242478, "52.4% (CHIKV value)")), function(x)
    srow("Symptomatic fraction (structure)", x[2],
         modifyList(BASE, list(psymp = as.numeric(x[1])))))),
  do.call(rbind, lapply(list(
      c(unname(BASE$latent), sprintf("%.2f d (Martins et al. 2020, 7-12 d; base case)", 7*BASE$latent)),
      c(sqrt(35)/7, "5.92 d (viraemia/symptom band 5-7 d)"),
      c(sqrt(6)/7,  "2.45 d (Diagne et al. 2020, ~0-6 d)"),
      c(3.0/7,      "3.0 d (Caicedo et al. 2021 estimated mean)"),
      c(7.5/7,      "7.5 d (CDC 1-14 d uniform median, earlier base case)"),
      c(12.0/7,     "12 d (earlier assumption)")), function(x)
    srow("Latent period (source)", x[2],
         modifyList(BASE, list(latent = as.numeric(x[1])))))),
  do.call(rbind, lapply(list(c(0.25,  "25% (base case)"),
                             c(0.389, "38.9% (alternative point estimate)")), function(x)
    srow("Reporting rate (structure)", x[2],
         modifyList(BASE, list(rho = as.numeric(x[1])))))))
rownames(struct) <- NULL

# ------------------------------------------------------------
# 4. Tornado
# ------------------------------------------------------------
b0 <- base_run$averted
sw <- owsa |> group_by(parameter) |> summarise(swing = diff(range(averted)), .groups = "drop")
owsa$parameter <- factor(owsa$parameter, levels = sw$parameter[order(sw$swing)])
p_t <- ggplot(owsa, aes(y = parameter)) +
  geom_vline(xintercept = b0, linetype = "dashed", colour = "grey45") +
  geom_segment(aes(x = b0, xend = averted, yend = parameter, colour = bound),
               linewidth = 6, alpha = .85) +
  scale_colour_manual(values = c(lower = "#d6604d", upper = "#4393c3"),
                      labels = c(lower = "Lower", upper = "Upper"), name = "Bound") +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Symptomatic cases averted", y = NULL) +
  theme_bw(11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("outputs/MAYV_owsa_symptomatic.png", p_t, width = 8.5, height = 4.6, dpi = 130)
print(p_t)

# ------------------------------------------------------------
# 4b. R0 response curve
# ------------------------------------------------------------
r0_response <- do.call(rbind, lapply(R0_SWEEP, function(x) {
  r <- run_scenario(modifyList(BASE, list(R0 = x)))
  data.frame(R0 = x, infections_base = r$infections_base,
             symptomatic_base = r$symptomatic_base, averted = r$averted,
             averted_per_100k_doses = r$averted_per_100k,
             attack_pct = r$attack,
             pct_reduction = 100 * r$averted / max(r$symptomatic_base, 1e-9),
             row.names = NULL) }))
cat("\n=== R0 response curve (deterministic, median inputs) ===\n")
print(transform(r0_response, infections_base = round(infections_base),
                symptomatic_base = round(symptomatic_base),
                averted = round(averted),
                averted_per_100k_doses = round(averted_per_100k_doses),
                attack_pct = round(attack_pct, 3),
                pct_reduction = round(pct_reduction, 1)), row.names = FALSE)

scen_pts <- subset(r0_response, R0 %in% round(R0_MARKS, 10))
p_r0 <- ggplot(r0_response, aes(R0)) +
  geom_line(aes(y = pmax(symptomatic_base, 0.5), colour = "Baseline symptomatic"), linewidth = .7) +
  geom_line(aes(y = pmax(averted, 0.5),          colour = "Averted by vaccine"),   linewidth = .7) +
  geom_point(data = scen_pts, aes(y = pmax(symptomatic_base, 0.5)), size = 1.8, colour = "#c0392b") +
  geom_vline(xintercept = unname(R0_MARKS), linetype = "dotted", colour = "grey40") +
  annotate("text", x = R0_MARKS["low"],  y = Inf, label = sprintf("low\n%.2f",  R0_MARKS["low"]),  vjust = 1.3, size = 2.5, colour = "grey30") +
  annotate("text", x = R0_MARKS["high"], y = Inf, label = sprintf("high\n%.2f", R0_MARKS["high"]), vjust = 1.3, size = 2.5, colour = "grey30") +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Baseline symptomatic" = "#c0392b",
                                 "Averted by vaccine"   = "#2c7fb8"), name = NULL) +
  labs(x = "Wet-season peak R0", y = "Symptomatic cases (log scale)") +
  theme_bw(9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
# Supplementary figure
ggsave("outputs/MAYV_r0_response.png", p_r0, width = 190 / 25.4, height = 3.8, dpi = 1000)
cat("Saved MAYV_r0_response.png\n")

# ------------------------------------------------------------
# 4c. probabilistic R0 increment
# ------------------------------------------------------------
if (!exists("R0_PSA_RUN")) R0_PSA_RUN <- TRUE
r0_response_psa <- NULL
if (R0_PSA_RUN) {
  nd      <- length(E$gamma)
  cov_v   <- as.numeric(M$cov_d)      # engine's sampled coverage of eligible 18-59
  ve_v    <- as.numeric(M$veb_d)      # engine's sampled disease-blocking VE
  delay_v <- as.numeric(M$delay_d)    # engine's sampled deployment delay (week)
  deliv_c <- M$deliv_median           # weekly delivery speed, held at median
  stopifnot(length(cov_v) == nd, length(ve_v) == nd, length(delay_v) == nd)
  season_psa <- env_of(BASE$env)
  q3 <- function(x) unname(quantile(x, c(.5, .025, .975), na.rm = TRUE))

  cat(sprintf("\n4c. Probabilistic R0 sweep: %d R0 values x %d draws = %d runs...\n",
              length(R0_SWEEP), nd, length(R0_SWEEP)*nd))
  t_start <- Sys.time()
  r0_response_psa <- do.call(rbind, lapply(seq_along(R0_SWEEP), function(k) {
    x <- R0_SWEEP[k]
    sb <- av <- at <- nf <- numeric(nd)
    dose_v <- target_pop_elig * cov_v  
    for (i in seq_len(nd)) {
      Rimm <- rep(E$immune_frac[i], A)
      sus  <- N * (1 - E$immune_frac[i])
      I0i  <- E$I0_total * sus / sum(sus)
      # VE_inf = 0 -> infections are vaccine-invariant
      r <- seirv_vaccinated_MAYV(T_weeks, A, N, Rimm, I0i, x * E$gamma[i] * season_psa,
             E$sigma[i], E$gamma[i], E$rho[i], target_age, cov_v[i], deliv_c,
             start_pre + delay_v[i], VE_inf = 0, VE_block = ve_v[i],
             immun_delay = BASE$immun, prop_symp = E$prop_symp[i],
             E0 = E$E0, seed_week = E$seed_week)
      ninf <- r$new_infections; psi <- E$prop_symp[i]
      sb[i] <- sum(colSums(psi * ninf)[EVAL_WIN])
      av[i] <- sb[i] - sum(colSums(psi * ninf * (1 - ve_v[i] * r$coverage_frac))[EVAL_WIN])
      nf[i] <- sum(ninf[, EVAL_WIN])
      at[i] <- 100 * nf[i] / sum(sus)
    }
    pr <- 100 * av / pmax(sb, 1e-12)
    p1e5 <- 1e5 * av / pmax(dose_v, 1e-12)   # averted per 100,000 doses
    qs <- q3(sb); qa <- q3(av); qp <- q3(pr); qt <- q3(at); qn <- q3(nf); q5 <- q3(p1e5)
    cat(sprintf("   R0 %.1f done (%2d/%2d)  symptomatic %.0f [%.0f, %.0f]\n",
                x, k, length(R0_SWEEP), qs[1], qs[2], qs[3]))
    data.frame(R0 = x,
               infections_med = qn[1], infections_lo = qn[2], infections_hi = qn[3],
               symptomatic_med = qs[1], symptomatic_lo = qs[2], symptomatic_hi = qs[3],
               averted_med = qa[1], averted_lo = qa[2], averted_hi = qa[3],
               averted_per_100k_med = q5[1], averted_per_100k_lo = q5[2],
               averted_per_100k_hi = q5[3],
               pct_reduction_med = qp[1], pct_reduction_lo = qp[2], pct_reduction_hi = qp[3],
               attack_med = qt[1], attack_lo = qt[2], attack_hi = qt[3],
               row.names = NULL) }))
  cat(sprintf("   sweep finished in %.1f min\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

  fmt3 <- function(m, lo, hi, d = 0)
    sprintf(paste0("%.", d, "f (%.", d, "f - %.", d, "f)"), m, lo, hi)
  r0_table <- data.frame(
    `Peak R0`                       = sprintf("%.2f", r0_response_psa$R0),
    `Attack rate (%)`               = fmt3(r0_response_psa$attack_med, r0_response_psa$attack_lo,
                                           r0_response_psa$attack_hi, 2),
    `Baseline infections`           = fmt3(r0_response_psa$infections_med, r0_response_psa$infections_lo,
                                           r0_response_psa$infections_hi, 0),
    `Baseline symptomatic cases`    = fmt3(r0_response_psa$symptomatic_med, r0_response_psa$symptomatic_lo,
                                           r0_response_psa$symptomatic_hi, 0),
    `Symptomatic cases averted`     = fmt3(r0_response_psa$averted_med, r0_response_psa$averted_lo,
                                           r0_response_psa$averted_hi, 0),
    `Averted per 100,000 doses`     = fmt3(r0_response_psa$averted_per_100k_med,
                                           r0_response_psa$averted_per_100k_lo,
                                           r0_response_psa$averted_per_100k_hi, 0),
    `Symptomatic cases reduced (%)` = fmt3(r0_response_psa$pct_reduction_med, r0_response_psa$pct_reduction_lo,
                                           r0_response_psa$pct_reduction_hi, 1),
    check.names = FALSE, stringsAsFactors = FALSE)
  write.csv(r0_table, "outputs/MAYV_r0_response_table.csv", row.names = FALSE)
  cat("\n=== R0 response TABLE: median (95% UI) over ", nd, " draws ===\n", sep = "")
  print(r0_table, row.names = FALSE)

  # figure 
  p_r0p <- ggplot(r0_response_psa, aes(R0)) +
    geom_ribbon(aes(ymin = pmax(symptomatic_lo, .5), ymax = pmax(symptomatic_hi, .5),
                    fill = "Baseline symptomatic"), alpha = .22) +
    geom_ribbon(aes(ymin = pmax(averted_lo, .5), ymax = pmax(averted_hi, .5),
                    fill = "Averted by vaccine"), alpha = .22) +
    geom_line(aes(y = pmax(symptomatic_med, .5), colour = "Baseline symptomatic"), linewidth = .7) +
    geom_line(aes(y = pmax(averted_med, .5),     colour = "Averted by vaccine"),   linewidth = .7) +
    geom_vline(xintercept = unname(R0_MARKS), linetype = "dotted", colour = "grey40") +
    annotate("text", x = R0_MARKS["low"],  y = Inf, label = sprintf("low\n%.2f",  R0_MARKS["low"]),  vjust = 1.3, size = 2.5, colour = "grey30") +
    annotate("text", x = R0_MARKS["high"], y = Inf, label = sprintf("high\n%.2f", R0_MARKS["high"]), vjust = 1.3, size = 2.5, colour = "grey30") +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(values = c("Baseline symptomatic" = "#c0392b",
                                   "Averted by vaccine"   = "#2c7fb8"), name = NULL,
                        aesthetics = c("colour", "fill")) +
    labs(x = "Wet-season peak R0", y = "Symptomatic cases (log scale)") +
    theme_bw(9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  # Supplementary figure
  ggsave("outputs/MAYV_r0_response_psa.png", p_r0p, width = 190 / 25.4, height = 3.8, dpi = 1000)
  cat("Saved MAYV_r0_response_psa.png and MAYV_r0_response_table.csv\n")
}

# ------------------------------------------------------------
# 5. Export
# ------------------------------------------------------------
write_xlsx(c(list(base_case = data.frame(symptomatic_base = base_run$symptomatic_base,
                                         averted = base_run$averted, attack_pct = base_run$attack,
                                         pct_reduction = 100*base_run$averted/base_run$symptomatic_base),
                owsa = owsa, structural_sensitivity = struct,
                r0_response = r0_response),
             if (!is.null(r0_response_psa))
               list(r0_response_psa = r0_response_psa, r0_response_table = r0_table)),
           "outputs/MAYV_owsa.xlsx")
saveRDS(list(owsa = owsa, struct = struct, base = base_run, BASE = BASE, BOUNDS = BOUNDS,
             R0_BASE = R0_BASE, R0_SWEEP = R0_SWEEP, r0_response = r0_response,
             r0_response_psa = r0_response_psa,
             R0_scenario = M$R0_scenario), "outputs/MAYV_owsa.rds")

cat("\n=== Swing in symptomatic cases averted ===\n")
print(owsa |> group_by(parameter) |>
        summarise(low = min(averted), high = max(averted), swing = high - low, .groups = "drop") |>
        arrange(desc(swing)) |> as.data.frame(), row.names = FALSE, digits = 5)
cat("\n=== Structural sensitivities (reported separately, not tornado rows) ===\n")
print(struct, row.names = FALSE, digits = 5)
cat("\nWrote MAYV_owsa.xlsx, MAYV_owsa.rds, MAYV_owsa_symptomatic.png\n")

# ------------------------------------------------------------
# 7. vaccine efficacy x coverage
# ------------------------------------------------------------
VE_GRID  <- seq(0.1, 1.0, by = 0.1)
COV_GRID <- seq(0.1, 1.0, by = 0.1)

sweep_season <- env_of(BASE$env)
sweep_Rimm   <- rep(BASE$imm, A)
sweep_sus    <- N * (1 - BASE$imm)
sweep_I0     <- E$I0_total * sweep_sus / sum(sweep_sus)
sweep_beta   <- BASE$R0 * GAMMA * sweep_season
sweep_sim <- function(cov) seirv_vaccinated_MAYV(T_weeks, A, N, sweep_Rimm, sweep_I0,
  sweep_beta, SIGMA, GAMMA, BASE$rho, target_age, cov, BASE$deliv,
  start_pre + BASE$delay, VE_inf = 0, VE_block = 0, immun_delay = BASE$immun,
  prop_symp = PSYMP, E0 = E$E0, seed_week = E$seed_week)

sw_base  <- sweep_sim(0)
symp_tot <- sum(colSums(PSYMP * sw_base$new_infections))
cf_list  <- lapply(COV_GRID, function(cv) sweep_sim(cv)$coverage_frac)

grid_df <- do.call(rbind, lapply(seq_along(COV_GRID), function(j) {
  cf <- cf_list[[j]]
  do.call(rbind, lapply(VE_GRID, function(ve) {
    symp_v <- sum(colSums(PSYMP * sw_base$new_infections * (1 - ve * cf)))
    data.frame(coverage = COV_GRID[j], ve = ve,
               averted = symp_tot - symp_v,
               pct = 100 * (symp_tot - symp_v) / symp_tot) })) }))

p_hm <- ggplot(grid_df, aes(factor(100*coverage), factor(100*ve), fill = pct)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.1f", pct)), size = 2.5, colour = "grey15") +
  scale_fill_distiller(palette = "Spectral", name = "% of total\nsymptomatic\ncases averted") +
  labs(x = "Vaccine coverage (% of eligible population aged 18-59)",
       y = "Vaccine efficacy (%)") +
  theme_bw(9) + theme(panel.grid = element_blank())
# Fig. S8
ggsave("outputs/MAYV_ve_coverage_heatmap.png", p_hm, width = 190 / 25.4, height = 5.0, dpi = 1000)
print(p_hm)

write_xlsx(list(ve_coverage_grid = grid_df), "outputs/MAYV_ve_coverage_sweep.xlsx")
near <- which.min(abs(grid_df$ve - BASE$ve) + abs(grid_df$coverage - BASE$cov))
cat(sprintf("\nVE x coverage sweep: %d cells | cell nearest the base case (VE %.0f%%, cov %.0f%%) = %.2f%% averted | max (VE 100%%, cov 100%%) = %.2f%%\n",
            nrow(grid_df), 100*grid_df$ve[near], 100*grid_df$coverage[near],
            grid_df$pct[near], max(grid_df$pct)))
cat("Saved MAYV_ve_coverage_heatmap.png and MAYV_ve_coverage_sweep.xlsx\n")
