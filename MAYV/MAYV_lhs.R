# ============================================================
# Latin hypercube sampling (LHS)
# ============================================================
suppressMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2)})
dir.create("outputs", showWarnings = FALSE) 

# ------------------------------------------------------------
# 1. SEIR
# ------------------------------------------------------------
seir_baseline_MAYV <- function(
    T_weeks, A, N, R_init_prop, I0, base_beta, sigma, gamma, rho, prop_symp,
    sub_steps = 7, E0 = rep(0, A), seed_week = 1) {
  pmax0 <- function(x) pmax(0, x); N_total <- sum(N); dt <- 1/sub_steps
  S <- E <- I <- R <- matrix(0, A, T_weeks)
  new_infections <- new_symptomatic <- matrix(0, A, T_weeks)
  S_now <- pmax0(N - E0 - R_init_prop*N); E_now <- E0; I_now <- rep(0, A); R_now <- R_init_prop*N
  for (t in 1:T_weeks) {
    new_I_week <- rep(0, A); beta_t <- base_beta[t]
    if (t == seed_week) { I_now <- I_now + I0; S_now <- pmax0(S_now - I0) }
    for (k in 1:sub_steps) {
      foi <- beta_t * sum(I_now)/N_total
      new_E <- foi*S_now*dt; new_I <- sigma*E_now*dt; new_R <- gamma*I_now*dt
      S_now <- pmax0(S_now - new_E); E_now <- pmax0(E_now + new_E - new_I)
      I_now <- pmax0(I_now + new_I - new_R); R_now <- pmax0(R_now + new_R)
      new_I_week <- new_I_week + new_I
    }
    S[,t]<-S_now; E[,t]<-E_now; I[,t]<-I_now; R[,t]<-R_now
    new_infections[,t] <- new_I_week; new_symptomatic[,t] <- prop_symp*new_I_week
  }
  list(new_infections=new_infections, new_symptomatic=new_symptomatic,
       new_reported=rho*new_symptomatic)
}

# ------------------------------------------------------------
# 2. Population in Caldas Novas)
# ------------------------------------------------------------
age_df <- read_excel("data/population.xlsx", sheet = "prop_immune")
age_df <- as.data.frame(age_df[tolower(age_df$municipality) == "caldas novas", ])
stopifnot(nrow(age_df) > 0)
A <- nrow(age_df)
pop_2022_total <- 98622; pop_2025_total <- 106820
growth_r <- log(pop_2025_total/pop_2022_total)/3
N        <- age_df$pop_num * exp(growth_r*3)

T_weeks <- 52 

# ------------------------------------------------------------
# 3. hybrid seasonal envelope
# ------------------------------------------------------------
season_mean1 <- readRDS("outputs/caldas_hybrid_season.rds")
stopifnot(length(season_mean1) == T_weeks, abs(mean(season_mean1) - 1) < 1e-6)
r0_is_peak <- TRUE
season <- if (r0_is_peak) season_mean1 / max(season_mean1) else season_mean1

# ------------------------------------------------------------
# 4. introduce infectious person
# ------------------------------------------------------------
I0_total  <- if (exists("MAYV_I0"))   MAYV_I0   else 1     # seed size
E0        <- rep(0, A)
seed_week <- if (exists("MAYV_SEED")) MAYV_SEED else 1     # seed week index; 1 = window open (2025-W24)

# ------------------------------------------------------------
# 5. Priors
# ------------------------------------------------------------
cal <- as.data.frame(read_excel("data/model_calibration.xlsx", sheet = 1))
cal <- cal[cal$Group == "MAYV" & !is.na(cal$Median), ]
row_for <- function(k) cal[grepl(k, cal$Parameter, ignore.case = TRUE), ][1, ]
sd_of   <- function(r) (r[["95% UI upper"]] - r[["95% UI lower"]]) / (2*1.96)

gr <- row_for("gamma"); sr <- row_for("sigma") # gamma sampled lognormal
g_ml <- as.numeric(gr[["Value 1"]]); g_sl <- as.numeric(gr[["Value 2"]])
stopifnot(tolower(gr[["Uncertainty distribution"]]) == "lognormal", is.finite(g_ml), is.finite(g_sl))
g_m <- exp(g_ml)
p_m <- sr$Median; p_sd <- sd_of(sr)
if (!exists("MAYV_LATENT")) MAYV_LATENT <- "mayv7_12"
stopifnot(MAYV_LATENT %in% c("mayv7_12", "mayv5_7", "mayv3_11", "diagne1_6", "cdc1_14", "caicedo3d", "long12d"))
# Latent-period priors (weeks); MAYV_LATENT picks one. Only the base case feeds the published results
LAT_SPEC <- list(
  # Base case: Martins et al. 2020, 7-12 d, as the 95% interval on the population mean
  mayv7_12  = list(dist = "lognormal", lo = 7/7,      hi = 12/7),
  # Diagne et al. 2020, ~0-6 d (read as 1-6 d)
  diagne1_6 = list(dist = "lognormal", lo = 1/7,      hi = 6/7),
  # "Usually 5-7 days" band (the 3-11 d full range is between-person variation, see mayv3_11)
  mayv5_7   = list(dist = "lognormal", lo = 5/7,      hi = 7/7),
  # Full 3-11 d range, lognormal fitted to both the 3-11 d range and the 5-7 d band
  mayv3_11  = list(dist = "lognormal", lo = 3/7,      hi = 11/7,
                   ml = mean(log(c(3,5,7,11)/7)),
                   sl = sum(qnorm(c(.025,.25,.75,.975)) * log(c(3,5,7,11))) /
                        sum(qnorm(c(.025,.25,.75,.975))^2)),
  # CDC 1-14 d, a hard range, so sampled uniformly
  cdc1_14   = list(dist = "uniform",   lo = 1/7,     hi = 14/7),
  # Caicedo et al. 2021 fitted mean, 3.0 d (95% CrI 2.4-4.1)
  caicedo3d = list(dist = "lognormal", lo = 2.4/7,   hi = 4.1/7),
  # Earlier 12-day assumption
  long12d   = list(dist = "lognormal", lo = 11.31/7, hi = 12.69/7))
lat <- LAT_SPEC[[MAYV_LATENT]]
lat_lo <- lat$lo; lat_hi <- lat$hi
lat_ml <- lat_sl <- NA_real_
if (lat$dist == "lognormal") {
  lat_ml <- if (!is.null(lat$ml)) lat$ml else (log(lat_lo) + log(lat_hi)) / 2
  lat_sl <- if (!is.null(lat$sl)) lat$sl else (log(lat_hi) - log(lat_lo)) / (2 * 1.96)
}
lat_q <- function(u) {
  if (lat$dist == "uniform") qunif(u, lat_lo, lat_hi) else qlnorm(u, lat_ml, lat_sl)
}
p_m <- lat_q(0.5) 
cat(sprintf("Latent period '%s': %s on [%.2f, %.2f] d -> median %.4f wk (%.2f d)\n",
            MAYV_LATENT, lat$dist, 7*lat_lo, 7*lat_hi, p_m, 7*p_m))
rab <- c(a = 20, b = 60)                        # rho ~ Beta(20,60)
# proportion symptomatic
ps_ml  <- (log(48/71) + log(33/36)) / 2
ps_sl  <- (log(33/36) - log(48/71)) / (2 * 1.96)
ps_cap <- plnorm(1, ps_ml, ps_sl)


if (!exists("R0_SCENARIO")) R0_SCENARIO <- "high"   # "high" = 2.1-2.9 | "low" = 1.1-1.3

R0_RANGE <- list(low  = c(1.1, 1.3),
                 high = c(2.1, 2.9)) 
stopifnot(R0_SCENARIO %in% names(R0_RANGE))
R0_LO <- R0_RANGE[[R0_SCENARIO]][1]; R0_HI <- R0_RANGE[[R0_SCENARIO]][2]
R0_meanlog <- (log(R0_LO) + log(R0_HI)) / 2
R0_sdlog   <- (log(R0_HI) - log(R0_LO)) / (2 * 1.96)
R0_median  <- exp(R0_meanlog) 

# assume population is naive, so 0% immunity
imm_base <- 0
cat("Prior immunity: 0% -- population fully susceptible\n")

cat(sprintf("SAMPLED priors:\n  gamma     ~ logN(%.4f, %.4f) /wk -> period %.2f d (%.2f-%.2f)\n  latent    ~ %s on [%.2f, %.2f] d  -> median %.2f d\n  rho       ~ Beta(%d, %d)\n  prop_symp ~ logN(%.4f, %.4f) trunc at 1 -> median %.3f, 95%% %.3f-%.3f\n",
            g_ml, g_sl, 7/exp(g_ml), 7/qlnorm(.975, g_ml, g_sl), 7/qlnorm(.025, g_ml, g_sl),
            lat$dist, 7*lat_lo, 7*lat_hi, 7*p_m,
            rab["a"], rab["b"],
            ps_ml, ps_sl, exp(ps_ml), qlnorm(.025, ps_ml, ps_sl), qlnorm(.975, ps_ml, ps_sl)))
cat(sprintf("R0 ~ lognormal on [%.1f, %.1f] (%s scenario, %s): median %.3f\n",
            R0_LO, R0_HI, R0_SCENARIO, if (r0_is_peak) "seasonal-peak R_eff" else "annual-mean", R0_median))

# ------------------------------------------------------------
# Seasonality summary
# ------------------------------------------------------------
season_mean <- mean(season)
cat(sprintf("Seasonality: season(t) peak %.3f, MEAN %.4f, min %.3f; %d of %d weeks above 0.5\n",
            max(season), season_mean, min(season), sum(season > 0.5), length(season)))
# S/N = 1 at the start (no prior immunity), so R_eff(t) = R0(t) until depletion sets in.
cat(sprintf("Implied MEAN R0(t) = %.3f (peak %.3f); %d of %d weeks with R0(t) > 1\n",
            R0_median*season_mean, R0_median*max(season),
            sum(R0_median*season > 1), length(season)))

# ------------------------------------------------------------
# 6. One forward run -> weekly reported / infections + summary scalars
#    imm = flat prior-immune fraction; sets S(0), susceptible_pop and the seed split
# ------------------------------------------------------------
run_draw <- function(R0, g, s, r, ps, imm) {
  Rimm <- rep(imm, A)   
  sus  <- N * (1 - Rimm)                         # susceptibles by age
  I0i  <- I0_total * sus / sum(sus)              # single seed split by susceptibility
  out  <- seir_baseline_MAYV(T_weeks, A, N, Rimm, I0i, R0*g*season,
                             sigma = s, gamma = g, rho = r, prop_symp = ps,
                             E0 = E0, seed_week = seed_week)
  rep_wk <- colSums(out$new_reported); inf_wk <- colSums(out$new_infections)
  list(rep = rep_wk, inf = inf_wk,
       tot_rep = sum(rep_wk), tot_inf = sum(inf_wk),
       peak_rep_wk = which.max(rep_wk), peak_rep = max(rep_wk),
       attack = 100 * sum(inf_wk) / sum(sus), immune = 100 * imm)
}

# Baseline at median inputs
PSYMP_MED <- exp(ps_ml)                          # median symptomatic fraction
base <- run_draw(R0_median, g_m, 1/p_m, rab["a"]/sum(rab), PSYMP_MED, imm_base)
cat(sprintf("Baseline R0=%.2f (%s scenario, median inputs, immune %.1f%%): total infections %.0f | total reported %.1f | peak reported wk %d (%.1f) | attack %.2f%%\n",
            R0_median, R0_SCENARIO, base$immune, base$tot_inf, base$tot_rep, base$peak_rep_wk, base$peak_rep, base$attack))
R0_other <- exp(mean(log(R0_RANGE[[setdiff(names(R0_RANGE), R0_SCENARIO)]])))
base_other <- run_draw(R0_other, g_m, 1/p_m, rab["a"]/sum(rab), PSYMP_MED, imm_base)
cat(sprintf("Reference R0=%.2f (other scenario, same scaling): total infections %.0f | total reported %.1f | attack %.3f%%\n",
            R0_other, base_other$tot_inf, base_other$tot_rep, base_other$attack))

# ------------------------------------------------------------
# 7. Latin Hypercube, sampled: gamma, sigma, rho, prop_symp and R0
# ------------------------------------------------------------
set.seed(2024); n <- 1000
lhs_col <- function(n) (sample.int(n) - runif(n)) / n
U   <- sapply(1:6, function(j) lhs_col(n))
gam <- qlnorm(U[,1], g_ml, g_sl) 
sig <- 1/lat_q(U[,2]) 
rho <- qbeta(U[,3], rab["a"], rab["b"])
psy <- qlnorm(U[,4]*ps_cap, ps_ml, ps_sl)
imm <- rep(0, n)
R0v <- qlnorm(U[,6], R0_meanlog, R0_sdlog) 

rep_mat <- inf_mat <- matrix(NA_real_, n, T_weeks)
tot_rep <- tot_inf <- peak_rep_wk <- peak_rep <- attack <- R0peak <- immune <- rep(NA_real_, n)
cat("Forward-simulating", n, "LHS draws...\n")
for (i in 1:n) {
  d <- tryCatch(run_draw(R0v[i], gam[i], sig[i], rho[i], psy[i], imm[i]), error = function(e) NULL)
  if (is.null(d)) next
  rep_mat[i,] <- d$rep; inf_mat[i,] <- d$inf
  tot_rep[i] <- d$tot_rep; tot_inf[i] <- d$tot_inf
  peak_rep_wk[i] <- d$peak_rep_wk; peak_rep[i] <- d$peak_rep; attack[i] <- d$attack
  immune[i] <- d$immune; R0peak[i] <- R0v[i] * max(season)
  if (i %% 100 == 0) cat("  ", i, "/", n, "\n")
}
ok <- which(is.finite(tot_inf) & is.finite(tot_rep))
cat(sprintf("Kept %d / %d finite draws.\n", length(ok), n))

# ------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------
q3   <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
band <- function(M) apply(M[ok, , drop = FALSE], 2, quantile, c(.025, .5, .975), na.rm = TRUE)
cat("\n=========== PROPAGATED MAYV RESULTS (", length(ok), " draws, R0 SAMPLED ",
    sprintf("%.1f-%.1f", R0_LO, R0_HI), ") ===========\n", sep = "")
cat("Bands below carry gamma / sigma / rho / prop_symp AND R0 uncertainty.\n")
cat(sprintf("Prior immunity:      propagated %.1f%% [%.1f%%, %.1f%%]\n",
            q3(immune[ok])[1], q3(immune[ok])[2], q3(immune[ok])[3]))
cat(sprintf("R0 at seasonal peak: SAMPLED %.2f [%.2f, %.2f] (lognormal on the %s range)\n",
            median(R0v[ok]), quantile(R0v[ok],.025), quantile(R0v[ok],.975), R0_SCENARIO))
cat(sprintf("Attack rate (of susceptibles): baseline %.2f%% -> propagated %.2f%% [%.2f%%, %.2f%%]\n",
            base$attack, q3(attack[ok])[1], q3(attack[ok])[2], q3(attack[ok])[3]))
cat(sprintf("Total infections:  baseline %.0f -> propagated %.0f [%.0f, %.0f]\n",
            base$tot_inf, q3(tot_inf[ok])[1], q3(tot_inf[ok])[2], q3(tot_inf[ok])[3]))
cat(sprintf("Total reported:    baseline %.1f -> propagated %.1f [%.1f, %.1f]\n",
            base$tot_rep, q3(tot_rep[ok])[1], q3(tot_rep[ok])[2], q3(tot_rep[ok])[3]))
cat(sprintf("Busiest reported week (count): baseline %.1f -> propagated %.1f [%.1f, %.1f]\n",
            base$peak_rep, q3(peak_rep[ok])[1], q3(peak_rep[ok])[2], q3(peak_rep[ok])[3]))

# ------------------------------------------------------------
# 9. Plots
# ------------------------------------------------------------
wk_num  <- function(idx) ifelse(idx <= 30, idx + 23, idx - 30) 
tick_idx <- c(7, 17, 27, 40, 50) 
x_ticks  <- data.frame(week_index = tick_idx, week = wk_num(tick_idx))
year_break <- 30.5    
wet_start <- which(season_mean1 >= 1)[1]   
wet_end   <- tail(which(season_mean1 >= 1), 1)
season_layers <- list(
  annotate("rect", xmin = 0.5,       xmax = wet_start,   ymin = -Inf, ymax = Inf, fill = "#f4ead7", alpha = 0.55),
  annotate("rect", xmin = wet_start, xmax = wet_end,     ymin = -Inf, ymax = Inf, fill = "#cfe6f2", alpha = 0.55),
  annotate("rect", xmin = wet_end,   xmax = T_weeks+0.5, ymin = -Inf, ymax = Inf, fill = "#f4ead7", alpha = 0.55))

weeks <- 1:T_weeks
rb <- band(rep_mat); ib <- band(inf_mat)

p_rep <- ggplot(data.frame(week = weeks, lo = rb[1,], med = rb[2,], hi = rb[3,], base = base$rep)) +
  season_layers +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(week, ymin = lo, ymax = hi), fill = "grey40", alpha = 0.30) +
  geom_line(aes(week, med), colour = "grey25", linewidth = 1) +
  geom_line(aes(week, base), colour = "#d6604d", linewidth = 0.8, linetype = "dashed") +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  labs(x = "Week", y = "Predicted reported MAYV cases",
       title = "Hypothetical MAYV outbreak: propagated 95% band",
       subtitle = sprintf("Median (solid) + 95%% band from gamma/sigma/rho/prop_symp/R0 priors (R0 ~ %.1f-%.1f); baseline at median inputs (dashed)", R0_LO, R0_HI)) +
  theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5),
                       plot.subtitle = element_text(hjust = 0.5, size = 9),
                       panel.grid.minor = element_blank())
ggsave("outputs/MAYV_lhs_reported.png", p_rep, width = 8, height = 4.5, dpi = 120)

p_inf <- ggplot(data.frame(week = weeks, lo = ib[1,], med = ib[2,], hi = ib[3,])) +
  season_layers +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(week, ymin = lo, ymax = hi), fill = "#a8d1e7", alpha = 0.5) +
  geom_line(aes(week, med), colour = "#3182bd", linewidth = 1) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  labs(x = "Week", y = "True MAYV infections (all)",
       title = "Hypothetical MAYV outbreak: true infections, propagated 95% band") +
  theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5),
                       panel.grid.minor = element_blank())
ggsave("outputs/MAYV_lhs_infections.png", p_inf, width = 8, height = 4.5, dpi = 120)

# ------------------------------------------------------------
# 10. Save
# ------------------------------------------------------------
if (!exists("DRAWS_DIR")) DRAWS_DIR <- "outputs/lhs_draws"
dir.create(DRAWS_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(draw = 1:n, R0 = R0v, gamma = gam, sigma = sig, rho = rho, prop_symp = psy,
                     immune_frac = imm, total_infections = tot_inf, total_reported = tot_rep,
                     peak_reported_wk = peak_rep_wk, attack_pct = attack,
                     finite = (seq_len(n) %in% ok)),
          file.path(DRAWS_DIR, "MAYV_lhs_draws.csv"), row.names = FALSE)
invisible(file.copy(file.path(DRAWS_DIR, "MAYV_lhs_draws.csv"),
                    file.path(DRAWS_DIR, sprintf("MAYV_lhs_draws_%s.csv", R0_SCENARIO)),
                    overwrite = TRUE))

mayv_lhs_ensemble <- list(
  rep = rep_mat[ok, , drop = FALSE], inf = inf_mat[ok, , drop = FALSE],
  R0 = R0v[ok], gamma = gam[ok], sigma = sig[ok], rho = rho[ok], prop_symp = psy[ok],
  immune_frac = imm[ok],  
  base_R0 = R0_median, R0_scenario = R0_SCENARIO, r0_is_peak = r0_is_peak,
  R0_median = R0_median, R0_sampled = TRUE, R0_lo = R0_LO, R0_hi = R0_HI,
  R0_meanlog = R0_meanlog, R0_sdlog = R0_sdlog,
  lat_meanlog = lat_ml, lat_sdlog = lat_sl, 
  lat_dist = lat$dist, lat_lo = lat_lo, lat_hi = lat_hi, MAYV_LATENT = MAYV_LATENT,
  ps_meanlog = ps_ml, ps_sdlog = ps_sl,
  base_gamma = g_m, base_sigma = 1/p_m, base_prop_symp = PSYMP_MED,
  base_immune = imm_base,
  season = season, N = N, A = A, age_df = age_df,
  I0_total = I0_total, E0 = E0, seed_week = seed_week,
  T_weeks = T_weeks, weeks = weeks, x_ticks = x_ticks, year_break = year_break)
saveRDS(mayv_lhs_ensemble, "outputs/MAYV_lhs_ensemble.rds")
saveRDS(mayv_lhs_ensemble, sprintf("outputs/MAYV_lhs_ensemble_%s.rds", R0_SCENARIO))

cat(sprintf("\nSaved MAYV_lhs_reported.png, MAYV_lhs_infections.png, MAYV_lhs_draws.csv,\n  MAYV_lhs_draws_%s.csv, MAYV_lhs_ensemble.rds and MAYV_lhs_ensemble_%s.rds (R0 median %.2f)\n",
            R0_SCENARIO, R0_SCENARIO, R0_median))
