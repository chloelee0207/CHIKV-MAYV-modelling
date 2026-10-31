# ============================================================
# MAYV_engine.R
# Run order: source("MAYV_lhs.R")  (ensemble, once);  source("MAYV_engine.R")
# ============================================================
library(dplyr); library(ggplot2)
source("common/common.R")   # fmtq, load_burden_params, load_daly_params

# ---- knobs --------------------------------------------------
PHASE_MODE             <- "hosp_severity" 
OUTBREAK_ATTACK_THRESH <- 1.0 
MAYV_ZERO_DEATHS       <- TRUE              # Mayaro deaths set to 0 (CFR = 0)
if (!exists("DEFS_ONLY")) DEFS_ONLY <- FALSE   # TRUE = load the functions only (used by MAYV_owsa.R and MAYV_threshold_sensitivity.R)
set.seed(2031)

# ------------------------------------------------------------
# SEIRV
# ------------------------------------------------------------
seirv_vaccinated_MAYV <- function(
    T_weeks, A, N, R_init_prop, I0, base_beta, sigma, gamma, rho,
    target_age, total_coverage, weekly_delivery_speed, delay,
    VE_inf = 0, VE_block = 0, immun_delay = 2, prop_symp = 0.5242478,
    sub_steps = 7, E0 = rep(0, A), seed_week = 1) {
  pmax0 <- function(x) pmax(0, x); N_total <- sum(N); dt <- 1/sub_steps
  S <- E <- I <- R <- V <- matrix(0, A, T_weeks)
  V_covered <- vacc_delayed <- coverage_frac <- matrix(0, A, T_weeks)
  new_infections <- new_symptomatic <- matrix(0, A, T_weeks)
  target_idx <- which(target_age == 1); target_pop <- sum(N[target_idx])
  total_supply <- target_pop * total_coverage; weekly_dose_total <- total_supply * weekly_delivery_speed
  total_avail_age <- rep(0, A); total_avail_age[target_idx] <- total_supply * (N[target_idx]/target_pop)
  total_used_age <- rep(0, A); unvaccinated <- N
  S_now <- pmax0(N - E0 - R_init_prop*N); E_now <- E0; I_now <- rep(0, A); R_now <- R_init_prop*N; V_now <- rep(0, A)
  for (t in 1:T_weeks) {
    prev_Vc <- if (t > 1) V_covered[, t-1] else rep(0, A)
    if (t - immun_delay >= 1) {
      eff_dose <- vacc_delayed[, t-immun_delay]; immunized <- round(VE_inf*eff_dose)
      V_covered[, t] <- prev_Vc + eff_dose
    } else { immunized <- rep(0, A); V_covered[, t] <- prev_Vc }
    S_now <- pmax0(S_now - immunized); V_now <- V_now + immunized
    coverage_frac[, t] <- V_covered[, t]/N
    if (t >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(ceiling(weekly_dose_total*(N[a]/target_pop)), rem, unvaccinated[a],
                     total_avail_age[a]-total_used_age[a])
        if (alloc > 0) {
          prop_S <- if (N[a] > 0) S_now[a]/N[a] else 0; vacc_to_S <- round(alloc*prop_S)
          vacc_delayed[a, t] <- vacc_to_S; total_used_age[a] <- total_used_age[a]+alloc
          unvaccinated[a] <- unvaccinated[a]-alloc; rem <- rem-alloc
        }
      }
    }
    if (t == seed_week) { I_now <- I_now + I0; S_now <- pmax0(S_now - I0) }
    new_I_week <- rep(0, A); beta_t <- base_beta[t]
    for (k in 1:sub_steps) {
      foi <- beta_t*sum(I_now)/N_total
      new_E <- foi*S_now*dt; new_I <- sigma*E_now*dt; new_R <- gamma*I_now*dt
      S_now <- pmax0(S_now-new_E); E_now <- pmax0(E_now+new_E-new_I)
      I_now <- pmax0(I_now+new_I-new_R); R_now <- pmax0(R_now+new_R); new_I_week <- new_I_week+new_I
    }
    S[,t]<-S_now; E[,t]<-E_now; I[,t]<-I_now; R[,t]<-R_now; V[,t]<-V_now
    new_infections[,t] <- new_I_week
    new_symptomatic[,t] <- prop_symp*new_I_week*(1 - VE_block*coverage_frac[,t])
  }
  list(new_infections = new_infections, coverage_frac = coverage_frac,
       total_used_age = total_used_age,        # doses administered
       V_covered = V_covered)                  # cumulative doses that reached susceptibles
}

# ------------------------------------------------------------
# 0. Load the MAYV LHS ensemble
# ------------------------------------------------------------
if (!file.exists("outputs/MAYV_lhs_ensemble.rds"))
  stop("MAYV_lhs_ensemble.rds not found -- run MAYV_lhs.R first.")
E <- readRDS("outputs/MAYV_lhs_ensemble.rds")
N <- E$N; A <- E$A; age_df <- E$age_df; season <- E$season
seed_week <- E$seed_week; E0 <- E$E0; I0_total <- E$I0_total
T_weeks <- E$T_weeks; weeks <- E$weeks; x_ticks <- E$x_ticks; year_break <- E$year_break
n_ens <- length(E$R0); N_DRAWS <- n_ens
cat(sprintf("Loaded MAYV ensemble: %d draws | R0 scenario '%s' (peak=%s) | seed wk %d\n",
            n_ens, E$R0_scenario, E$r0_is_peak, seed_week))

EVAL_WIN <- 1:T_weeks

# ------------------------------------------------------------
# 1. Severity & DALY parameters
# ------------------------------------------------------------
invisible(list2env(load_burden_params(A), globalenv()))

{
  .dp <- as.data.frame(readxl::read_excel("data/disease_progression.xlsx",
                                          sheet = "disease_progression"))
  .r  <- .dp[.dp[[1]] == "Probability of hospitalisation among symptomatic cases (MAYV)", ]
  if (nrow(.r) != 1)
    stop("disease_progression.xlsx: expected exactly one MAYV hospitalisation row, found ", nrow(.r))
  hosp_a    <<- as.numeric(.r[[8]]); hosp_b <<- as.numeric(.r[[10]])
  hosp_rate <<- hosp_a / (hosp_a + hosp_b)
  cat(sprintf("MAYV hospitalisation override: Beta(%.3f, %.3f) -> mean %.4f (95%% UI %.4f-%.4f)\n",
              hosp_a, hosp_b, hosp_rate, qbeta(.025, hosp_a, hosp_b), qbeta(.975, hosp_a, hosp_b)))
}
dp <- load_daly_params()
age_weight <- rep(1, A) 
young_idx  <- which(age_to_band <= 4); old_idx <- which(age_to_band >= 5)

target_age <- rep(0, A); target_age[c(4,5,6,7,8)] <- 1     # eligible adults 18-59
target_pop_elig <- sum(N[target_age == 1])
immun_delay <- 2
start_pre   <- 17L      # 2025-W40, the same pre-outbreak campaign week as the CHIKV engine
stopifnot(start_pre >= 1, start_pre <= T_weeks, seed_week >= 1, seed_week <= T_weeks)

# ------------------------------------------------------------
# 2. Outcome
# ------------------------------------------------------------
OUTCOMES <- c("infections","reported","symptomatic","hospitalisations","deaths",
              "n_nonhosp","n_subacute","n_chronic",
              "yld_acute","yld_subacute","yld_chronic","yld","yll","daly","doses")
NNV_OUT  <- c("reported","symptomatic","hospitalisations","deaths","daly")

outcome_one <- function(symp_age, infections, doses, rho, hosp_j, cfr_j, le_band,
                        dwmm, dwsv, dwch, dumm, dusv, dusb, duch, acy, aco, sby, sbo, chy, cho) {
  symp_w  <- symp_age * age_weight
  symp_dw <- if (sum(symp_w) > 0) symp_w * (sum(symp_age)/sum(symp_w)) else symp_age
  st <- sum(symp_dw); sy <- sum(symp_dw[young_idx]); so <- sum(symp_dw[old_idx])
  reported <- rho * st                                     # rho * symptomatic
  n_hosp    <- st * hosp_j                                 # hospitalised (severe acute)
  n_nonhosp <- st * (1 - hosp_j)                           # non-hospitalised (mild/moderate acute)
  # Acute / sub-acute / chronic shares come from the Halsey 2013 Dirichlet(1, 4, 11) draw
  # dividing by their sum keeps every symptomatic case in exactly one phase
  sy_t <- acy + sby + chy; so_t <- aco + sbo + cho
  n_acute    <- sy*acy/sy_t + so*aco/so_t
  n_subacute <- sy*sby/sy_t + so*sbo/so_t
  n_chronic  <- sy*chy/sy_t + so*cho/so_t
  f_acute    <- if (st > 0) n_acute / st else 0
  if (PHASE_MODE != "hosp_severity")
    stop("MAYV engine implements PHASE_MODE 'hosp_severity' only (the CHIKV default).")
  yld_ac  <- (n_nonhosp*dwmm*dumm + n_hosp*dwsv*dusv) * f_acute
  yld_sub <- n_subacute * dwch * dusb
  yld_chr <- n_chronic  * dwch * duch
  yld_tot <- unname(yld_ac + yld_sub + yld_chr)
  deaths  <- sum(symp_dw * cfr_j)
  yll     <- sum(symp_dw * cfr_j * le_band[age_to_band])
  c(infections = infections, reported = reported, symptomatic = st,
    hospitalisations = n_hosp, deaths = deaths,
    n_nonhosp = n_nonhosp, n_subacute = n_subacute, n_chronic = n_chronic,
    yld_acute = unname(yld_ac), yld_subacute = unname(yld_sub),
    yld_chronic = unname(yld_chr), yld = yld_tot,
    yll = yll, daly = yld_tot + yll, doses = doses)
}

if (!isTRUE(DEFS_ONLY)) {

# ------------------------------------------------------------
# 3. Latin-hypercube
# ------------------------------------------------------------
beta_from_ci <- function(m, lo, hi) { v <- ((hi-lo)/(2*1.96))^2; k <- m*(1-m)/v-1; c(a=m*k, b=(1-m)*k) }
cov_ab <- beta_from_ci(0.30, 0.20, 0.40) 
veb_ab <- c(a = 4, b = 8)                   # mean 0.3333 (4/12), 95% UI 0.109-0.610
del_ab <- beta_from_ci(0.10, 0.09, 0.11)     # weekly delivery speed

lhs_col <- function(n) (sample.int(n) - runif(n)) / n
K <- 42; U <- sapply(1:K, function(j) lhs_col(N_DRAWS)); col <- 0 
nextU <- function(w = 1) { idx <- (col+1):(col+w); col <<- col + w; U[, idx, drop = FALSE] }
cov_d  <- qbeta(nextU(), cov_ab["a"], cov_ab["b"])
veb_d  <- qbeta(nextU(), veb_ab["a"], veb_ab["b"])
del_d  <- qbeta(nextU(), del_ab["a"], del_ab["b"])
# Deployment delay, 1-3 weeks (median 2), as in the CHIKV engine
delay_d <- 1 + round(2 * nextU())
hosp_d <- qbeta(nextU(), hosp_a, hosp_b)
cfrH_d <- qbeta(nextU(9), matrix(cfr_hosp_a, N_DRAWS, 9, byrow=TRUE), matrix(cfr_hosp_b, N_DRAWS, 9, byrow=TRUE))
cfrN_d <- qbeta(nextU(9), matrix(cfr_nonh_a, N_DRAWS, 9, byrow=TRUE), matrix(cfr_nonh_b, N_DRAWS, 9, byrow=TRUE))
dwMM_d <- qbeta(nextU(), dp$dw_mm$a, dp$dw_mm$b); dwSV_d <- qbeta(nextU(), dp$dw_sev$a, dp$dw_sev$b)
dwCH_d <- qbeta(nextU(), dp$dw_chr$a, dp$dw_chr$b)
{
  .d2 <- as.data.frame(readxl::read_excel("data/disease_progression.xlsx",
                                          sheet = "disease_progression"))
  .r2 <- .d2[.d2[[1]] == "Duration of illness for mild and moderate Mayaro (years)", ]
  if (nrow(.r2) != 1)
    stop("disease_progression.xlsx: expected exactly one MAYV mild/moderate duration row, found ", nrow(.r2))
  MAYV_ACUTE_DUR <<- list(m = as.numeric(.r2[[8]]), s = as.numeric(.r2[[10]]))
  .gs <- as.numeric(as.data.frame(readxl::read_excel("data/model_calibration.xlsx", 1)) |>
           subset(Group == "MAYV" & grepl("gamma", Parameter)) |> (\(z) z[["Value 2"]])())
  if (!isTRUE(all.equal(MAYV_ACUTE_DUR$s, .gs, tolerance = 1e-8)))
    stop("MAYV acute duration sdlog (", MAYV_ACUTE_DUR$s, ") does not match gamma sdlog (", .gs,
         "). They are the same quantity -- fix the workbooks before running.")
  cat(sprintf("MAYV acute duration from workbook: %.2f d (95%% UI %.2f-%.2f), sdlog matches gamma\n",
              365*exp(MAYV_ACUTE_DUR$m), 365*qlnorm(.025, MAYV_ACUTE_DUR$m, MAYV_ACUTE_DUR$s),
              365*qlnorm(.975, MAYV_ACUTE_DUR$m, MAYV_ACUTE_DUR$s)))
}
duMM_d <- qlnorm(nextU(), MAYV_ACUTE_DUR$m, MAYV_ACUTE_DUR$s)
duSV_d <- qlnorm(nextU(), dp$du_sev$m, dp$du_sev$s)   # severe, borrowed from CHIKV
duSB_d <- qlnorm(nextU(), dp$du_sub$m, dp$du_sub$s)   # sub-acute duration (chronic DW applied)
duCH_d <- qlnorm(nextU(), dp$du_chr$m, dp$du_chr$s)
le_d   <- qlnorm(nextU(9), matrix(dp$le$m, N_DRAWS, 9, byrow=TRUE), matrix(dp$le$s, N_DRAWS, 9, byrow=TRUE))

HALSEY_COUNTS <- c(acute = 1, subacute = 4, chronic = 11)
gam <- sapply(HALSEY_COUNTS, function(a) qgamma(nextU(), shape = a, rate = 1))
dir_d  <- gam / rowSums(gam)                                   # Dirichlet(1, 4, 11)
acy_d <- aco_d <- dir_d[, 1]
sby_d <- sbo_d <- dir_d[, 2]
chy_d <- cho_d <- dir_d[, 3]
stopifnot(col == K)

if (MAYV_ZERO_DEATHS) { cfrH_d[] <- 0; cfrN_d[] <- 0; cfr_vec[] <- 0 }

# ------------------------------------------------------------
# 4. Monte Carlo
# ------------------------------------------------------------
scen_names <- c("No vaccine (baseline)", "Pre-outbreak | Disease-blocking")
vac_name   <- "Pre-outbreak | Disease-blocking"
per_draw <- setNames(lapply(scen_names, function(x)
  matrix(NA_real_, N_DRAWS, length(OUTCOMES), dimnames = list(NULL, OUTCOMES))), scen_names)
attack_base <- numeric(N_DRAWS)
# Susceptible pool at t = 0, per draw
sus_pool <- numeric(N_DRAWS)
wk_base <- wk_vacc <- matrix(NA_real_, N_DRAWS, T_weeks)   # weekly symptomatic
doses_deliv <- doses_ontarget <- numeric(N_DRAWS)

cat(sprintf("Running %d draws (baseline + pre-outbreak disease-blocking), R0 ~ %.1f-%.1f...\n",
            N_DRAWS, E$R0_lo, E$R0_hi))
for (i in 1:N_DRAWS) {
  R0i<-E$R0[i]; gi<-E$gamma[i]; si<-E$sigma[i]; ri<-E$rho[i]; psi<-E$prop_symp[i]; immi<-E$immune_frac[i]
  Rimm <- rep(immi, A); sus <- N*(1-Rimm); I0i <- I0_total * sus/sum(sus)
  sus_pool[i] <- sum(sus)
  base_beta <- R0i * gi * season
  run <- seirv_vaccinated_MAYV(T_weeks, A, N, Rimm, I0i, base_beta, si, gi, ri,
           target_age, cov_d[i], del_d[i], start_pre + delay_d[i], VE_inf = 0, VE_block = veb_d[i],
           immun_delay = immun_delay, prop_symp = psi, E0 = E0, seed_week = seed_week)
  ninf <- run$new_infections; covf <- run$coverage_frac
  doses_deliv[i]    <- sum(run$total_used_age)              # actually administered
  doses_ontarget[i] <- sum(run$V_covered[, T_weeks])       # reached susceptibles
  inf_tot <- sum(ninf[, EVAL_WIN, drop = FALSE])
  attack_base[i] <- 100 * inf_tot / sum(sus)

  symp_base <- rowSums((psi * ninf)[, EVAL_WIN, drop = FALSE])                      # no vaccine
  symp_vacc <- rowSums((psi * ninf * (1 - veb_d[i]*covf))[, EVAL_WIN, drop = FALSE]) # disease-blocked
  wk_base[i, ] <- colSums(psi * ninf)                                              # weekly
  wk_vacc[i, ] <- colSums(psi * ninf * (1 - veb_d[i]*covf))
  cfr_j <- (hosp_d[i]*cfrH_d[i,] + (1-hosp_d[i])*cfrN_d[i,])[age_to_band]
  args_daly <- list(le_d[i,], dwMM_d[i], dwSV_d[i], dwCH_d[i], duMM_d[i], duSV_d[i], duSB_d[i],
                    duCH_d[i], acy_d[i], aco_d[i], sby_d[i], sbo_d[i], chy_d[i], cho_d[i])
  per_draw[["No vaccine (baseline)"]][i, ] <-
    do.call(outcome_one, c(list(symp_base, inf_tot, 0, ri, hosp_d[i], cfr_j), args_daly))
  per_draw[[vac_name]][i, ] <-
    do.call(outcome_one, c(list(symp_vacc, inf_tot, target_pop_elig*cov_d[i], ri, hosp_d[i], cfr_j), args_daly))
  if (i %% 100 == 0) cat("  ", i, "/", N_DRAWS, "\n")
}

# ------------------------------------------------------------
# 5. averted + NNV per draw
# ------------------------------------------------------------
# `outbreak` is the row index used downstream (MAYV_outputs.R reads G$outbreak)
outbreak   <- seq_len(N_DRAWS)
p_outbreak <- 1
frac_over_thresh <- mean(attack_base > OUTBREAK_ATTACK_THRESH)
cat(sprintf("\nR0 SAMPLED %.1f-%.1f (%s scenario, median %.2f) -> reporting over ALL %d draws, no conditioning.\n",
            E$R0_lo, E$R0_hi, E$R0_scenario, E$R0_median, N_DRAWS))
cat(sprintf("  Diagnostic: %.1f%% of draws exceed the legacy %.1f%% attack-rate filter (NOT used to condition).\n",
            100*frac_over_thresh, OUTBREAK_ATTACK_THRESH))
cat(sprintf("  Baseline attack rate: median %.3f%% [%.3f%%, %.3f%%]\n",
            median(attack_base), quantile(attack_base, .025), quantile(attack_base, .975)))
cat("  Baseline symptomatic deciles: ",
    paste(round(quantile(per_draw[["No vaccine (baseline)"]][, "symptomatic"], seq(0, 1, .1))),
          collapse = " | "), "\n")

base_pd <- per_draw[["No vaccine (baseline)"]]
vac_pd  <- per_draw[[vac_name]]
averted <- base_pd[, NNV_OUT, drop=FALSE] - vac_pd[, NNV_OUT, drop=FALSE]
nnv     <- vac_pd[, "doses"] / averted 
nnv[!is.finite(nnv) | nnv < 0] <- NA

q3  <- function(x) quantile(x, c(.5,.025,.975), na.rm = TRUE)
aggc <- function(mat) do.call(rbind, lapply(colnames(mat), function(o) {
  q <- q3(mat[outbreak, o]); data.frame(outcome=o, median=q[1], lo=q[2], hi=q[3], row.names=NULL) }))
agg_burden_cond  <- cbind(scenario = "baseline",           aggc(base_pd))
agg_averted_cond <- cbind(scenario = vac_name,             aggc(averted))
agg_nnv_cond     <- cbind(scenario = vac_name,             aggc(nnv))

# ------------------------------------------------------------
# 6. summary
# ------------------------------------------------------------
fmtq_ <- function(v, d=0) fmtq(v[outbreak], d)
cat("\n=== Baseline burden (median, 95% UI over all draws) ===\n")
for (o in c("infections","symptomatic","hospitalisations","deaths","daly"))
  cat(sprintf("  %-16s %s\n", o, fmtq_(base_pd[, o], if (o=="deaths") 1 else 0)))
cat("\n=== Pre-outbreak disease-blocking vaccine, AVERTED (all draws) ===\n")
for (o in NNV_OUT) cat(sprintf("  %-16s %s\n", o, fmtq_(averted[, o], if (o=="deaths") 1 else 0)))
cat("\n=== NNV (doses per burden averted, all draws) ===\n")
for (o in NNV_OUT) cat(sprintf("  %-16s %s\n", o, fmtq(nnv[outbreak, o], 0)))

# ------------------------------------------------------------
# 6b. Epidemic curve
# ------------------------------------------------------------
draw_set <- seq_len(N_DRAWS)
bandq    <- function(M) apply(M[draw_set, , drop = FALSE], 2, quantile, c(.025,.5,.975), na.rm = TRUE)
bb <- bandq(wk_base); bv <- bandq(wk_vacc)
wk_num   <- function(idx) ifelse(idx <= 30, idx + 23, idx - 30)
tick_idx <- c(7, 17, 27, 40, 50); xt <- data.frame(i = tick_idx, w = wk_num(tick_idx))
red <- 100 * (base_pd[draw_set,"symptomatic"] - vac_pd[draw_set,"symptomatic"]) / base_pd[draw_set,"symptomatic"]
rq  <- quantile(red, c(.5,.025,.975), na.rm = TRUE)
lab <- sprintf("%% symptomatic reduction\nDisease-blocking: %.1f%% (%.1f-%.1f%%)", rq[1], rq[2], rq[3])
roll_beg <- start_pre + median(delay_d)                            # vaccination begins after the delay
roll_end <- roll_beg + max(1, round(1 / mean(del_d))) - 0.5
pdf_df <- data.frame(week = 1:T_weeks, b_lo=bb[1,], b_md=bb[2,], b_hi=bb[3,],
                     v_lo=bv[1,], v_md=bv[2,], v_hi=bv[3,])
ttl <- sprintf("MAYV symptomatic cases (2025-W24 - 2026-W22) | %s scenario, R0 %.1f-%.1f (all %d draws)",
               E$R0_scenario, E$R0_lo, E$R0_hi, N_DRAWS)
p_epi <- ggplot(pdf_df, aes(week)) +
  annotate("rect", xmin = roll_beg-0.5, xmax = roll_end, ymin = -Inf, ymax = Inf, fill = "#cdebc5", alpha = 0.5) +
  geom_vline(xintercept = E$year_break, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = b_lo, ymax = b_hi), fill = "grey55",  alpha = 0.25) +
  geom_ribbon(aes(ymin = v_lo, ymax = v_hi), fill = "#4292c6", alpha = 0.22) +
  geom_line(aes(y = b_md), colour = "grey30",  linewidth = 1) +
  geom_line(aes(y = v_md), colour = "#2171b5", linewidth = 1) +
  annotate("text", x = T_weeks, y = Inf, label = lab, hjust = 1, vjust = 1.2, size = 3.1) +
  scale_x_continuous(breaks = xt$i, labels = xt$w) +
  labs(x = "Week", y = "Predicted symptomatic cases", title = ttl) +
  theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
                       panel.grid.minor = element_blank())
epi_fn <- sprintf("outputs/MAYV_symptomatic_%s.png", E$R0_scenario)
ggsave(epi_fn, p_epi, width = 8, height = 4.5, dpi = 120)
cat(sprintf("Saved epidemic-curve plot: %s (grey = no vaccine, blue = disease-blocking)\n", epi_fn))

# no vaccine
base_curve_plot <- function(M, draws, ytitle, ttl, fn) {
  bq <- apply(M[draws, , drop = FALSE], 2, quantile, c(.025,.5,.975), na.rm = TRUE)
  dfp <- data.frame(week = 1:T_weeks, lo = bq[1,], md = bq[2,], hi = bq[3,])
  p <- ggplot(dfp, aes(week)) +
    geom_vline(xintercept = E$year_break, linetype = "dashed", colour = "grey55") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey55", alpha = 0.28) +
    geom_line(aes(y = md), colour = "grey20", linewidth = 1) +
    scale_x_continuous(breaks = xt$i, labels = xt$w) +
    labs(x = "Week", y = ytitle, title = ttl) +
    theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
                         panel.grid.minor = element_blank())
  ggsave(fn, p, width = 8, height = 4.5, dpi = 120); cat("Saved baseline (no-vaccine) plot:", fn, "\n")
}
base_curve_plot(wk_base, draw_set, "Predicted symptomatic cases (no vaccine)",
  sprintf("MAYV symptomatic, NO vaccine (2025-W24 - 2026-W22) | %s scenario, R0 %.1f-%.1f (all %d draws)",
          E$R0_scenario, E$R0_lo, E$R0_hi, N_DRAWS),
  sprintf("outputs/MAYV_baseline_%s.png", E$R0_scenario))


# ------------------------------------------------------------
# 7. save results
# ------------------------------------------------------------
{
  mcols <- function(m, stem) setNames(as.data.frame(m), sprintf("%s_%d", stem, seq_len(ncol(m))))
  v <- as.vector
  burden_draws <- cbind(
    data.frame(draw = seq_len(N_DRAWS), transmission_row = seq_len(N_DRAWS),
               R0_scenario = E$R0_scenario,
               coverage = v(cov_d), ve_cross_protection = v(veb_d), delivery = v(del_d),
               delay_weeks = v(delay_d), hosp_rate = v(hosp_d)),
    mcols(cfrH_d, "cfr_hosp_band"), mcols(cfrN_d, "cfr_nonhosp_band"),
    data.frame(dw_mild_mod = v(dwMM_d), dw_severe = v(dwSV_d), dw_chronic = v(dwCH_d),
               dur_mild_mod = v(duMM_d), dur_severe = v(duSV_d),
               dur_subacute = v(duSB_d), dur_chronic = v(duCH_d)),
    mcols(le_d, "life_expectancy_band"),
    data.frame(rec_acute = v(acy_d), rec_subacute = v(sby_d), rec_chronic = v(chy_d)))
  dir.create(DRAWS_DIR, recursive = TRUE, showWarnings = FALSE)
  bd_path <- file.path(DRAWS_DIR, "MAYV_burden_draws.csv")
  write.csv(burden_draws, bd_path, row.names = FALSE)
  file.copy(bd_path, file.path(DRAWS_DIR, sprintf("MAYV_burden_draws_%s.csv", E$R0_scenario)),
            overwrite = TRUE)
  cat(sprintf("Saved %s (%d draws x %d parameters, scenario '%s').\n",
              bd_path, nrow(burden_draws), ncol(burden_draws) - 3, E$R0_scenario))
}

saveRDS(list(
  per_draw = per_draw, averted = averted, nnv = nnv,
  attack_base = attack_base, outbreak = outbreak, p_outbreak = p_outbreak,
  sus_pool = sus_pool, pop_total = sum(N),
  R0_median = E$R0_median, R0_sampled = TRUE, R0_lo = E$R0_lo, R0_hi = E$R0_hi,
  conditioning = "none (all draws; R0 sampled within the scenario range)",
  frac_over_legacy_thresh = frac_over_thresh,
  OUTBREAK_ATTACK_THRESH = OUTBREAK_ATTACK_THRESH,
  rho_draw = E$rho,
  doses_deliv = doses_deliv, doses_ontarget = doses_ontarget,
  agg_burden_cond = agg_burden_cond, agg_averted_cond = agg_averted_cond, agg_nnv_cond = agg_nnv_cond,
  scen_names = scen_names, vac_name = vac_name, OUTCOMES = OUTCOMES, NNV_OUT = NNV_OUT,
  N_DRAWS = N_DRAWS, PHASE_MODE = PHASE_MODE, R0_scenario = E$R0_scenario,
  target_pop_elig = target_pop_elig, cov_d = cov_d, veb_d = veb_d, del_d = del_d,
  EVAL_WIN = EVAL_WIN, 
  wk_base = wk_base, wk_vacc = wk_vacc, T_weeks = T_weeks,    # weekly symptomatic
  start_pre = start_pre, immun_delay = immun_delay, deliv_median = median(del_d),
  delay_d = delay_d, dose_start = start_pre + median(delay_d)),
  "outputs/MAYV_engine_results.rds")
file.copy("outputs/MAYV_engine_results.rds",
          sprintf("outputs/MAYV_engine_results_%s.rds", E$R0_scenario), overwrite = TRUE)
cat(sprintf("\nSaved MAYV_engine_results.rds and MAYV_engine_results_%s.rds\n  (per-draw + aggregates over ALL draws, R0 %.1f-%.1f; severity-phase counts included).\n",
            E$R0_scenario, E$R0_lo, E$R0_hi))


}  # end !DEFS_ONLY
