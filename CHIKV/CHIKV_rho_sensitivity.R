# ============================================================
# CHIKV_rho_sensitivity.R -- post hoc reporting-rate sensitivity (Discussion).
# ------------------------------------------------------------
# The main analysis (CHIKV_lhs.R) samples the reporting rate rho ~ Beta(20, 60), the
# national estimate of 25% (16-35%). A post hoc comparison of long-term average reported
# cases in Caldas Novas (2016-2025) with long-term average model-predicted symptomatic
# cases gave a municipal reporting rate of 15.40% (95% UI 7.02-38.87%). This script asks
# what attack rate the fit implies if rho is instead fixed at that estimate's upper
# bound, 38.87%.
#
# LIKE-FOR-LIKE with the main analysis. It re-fits the SAME 1000 Latin hypercube draws:
# FOI, gamma, sigma and prop_symp are read per draw from the main run's
# lhs_draws/CHIKV_lhs_draws.csv, and the spline window, penalty, peak weighting,
# warm start and attack < 100% guard are all those of CHIKV_lhs.R. ONLY rho changes,
# so each draw is paired with its main-analysis counterpart and the difference between
# the two attack rates is due to the reporting rate alone.
#
# Run after: CHIKV_lhs.R  (which writes the draws file read here)
# Output:    CHIKV_rho_sensitivity.xlsx
# ============================================================
DEFS_ONLY <- TRUE
suppressMessages({library(writexl); source("CHIKV_lhs.R")})

RHO_UPPER <- 0.3887               # upper 95% bound of the municipal estimate 15.40% (7.02-38.87%)

draws_csv <- file.path(DRAWS_DIR, "CHIKV_lhs_draws.csv")
if (!file.exists(draws_csv))
  stop(draws_csv, " not found -- run CHIKV_lhs.R first (it writes the draws re-fitted here).")
D <- read.csv(draws_csv)
n <- nrow(D)

# Warm start from the point-estimate fit, exactly as CHIKV_lhs.R does.
warm <- refit(0.008, 0.54, 1/0.60, rho_pt, prop_symp, gen_start)$par

cat(sprintf("Re-fitting %d LHS draws with rho fixed at %.2f%%...\n", n, 100*RHO_UPPER))
attack <- infections <- totrep <- rep(NA_real_, n)
for (i in seq_len(n)) {
  f <- tryCatch(refit(D$FOI[i], D$gamma[i], D$sigma[i], RHO_UPPER, D$prop_symp[i], warm),
                error = function(e) NULL)
  if (is.null(f)) next
  attack[i] <- f$attack; infections[i] <- sum(f$inf); totrep[i] <- f$total
  if (i %% 100 == 0) cat("  ", i, "/", n, "\n")
}
# Same guard as the main analysis: drop a failed re-fit or an impossible attack rate.
ok_up   <- which(!is.na(totrep) & attack < 100)
ok_main <- which(D$feasible)

q3  <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
fmt <- function(x, d = 1) { q <- q3(x); sprintf("%s (%s-%s)",
        formatC(round(q[1],d), big.mark=",", format="f", digits=d),
        formatC(round(q[2],d), big.mark=",", format="f", digits=d),
        formatC(round(q[3],d), big.mark=",", format="f", digits=d)) }
both <- intersect(ok_main, ok_up)
summary_tbl <- data.frame(
  reporting_rate  = c("Main analysis: rho ~ Beta(20,60), 25% (16-35%)",
                      sprintf("Fixed at %.2f%% (upper bound of 15.40%%, 7.02-38.87%%)", 100*RHO_UPPER)),
  draws_retained  = c(sprintf("%d / %d", length(ok_main), n), sprintf("%d / %d", length(ok_up), n)),
  attack_rate_pct = c(fmt(D$attack_pct[ok_main]), fmt(attack[ok_up])),
  stringsAsFactors = FALSE)
print(summary_tbl, row.names = FALSE)
dq <- q3(attack[both] - D$attack_pct[both])
cat(sprintf("Paired change in attack rate (fixed %.2f%% minus main), %d draws retained in both: %.1f pp (%.1f to %.1f)\n",
            100*RHO_UPPER, length(both), dq[1], dq[2], dq[3]))

notes <- data.frame(item = c("Question", "Pairing", "Held as in the main analysis", "Guard",
                             "Attack rate"),
  detail = c(
  sprintf("Attack rate implied by the fit if the reporting rate is fixed at %.2f%%, the upper 95%% bound of the post hoc municipal estimate 15.40%% (7.02-38.87%%), instead of sampled from Beta(20,60).", 100*RHO_UPPER),
  "The same 1000 LHS draws as CHIKV_lhs.R (read from lhs_draws/CHIKV_lhs_draws.csv), re-fitted with rho replaced. Every other input is the draw's own value, so the two rows differ only in the reporting rate.",
  "FOI, gamma, sigma and prop_symp per draw; beta_t spline over weeks 1-49 held flat to week 52; ridge penalty; peak weighting; warm start from the point-estimate fit.",
  "A draw is dropped only if its re-fit fails or its attack rate reaches 100%, as in CHIKV_lhs.R. No goodness-of-fit filter.",
  "Infections over the 52-week window as a % of the susceptible pool at t = 0, median (95% UI) across retained draws."),
  stringsAsFactors = FALSE)
per_draw <- data.frame(draw = D$draw, FOI = D$FOI, gamma = D$gamma, sigma = D$sigma,
                       prop_symp = D$prop_symp,
                       rho_main = D$rho, attack_pct_main = D$attack_pct, retained_main = D$feasible,
                       rho_fixed = RHO_UPPER, attack_pct_fixed = attack, infections_fixed = infections,
                       total_reported_fixed = totrep, retained_fixed = seq_len(n) %in% ok_up)
write_xlsx(list(summary = summary_tbl, notes = notes, per_draw = per_draw),
           "CHIKV_rho_sensitivity.xlsx")
cat("\nWrote CHIKV_rho_sensitivity.xlsx\n")
