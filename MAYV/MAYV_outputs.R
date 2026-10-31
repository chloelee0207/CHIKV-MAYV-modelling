# ============================================================
# MAYV_outputs.R
# Run order: source("MAYV_lhs.R"); source("MAYV_engine.R"); source(this)
# ============================================================
library(writexl)
if (!exists("fmtq")) source("common/common.R") 

if (!file.exists("outputs/MAYV_engine_results.rds"))
  stop("MAYV_engine_results.rds not found -- run MAYV_engine.R first.")
G <- readRDS("outputs/MAYV_engine_results.rds")

ok        <- G$outbreak 
rho_draw  <- G$rho_draw
scen_names <- G$scen_names
vac_names  <- setdiff(scen_names, "No vaccine (baseline)")
base_pd   <- G$per_draw[["No vaccine (baseline)"]]
outcomes  <- c("infections","symptomatic","hospitalisations","deaths") 

lab_timing <- function(nm) if (nm == "No vaccine (baseline)") "No vaccine" else sub(" \\|.*","",nm)
lab_arm    <- function(nm) if (nm == "No vaccine (baseline)") "No vaccine" else sub(".*\\| ","",nm)
d_death    <- function(o) if (o == "deaths") 1 else 0

# % reduction
fmtpct <- function(ratio) {
  q <- quantile(ratio, c(.5, .025, .975), na.rm = TRUE)
  sprintf("%.1f%% (%.1f - %.1f%%)", 100*q[1], 100*q[2], 100*q[3])
}

# averted (baseline - scenario)
av <- setNames(lapply(vac_names, function(nm)
  base_pd[ok, outcomes, drop = FALSE] - G$per_draw[[nm]][ok, outcomes, drop = FALSE]), vac_names)
base_symp_ok <- base_pd[ok, "symptomatic"]

# ============================================================
# MAYV_vacc_outputs.xlsx
# ============================================================
# baseline true-vs-reported: reported = rho x true
base_tbl <- do.call(rbind, lapply(outcomes, function(o) {
  tv <- base_pd[ok, o]; d <- d_death(o)
  data.frame(outcome = o, true = fmtq(tv, d), reported = fmtq(rho_draw[ok]*tv, d),
             row.names = NULL)
}))

# per-scenario totals, true & reported
vtr <- do.call(rbind, lapply(scen_names, function(nm) {
  m <- G$per_draw[[nm]]
  data.frame(timing = lab_timing(nm), arm = lab_arm(nm), outcome = outcomes,
             true     = sapply(outcomes, function(o) fmtq(m[ok,o],             d_death(o))),
             reported = sapply(outcomes, function(o) fmtq(rho_draw[ok]*m[ok,o], d_death(o))),
             row.names = NULL)
}))

# averted
mc_tbl <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av[[nm]]
  data.frame(timing = lab_timing(nm), arm = lab_arm(nm),
             Infections = fmtq(m[,"infections"]), Symptomatic = fmtq(m[,"symptomatic"]),
             Hospitalisations = fmtq(m[,"hospitalisations"], 1), Deaths = fmtq(m[,"deaths"], 2),
             pct_symp = fmtpct(m[,"symptomatic"] / base_symp_ok),
             row.names = NULL)
}))

# outcomes averted per 100,000 doses
per100k_outcomes <- c(infections="Infections", symptomatic="Symptomatic",
                      hospitalisations="Hospitalisations", deaths="Deaths",
                      yld="YLDs", yll="YLLs", daly="DALYs")
mc_per100k <- do.call(rbind, lapply(vac_names, function(nm) {
  doses <- G$per_draw[[nm]][ok, "doses"]
  b <- base_pd[ok, , drop = FALSE]; s <- G$per_draw[[nm]][ok, , drop = FALSE]
  cells <- lapply(names(per100k_outcomes), function(o) {
    r <- 1e5 * (b[,o] - s[,o]) / doses; r[!is.finite(r) | r < 0] <- NA; fmtq(r, 0)
  })
  setNames(data.frame(timing = lab_timing(nm), arm = lab_arm(nm), cells,
                      row.names = NULL, check.names = FALSE),
           c("timing","arm", unname(per100k_outcomes)))
}))

scenario_totals <- data.frame(
  scenario = scen_names,
  do.call(rbind, lapply(scen_names, function(nm)
    round(apply(G$per_draw[[nm]][ok, outcomes, drop=FALSE], 2, median, na.rm=TRUE), 1))),
  row.names = NULL, check.names = FALSE)


# doses actually delivered to the eligible 18-59, and dose wastage
doses_wastage <- data.frame(
  timing = lab_timing(G$vac_name), arm = lab_arm(G$vac_name),
  `eligible population (18-59)` = round(G$target_pop_elig),
  `coverage of eligible`        = fmtpct(G$cov_d),
  `doses delivered`             = fmtq(G$doses_deliv, 0),
  `doses on-target (to susceptibles)` = fmtq(G$doses_ontarget, 0),
  `wastage %`                   = fmtpct(1 - G$doses_ontarget / G$doses_deliv),
  check.names = FALSE, row.names = NULL)


# ------------------------------------------------------------
# baseline attack rates
# Denominator is the per-draw susceptible pool at t = 0
# ------------------------------------------------------------
stopifnot(!is.null(G$sus_pool)) 
.sus <- G$sus_pool[ok]
attack_rates <- data.frame(
  measure = c("Baseline immunity (% of population)",
              "Infections / susceptibles (%)", "Symptomatic / susceptibles (%)",
              "Susceptible pool at t = 0", "Total population"),
  baseline = c(
    fmtq(100 * (1 - .sus / G$pop_total), 1),
    fmtq(100 * base_pd[ok, "infections"]  / .sus, 1),
    fmtq(100 * base_pd[ok, "symptomatic"] / .sus, 1),
    fmtq(.sus, 0),
    format(round(G$pop_total), big.mark = ",")),
  row.names = NULL)

write_xlsx(list(baseline_true_reported = base_tbl,
                vaccinated_true_reported = vtr, averted_MC_95UI = mc_tbl,
                averted_per_100k_doses = mc_per100k, doses_wastage = doses_wastage,
                scenario_totals = scenario_totals, attack_rates = attack_rates),
           "outputs/MAYV_vacc_outputs.xlsx")

# ============================================================
# DALY workbook
# ============================================================
daly_by_scenario <- do.call(rbind, lapply(scen_names, function(nm) {
  m <- G$per_draw[[nm]]
  data.frame(timing = lab_timing(nm), arm = lab_arm(nm),
             YLD = fmtq(m[ok,"yld"]), YLL = fmtq(m[ok,"yll"]), DALY = fmtq(m[ok,"daly"]),
             YLD_acute = fmtq(m[ok,"yld_acute"]), YLD_subacute = fmtq(m[ok,"yld_subacute"]),
             YLD_chronic = fmtq(m[ok,"yld_chronic"]), row.names = NULL, check.names = FALSE)
}))
base_daly_ok <- base_pd[ok, "daly"]
# YLD and YLL averted, absolute and per 100,000 doses
daly_averted <- do.call(rbind, lapply(vac_names, function(nm) {
  s <- G$per_draw[[nm]][ok, , drop = FALSE]; doses <- s[,"doses"]
  ad <- base_daly_ok - s[,"daly"]
  a_yld <- base_pd[ok,"yld"] - s[,"yld"]; a_yll <- base_pd[ok,"yll"] - s[,"yll"]
  p100 <- function(v) { r <- 1e5 * v / doses; r[!is.finite(r) | r < 0] <- NA; fmtq(r, 0) }
  data.frame(timing = lab_timing(nm), arm = lab_arm(nm),
             YLD_averted = fmtq(a_yld), YLL_averted = fmtq(a_yll),
             DALY_averted = fmtq(ad),
             pct_DALY = fmtpct(ad / base_daly_ok),
             YLD_averted_per_100k_doses  = p100(a_yld),
             YLL_averted_per_100k_doses  = p100(a_yll),
             DALY_averted_per_100k_doses = p100(ad),
             row.names = NULL, check.names = FALSE)
}))

# ------------------------------------------------------------
# yld by phase
# ------------------------------------------------------------
PHASE <- c(yld_acute = "Acute (<=14 d)", yld_subacute = "Sub-acute (14 d - 3 mo)",
           yld_chronic = "Chronic (>3 mo)", yld = "TOTAL YLD")
yld_by_phase <- do.call(rbind, lapply(vac_names, function(nm) {
  s_ <- G$per_draw[[nm]][ok, , drop = FALSE]; doses <- s_[, "doses"]
  tot <- base_pd[ok, "yld"]
  do.call(rbind, lapply(names(PHASE), function(k) {
    bl <- base_pd[ok, k]; av <- bl - s_[, k]
    r  <- 1e5 * av / doses; r[!is.finite(r) | r < 0] <- NA
    data.frame(timing = lab_timing(nm), arm = lab_arm(nm), phase = PHASE[[k]],
               YLD_baseline = fmtq(bl), YLD_vaccinated = fmtq(s_[, k]),
               YLD_averted = fmtq(av), YLD_averted_per_100k_doses = fmtq(r, 1),
               pct_of_total_YLD = sprintf("%.1f%%", 100 * median(bl) / median(tot)),
               row.names = NULL, check.names = FALSE) }))
}))

write_xlsx(list(daly_by_scenario = daly_by_scenario, daly_averted = daly_averted,
                yld_by_phase = yld_by_phase),
           "outputs/MAYV_daly_outputs.xlsx")

# ============================================================
# NNV workbook
# ============================================================
# NNV = doses / burden averted
nnv_outcomes <- c(infections = "Infection", symptomatic = "Symptomatic case",
                  hospitalisations = "Hospitalisation", deaths = "Death", daly = "DALY")
nnv_tbl <- do.call(rbind, lapply(vac_names, function(nm) {
  doses <- G$per_draw[[nm]][ok, "doses"]
  cells <- lapply(names(nnv_outcomes), function(o) {
    avo <- base_pd[ok, o] - G$per_draw[[nm]][ok, o]
    n   <- doses / avo; n[!is.finite(n) | n < 0] <- NA
    fmtq(n, 0)
  })
  setNames(data.frame(timing = lab_timing(nm), arm = lab_arm(nm), cells,
                      row.names = NULL, check.names = FALSE),
           c("timing","arm", unname(nnv_outcomes)))
}))
write_xlsx(list(nnv = nnv_tbl), "outputs/MAYV_nnv_outputs.xlsx")

# ------------------------------------------------------------
# Residual burden
# ------------------------------------------------------------
library(ggplot2)
stopifnot(max(base_pd[ok, "deaths"]) == 0)
stopifnot(length(vac_names) == 1)

pct_of_base <- function(num, den) { r <- 100 * num / den; r[!is.finite(r)] <- NA; r }
res_rows <- list()
add_res <- function(outcome, scen, v) res_rows[[length(res_rows)+1]] <<- data.frame(
  outcome = outcome, scenario = scen,
  med = median(v, na.rm = TRUE), lo = quantile(v, .025, na.rm = TRUE),
  hi = quantile(v, .975, na.rm = TRUE), row.names = NULL)

# Outcomes shown in the residual-burden panels (Fig. 1 B-E)
resid_outcomes <- list(
  list(lab = "Cumulative symptomatic cases", col = "symptomatic"),
  list(lab = "Cumulative hospitalisations",  col = "hospitalisations"),
  list(lab = "Cumulative DALYs",             col = "daly"))
PLOT_OUTCOMES <- sapply(resid_outcomes, `[[`, "lab")
nm <- vac_names[1]
for (o in resid_outcomes) {
  v <- pct_of_base(G$per_draw[[nm]][ok, o$col], base_pd[ok, o$col])
  add_res(o$lab, "No vaccination", 100)
  add_res(o$lab, "Vaccination",    v)
}
resid <- do.call(rbind, res_rows)
resid$scenario <- factor(resid$scenario, levels = c("No vaccination", "Vaccination"))
resid$outcome  <- factor(resid$outcome, levels = sapply(resid_outcomes, `[[`, "lab"))
resid$arm      <- factor("Disease-blocking")      # MAYV has no infection-blocking arm
resid_plot <- resid[resid$outcome %in% PLOT_OUTCOMES, ]
resid_plot$outcome <- droplevels(resid_plot$outcome)

p_resid <- ggplot(resid_plot, aes(scenario, med, fill = scenario)) +
  geom_col(width = .6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .35) +
  facet_grid(outcome ~ arm) +
  scale_fill_manual(values = c("No vaccination" = "grey60", "Vaccination" = "#4e79a7"),
                    name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA)) +
  labs(x = NULL, y = "Cumulative burden (% of no vaccination)") +
  theme_bw(11) +
  theme(axis.text.x = element_text(size = 9), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank())
print(p_resid)
ggsave("outputs/MAYV_residual_burden.png", p_resid, width = 3.6, height = 8.0, dpi = 130)
# % reduction from no vaccination
resid$red_med <- 100 - resid$med
resid$red_lo  <- 100 - resid$hi
resid$red_hi  <- 100 - resid$lo
write_xlsx(list(residual_burden_pct = resid), "outputs/MAYV_residual_burden.xlsx")
attr(resid_plot, "R0_scenario") <- G$R0_scenario
saveRDS(resid_plot, "outputs/MAYV_residual_burden.rds") 
cat("Saved MAYV_residual_burden.png and .xlsx (burden as % of no vaccination;\n",
    "     deaths panel omitted -- MAYV CFR is fixed at 0).\n", sep = "")

# ------------------------------------------------------------
cat("Wrote MAYV_vacc_outputs.xlsx (baseline_true_reported, vaccinated_true_reported,\n",
    "     averted_MC_95UI, averted_per_100k_doses, scenario_totals)\n", sep = "")
cat("Wrote MAYV_daly_outputs.xlsx (daly_by_scenario, daly_averted, yld_by_phase)\n")
cat("Wrote MAYV_nnv_outputs.xlsx  (nnv)\n\n")
cat(sprintf("R0 %.1f-%.1f (median %.2f), all %d draws (no conditioning).  pct symptomatic reduced: %s\n",
            G$R0_lo, G$R0_hi, G$R0_median, G$N_DRAWS, mc_tbl$pct_symp[1]))
print(mc_tbl, row.names = FALSE)
