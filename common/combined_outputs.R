# ============================================================
# combined_outputs.R
#
# Run order:  CHIKV_engine.R -> CHIKV_outputs.R
#             MAYV_engine.R  -> MAYV_outputs.R
#             (CHIKV_owsa.R, MAYV_owsa.R for figure 3)  ->  this file
# ============================================================
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(writexl)
                  library(readxl)})
dir.create("outputs", showWarnings = FALSE)

if (!exists("MAYV_EPI_SCENARIO")) MAYV_EPI_SCENARIO <- "high"   # "high" R0 2.1-2.9 | "low" R0 1.1-1.3

ARMS      <- c("Disease-blocking", "Disease + infection blocking")
OUT_LV    <- c("Cumulative symptomatic cases", "Cumulative hospitalisations",
               "Cumulative DALYs", "Cumulative deaths")
OUT_LAB   <- c("Symptomatic\ncases", "Hospitalisations", "DALYs", "Deaths")
FILL      <- c("No vaccination" = "grey60", "Vaccination" = "#4e79a7")
scen_cols <- c("No vaccination" = "grey55", "Disease-blocking" = "#4393c3",
               "Disease + infection blocking" = "#d6604d")

# ------------------------------------------------------------
# Font size
# ------------------------------------------------------------
FIG_W <- 190 / 25.4
FS <- list(
  panel_letter   = 10,
  disease_header = 9,
  block_ylab     = 8,

  # panel A
  epi_base       = 8,
  epi_axis_text  = 7,
  epi_axis_title = 8,
  epi_strip      = 8,   # grey strips
  epi_annot_chik = 2.5, # annotation
  epi_annot_mayv = 2.5,
  epi_legend     = 7,   # legend under panel A

  # panels B-E
  bur_base       = 8,
  bur_axis_x     = 7,
  bur_axis_y     = 7,
  bur_strip_x    = 7.5, # top grey strips
  bur_strip_y    = 7.5, # right grey strips

  # tornado figure
  owsa_base      = 9,
  owsa_title     = 9,
  owsa_strip     = 8
)

need1 <- function(f, how) if (!file.exists(f)) stop("missing ", f, " -- ", how) else f
G <- readRDS(need1("outputs/CHIKV_engine_results.rds", "run CHIKV_engine.R"))
M <- readRDS(need1(sprintf("outputs/MAYV_engine_results_%s.rds", MAYV_EPI_SCENARIO),
                   sprintf("run MAYV_lhs.R + MAYV_engine.R with R0_SCENARIO <- '%s'",
                           MAYV_EPI_SCENARIO)))
chik <- readRDS(need1("outputs/CHIKV_residual_burden.rds", "run CHIKV_outputs.R"))
mayv <- readRDS(need1("outputs/MAYV_residual_burden.rds",  "run MAYV_outputs.R"))
stopifnot(!is.null(M$wk_base))

check_mayv_scenario <- function(x, file) {
  s <- if (!is.null(attr(x, "R0_scenario"))) attr(x, "R0_scenario") else x[["R0_scenario"]]
  if (is.null(s))
    stop(file, " records no R0 scenario -- re-run the MAYV script that writes it.")
  if (!identical(s, MAYV_EPI_SCENARIO))
    stop(sprintf("%s is from the '%s' R0 scenario but MAYV_EPI_SCENARIO is '%s' -- re-run the MAYV chain with R0_SCENARIO <- '%s'.",
                 file, s, MAYV_EPI_SCENARIO, MAYV_EPI_SCENARIO))
}
check_mayv_scenario(M,    sprintf("outputs/MAYV_engine_results_%s.rds", MAYV_EPI_SCENARIO))
check_mayv_scenario(mayv, "outputs/MAYV_residual_burden.rds")

T_sim <- G$T_sim
mayv_lab <- sprintf("Mayaro (R0 %.1f-%.1f)", M$R0_lo, M$R0_hi)

# ------------------------------------------------------------
# 1. Epidemic curves
# ------------------------------------------------------------
qband <- function(mat, rows, lab, measure) {
  q <- apply(mat[rows, , drop = FALSE], 2, quantile, c(.025, .5, .975), na.rm = TRUE)
  data.frame(week = seq_len(ncol(mat)), lo = q[1,], med = q[2,], hi = q[3,],
             scenario = lab, measure = measure)
}
pair <- function(mat, rho, rows, lab)
  rbind(qband(mat, rows, lab, "True symptomatic"),
        qband(mat * rho, rows, lab, "Reported"))

ch_rows <- seq_len(nrow(G$wk_symp[["No vaccine (baseline)"]]))
mv_rows <- seq_len(M$N_DRAWS)
epi_ch <- rbind(
  pair(G$wk_symp[["No vaccine (baseline)"]],                    G$rho_i, ch_rows, "No vaccination"),
  pair(G$wk_symp[["pre-outbreak | Disease-blocking"]],          G$rho_i, ch_rows, ARMS[1]),
  pair(G$wk_symp[["pre-outbreak | Disease + infection blocking"]], G$rho_i, ch_rows, ARMS[2]))
epi_mv <- rbind(pair(M$wk_base, M$rho_draw, mv_rows, "No vaccination"),
                pair(M$wk_vacc, M$rho_draw, mv_rows, ARMS[1]))
for (nm in c("epi_ch", "epi_mv")) {
  d <- get(nm); d$measure  <- factor(d$measure, levels = c("True symptomatic", "Reported"))
  d$scenario <- factor(d$scenario, levels = names(scen_cols)); assign(nm, d)
}

pct_str <- function(b, v) { p <- 100*(b - v)/b; q <- quantile(p, c(.5, .025, .975), na.rm = TRUE)
  sprintf("%.1f %% (95%% UI %.1f-%.1f %%)", q[1], q[2], q[3]) }
cb <- G$per_draw[["No vaccine (baseline)"]][, "symptomatic"]
mb <- M$per_draw[["No vaccine (baseline)"]][mv_rows, "symptomatic"]
ann_ch <- paste0("% Reduction in predicted symptomatic cases\nDisease blocking only: ",
                 pct_str(cb, G$per_draw[["pre-outbreak | Disease-blocking"]][, "symptomatic"]),
                 "\nDisease & infection blocking: ",
                 pct_str(cb, G$per_draw[["pre-outbreak | Disease + infection blocking"]][, "symptomatic"]))
ann_mv <- paste0("% Reduction in predicted symptomatic cases\nDisease blocking only: ",
                 pct_str(mb, M$per_draw[[M$vac_name]][mv_rows, "symptomatic"]))

# vaccine rollout
ch_dose <- G$timings[["pre-outbreak"]] + 2 # median delay 2 weeks
mv_dose <- if (!is.null(M$dose_start))   M$dose_start            else ch_dose
mv_len  <- if (!is.null(M$deliv_median)) round(1/M$deliv_median) else 10

DPI_SUPP <- 1000
save_fig <- function(file, plot, width, height, dpi) {
  ggsave(file, plot, width = width, height = height, dpi = dpi)
  cat(sprintf("  saved %s (%.2f x %.2f in, %d dpi, %.1f MB)\n",
              file, width, height, dpi, file.size(file)/1e6))
  invisible(NULL)
}
save_pdf <- function(file, plot, width, height) {
  if (!isTRUE(capabilities("cairo")))
    stop("combined_master.pdf needs R's Cairo graphics (on macOS install XQuartz, then restart R).")
  ggsave(file, plot, width = width, height = height, device = cairo_pdf)
  cat(sprintf("  saved %s (%.2f x %.2f in, vector)\n", file, width, height))
  invisible(NULL)
}

epicurve <- function(d, strip, ann, xmin, xmax, ylab, ann_size = FS$epi_annot_chik) {
  d$panel <- factor(strip)
  TOP_EXP <- 0.20
  ann <- paste0("\n", ann)
  ggplot(d, aes(week, med, colour = scenario, fill = scenario)) +
    annotate("rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
             fill = "#3a7d3a", alpha = .12) +
    geom_area(data = subset(d, measure == "True symptomatic"),
              position = "identity", alpha = .25, colour = NA) +
    geom_line(data = subset(d, measure == "Reported"),
              aes(linetype = measure), linewidth = .8) +
    annotate("text", x = 1, y = Inf, label = ann, hjust = 0, vjust = 1,
             size = ann_size, lineheight = 1.15) +
    expand_limits(y = 0) +
    facet_wrap(~ panel) +
    scale_colour_manual(values = scen_cols, aesthetics = c("colour", "fill"), drop = FALSE) +
    scale_linetype_manual(values = c("Reported" = "dotted"),
                          labels = c("Reported" = "Reported symptomatic cases"), name = NULL) +
    guides(fill   = guide_legend(order = 1, nrow = 1,
                                 override.aes = list(linetype = 0, alpha = .55)),
           colour = guide_legend(order = 1, nrow = 1,
                                 override.aes = list(linetype = 0)),
           linetype = guide_legend(order = 2, nrow = 1,
                                   override.aes = list(colour = "grey25"))) +
    scale_x_continuous(breaks = c(1, seq(10, T_sim, by = 10))) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0.05, TOP_EXP))) +
    labs(x = "Week (index, 1 = 2025-W24)", y = ylab, colour = NULL, fill = NULL) +
    theme_bw(FS$epi_base) +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.box.just = "center",
          legend.spacing.y = unit(1, "pt"),
          strip.text = element_text(face = "bold", size = FS$epi_strip),
          legend.text = element_text(size = FS$epi_legend),
          legend.key.size = unit(.9, "lines"),
          axis.text  = element_text(size = FS$epi_axis_text),
          axis.title = element_text(size = FS$epi_axis_title),
          legend.justification = "center",
          panel.grid.minor = element_blank())
}

wrap_ann <- function(x, w) paste(vapply(strsplit(x, "\n")[[1]],
  function(l) paste(strwrap(l, width = w), collapse = "\n"), character(1)), collapse = "\n")

pE_c <- epicurve(epi_ch, "Chikungunya", wrap_ann(ann_ch, 60), ch_dose, ch_dose + 10,
                 "Predicted symptomatic cases")
pE_m <- epicurve(epi_mv, mayv_lab, wrap_ann(ann_mv, 60), mv_dose, mv_dose + mv_len, NULL)

p_epi <- pE_c + (pE_m + guides(fill = "none", colour = "none", linetype = "none")) +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box = "vertical", legend.box.just = "center",
        legend.spacing.y = unit(1, "pt"))
save_fig("outputs/CHIKV_MAYV_epicurves.png", p_epi, width = FIG_W, height = 3.4, dpi = 300)
cat(sprintf("Saved CHIKV_MAYV_epicurves.png (MAYV = '%s', R0 %.1f-%.1f, all %d draws).\n",
            MAYV_EPI_SCENARIO, M$R0_lo, M$R0_hi, M$N_DRAWS))

# ------------------------------------------------------------
# 2. Residual burden as a % of no vaccination
# ------------------------------------------------------------
chik$outcome <- factor(as.character(chik$outcome), levels = OUT_LV, labels = OUT_LAB)
mayv$outcome <- factor(as.character(mayv$outcome), levels = OUT_LV, labels = OUT_LAB)
mayv$arm     <- factor(ARMS[1], levels = ARMS[1])
chik$arm     <- factor(as.character(chik$arm), levels = ARMS)

burden <- function(d, strip, ylab, x_axis, y_axis, tag, row_strip = FALSE) {
  ggplot(d, aes(scenario, med, fill = scenario)) +
    geom_col(width = .6, na.rm = TRUE) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .35, na.rm = TRUE) +
    facet_grid(outcome ~ arm, drop = FALSE) +
    scale_fill_manual(values = FILL, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 105), breaks = seq(0, 100, 25)) +
    labs(tag = tag, x = NULL, y = ylab) +
    theme_bw(FS$bur_base) +
    theme(plot.tag = element_text(face = "bold", size = FS$panel_letter),
          plot.tag.position = "topleft", plot.tag.location = "margin",
          plot.margin = margin(t = 4, r = 5.5, b = 4, l = 2),
          axis.text.x  = if (x_axis) element_text(size = FS$bur_axis_x) else element_blank(),
          axis.ticks.x = if (x_axis) element_line() else element_blank(),
          axis.text.y  = if (y_axis) element_text(size = FS$bur_axis_y) else element_blank(),
          axis.ticks.y = if (y_axis) element_line() else element_blank(),
          strip.text.x = if (is.null(strip)) element_blank()
                         else element_text(face = "bold", size = FS$bur_strip_x),
          strip.text.y = if (row_strip) element_text(face = "bold", size = FS$bur_strip_y,
                                                    angle = -90)
                         else element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
}

one_row <- function(d, o) { d <- subset(d, outcome == o); d$outcome <- droplevels(d$outcome); d }

# The Mayaro model assumes CFR is 0, so blanked
blank_row <- function(o) {
  d <- data.frame(outcome = factor(o), arm = factor(ARMS[1], levels = ARMS[1]),
                  scenario = factor("No vaccination", levels = names(FILL)),
                  med = NA_real_, lo = NA_real_, hi = NA_real_)
  burden(d, NULL, NULL, FALSE, FALSE, NULL, row_strip = TRUE) +
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank())
}

row_burden <- function(o, tag, strip = FALSE, x_axis = FALSE) {
  st <- if (strip) TRUE else NULL
  c_blk <- burden(one_row(chik, o), st, NULL, x_axis, TRUE, tag)
  m_blk <- if (o == "Deaths") blank_row(o)
           else burden(one_row(mayv, o), st, NULL, x_axis, FALSE, NULL, row_strip = TRUE)
  list(c = c_blk, m = m_blk)
}

head_lab <- function(txt) ggplot() + labs(title = txt) + theme_void() +
  theme(plot.title = element_text(face = "bold", size = FS$disease_header, hjust = .5,
                                  margin = margin(t = 0, b = 2)),
        plot.margin = margin(0, 0, 0, 0))

YLAB <- "Cumulative burden (% of no vaccination)"
rB <- row_burden(OUT_LAB[1], "B", strip = TRUE)
rC <- row_burden(OUT_LAB[2], "C")
rD <- row_burden(OUT_LAB[3], "D")
rD$m <- burden(one_row(mayv, OUT_LAB[3]), NULL, NULL, TRUE, FALSE, NULL, row_strip = TRUE)
rE <- row_burden(OUT_LAB[4], "E", x_axis = TRUE)

# standalone figure
pair_row <- function(r) r$c + r$m + plot_layout(widths = c(2, 1))
burden_rows <- (head_lab("Chikungunya") + head_lab("Mayaro") +
                  plot_layout(widths = c(2, 1))) /
               pair_row(rB) / pair_row(rC) / pair_row(rD) / pair_row(rE) +
               plot_layout(heights = c(.01, 1, 1, 1, 1))

with_ylab <- function(blk) wrap_elements(patchworkGrob(blk)) +
  labs(tag = YLAB) +
  theme(plot.tag = element_text(size = FS$block_ylab, angle = 90), plot.tag.position = "left",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

p_burden <- with_ylab(burden_rows)
save_fig("outputs/combined_residual_burden.png", p_burden, width = FIG_W, height = 7.0, dpi = 300)
cat("Saved combined_residual_burden.png (B symptomatic / C hospitalisations / D DALYs / E deaths).\n")

# ------------------------------------------------------------
# 3. One-way sensitivity
# ------------------------------------------------------------
if (!all(file.exists("outputs/CHIKV_owsa.rds", "outputs/MAYV_owsa.rds"))) {
  cat("Skipped CHIKV_MAYV_owsa.png (run CHIKV_owsa.R and MAYV_owsa.R).\n")
} else {
  CO <- readRDS("outputs/CHIKV_owsa.rds"); MO <- readRDS("outputs/MAYV_owsa.rds")
  check_mayv_scenario(MO, "outputs/MAYV_owsa.rds")

  ca <- CO$owsa; ca$val <- ca$symptomatic
  bA <- setNames(CO$base$symptomatic, CO$base$arm)
  ca$base <- bA[ca$arm]; ca$arm <- factor(ca$arm, levels = ARMS)
  swA <- ca |> group_by(parameter) |> summarise(s = max(abs(val - base)), .groups = "drop")
  ca$parameter <- factor(ca$parameter, levels = swA$parameter[order(swA$s)])
  blA <- data.frame(arm = factor(names(bA), levels = ARMS), base = as.numeric(bA))

  mb <- MO$owsa; mb$val <- mb$averted; mb$base <- MO$base$averted
  mb$arm <- factor(ARMS[1], levels = ARMS[1])
  swB <- mb |> group_by(parameter) |> summarise(s = max(abs(val - base)), .groups = "drop")
  mb$parameter <- factor(mb$parameter, levels = swB$parameter[order(swB$s)])
  blB <- data.frame(arm = factor(ARMS[1], levels = ARMS[1]), base = MO$base$averted)

  tor <- function(d, bl, ttl, legend) {
    ggplot(d, aes(y = parameter)) +
      geom_vline(data = bl, aes(xintercept = base), linetype = "dashed", colour = "grey45") +
      geom_segment(aes(x = base, xend = val, yend = parameter, colour = bound),
                   linewidth = 3.2, alpha = .85) +
      facet_wrap(~ arm, nrow = 1, scales = "free_x") +
      scale_colour_manual(values = c(lower = "#d6604d", upper = "#4393c3"),
                          labels = c(lower = "Lower", upper = "Upper"), name = "Bound") +
      scale_x_continuous(labels = scales::comma) +
      labs(title = ttl, x = "Symptomatic cases averted", y = NULL) +
      theme_bw(FS$owsa_base) +
      theme(text = element_text(size = FS$owsa_base),
            legend.position = if (legend) "bottom" else "none",
            axis.text = element_text(size = 7.2), legend.text = element_text(size = 7.2),
            plot.title = element_text(face = "bold", size = FS$owsa_title),
            strip.text = element_text(face = "bold", size = FS$owsa_strip),
            panel.grid.minor = element_blank())
  }
  g <- tor(ca, blA, "A   Chikungunya", FALSE) + tor(mb, blB, "B   Mayaro", TRUE) +
    plot_layout(design = "AA\nB#", heights = c(1, 1.05))
  save_fig("outputs/CHIKV_MAYV_owsa.png", g, width = FIG_W, height = 6.2, dpi = DPI_SUPP)
  cat("Saved CHIKV_MAYV_owsa.png (A = CHIKV, B = MAYV).\n")
}

# ------------------------------------------------------------
# 4. Fig. 1
# ------------------------------------------------------------
pA_c <- pE_c + labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = FS$panel_letter),
        plot.tag.position = c(0, 1))

pA_m <- epicurve(epi_mv, "Mayaro", wrap_ann(ann_mv, 60), mv_dose, mv_dose + mv_len, NULL) +
  guides(fill = "none", colour = "none", linetype = "none")

row_A <- pA_c + pA_m + plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box = "vertical", legend.box.just = "center",
        legend.spacing.y = unit(1, "pt"),
        legend.text = element_text(size = FS$epi_legend), legend.key.size = unit(1.1, "lines"),
        plot.margin = margin(t = 5.5, r = 5.5, b = 2, l = 5.5))

row_BE <- plot_spacer() + with_ylab(burden_rows) + plot_spacer() +
  plot_layout(widths = c(.08, 1, .08))

master <- row_A / row_BE + plot_layout(heights = c(1.15, 3.6))

save_pdf("outputs/combined_master.pdf", master, width = FIG_W, height = 9.3)
cat("Saved combined_master.pdf (A epicurves, B symptomatic, C hospitalisations, D DALYs, E deaths).\n")


# ------------------------------------------------------------
# 5. Benefit per 100,000 doses
# ------------------------------------------------------------
per100k <- function(pd, base, o, doses) {
  r <- 1e5 * (base[, o] - pd[, o]) / doses; r[!is.finite(r) | r < 0] <- NA; r
}
OUT100 <- c(infections = "Infections", symptomatic = "Symptomatic",
            hospitalisations = "Hospitalisations", deaths = "Deaths", daly = "DALYs")
gb <- G$per_draw[["No vaccine (baseline)"]]; mbase <- M$per_draw[["No vaccine (baseline)"]]
ch_db <- G$per_draw[["pre-outbreak | Disease-blocking"]]
ch_di <- G$per_draw[["pre-outbreak | Disease + infection blocking"]]
mv_db <- M$per_draw[[M$vac_name]]
set.seed(20260802); perm <- sample.int(nrow(mv_db))

fmt100 <- function(v, d = 0) { q <- quantile(v, c(.5,.025,.975), na.rm = TRUE)
  sprintf("%s (%s-%s)", formatC(round(q[1],d), big.mark=",", format="f", digits=d),
          formatC(round(q[2],d), big.mark=",", format="f", digits=d),
          formatC(round(q[3],d), big.mark=",", format="f", digits=d)) }

rows <- lapply(names(OUT100), function(o) {
  a <- per100k(ch_db, gb, o, ch_db[, "doses"])
  b <- per100k(ch_di, gb, o, ch_di[, "doses"])
  m <- per100k(mv_db, mbase, o, mv_db[, "doses"])
  data.frame(outcome = OUT100[[o]],
             CHIKV_disease_blocking = fmt100(a),
             CHIKV_disease_and_infection = fmt100(b),
             MAYV_disease_blocking = fmt100(m),
             combined_CHIKV_db_plus_MAYV = fmt100(a + m[perm]),
             stringsAsFactors = FALSE) })

per100k_tbl <- do.call(rbind, rows)
p100_notes <- data.frame(item = c("Unit", "Why doses", "Why not a ratio", "Combined column",
                                  "Deaths", "Doses"),
  detail = c(
  "Outcomes averted per 100,000 doses administered, median (95% UI), computed per draw.",
  "Ixchiq is deployed once, so both models consume the same doses. Doses are the shared input, which makes per-dose benefit comparable across pathogens without either model being divided by the other.",
  "Dividing Mayaro cases averted by chikungunya cases averted compares two different models on a scale with no natural zero. Per-dose benefit has a natural zero (no benefit) and is additive.",
  "CHIKV disease-blocking + MAYV, added per draw with MAYV randomly permuted, since the two models' draws are not paired. Assumes the two epidemics are independent, and is conditional on BOTH occurring as modelled -- the Mayaro outbreak is hypothetical.",
  "MAYV deaths are identically zero (MAYV_ZERO_DEATHS = TRUE fixes the CFR at 0), so the Mayaro column is 0 by construction.",
  sprintf("CHIKV %s, MAYV %s (median), a %.2f%% difference from sampling alone.",
          format(round(median(ch_db[, "doses"])), big.mark=","),
          format(round(median(mv_db[, "doses"])), big.mark=","),
          100*abs(median(mv_db[,"doses"])/median(ch_db[,"doses"]) - 1))),
  stringsAsFactors = FALSE)
writexl::write_xlsx(list(averted_per_100k_doses = per100k_tbl, notes = p100_notes),
                    "outputs/combined_per_100k_doses.xlsx")
cat("Saved combined_per_100k_doses.xlsx (CHIKV beside MAYV, per 100,000 doses).\n")
print(per100k_tbl, row.names = FALSE)


# ------------------------------------------------------------
# 6. Proportion still symptomatic at three months (supplementary)
# ------------------------------------------------------------
source("common/common.R")
dpar <- load_daly_params()
set.seed(20260815); N_MC <- 2e5
rb <- function(p) rbeta(N_MC, p$a, p$b)

chronic_share <- function(p14, p90, p6, p12, p30) {
  a <- rb(p14); s <- rb(p90); ch <- rb(p6) + rb(p12) + rb(p30)
  tot <- a + s + ch
  list(acute = a/tot, subacute = s/tot, chronic = ch/tot, raw_total = tot)
}
cs_y <- chronic_share(dpar$p14_y, dpar$p90_y, dpar$p6_y, dpar$p12_y, dpar$p30_y)
cs_o <- chronic_share(dpar$p14_o, dpar$p90_o, dpar$p6_o, dpar$p12_o, dpar$p30_o)
mv_chronic <- rbeta(N_MC, 11, 5)

pc <- function(v, d = 1) { q <- quantile(v, c(.5, .025, .975))
  sprintf("%.*f%% (%.*f-%.*f)", d, 100*q[1], d, 100*q[2], d, 100*q[3]) }

chronic_tbl <- data.frame(
  pathogen  = c("Chikungunya", "Chikungunya", "Mayaro"),
  age_group = c("< 40 years", ">= 40 years", "All ages (not stratified)"),
  still_symptomatic_3m = c(pc(cs_y$chronic), pc(cs_o$chronic), pc(mv_chronic)),
  distribution = c("Sum of Beta(799.3, 3307.0), Beta(77.6, 836.7) and Beta(7.6, 1075.0), renormalised",
                   "Sum of Beta(5998.9, 21654.9), Beta(184.7, 1216.1) and Beta(15.8, 516.3), renormalised",
                   "Beta(11, 5), the chronic margin of Dirichlet(1, 4, 11)"),
  source = c("O'Driscoll et al. 2021 IJID (pooled cohorts)",
             "O'Driscoll et al. 2021 IJID (pooled cohorts)",
             "Halsey et al. 2013 (Peru, n = 16)"),
  stringsAsFactors = FALSE)

partition_tbl <- data.frame(
  pathogen  = c(rep("Chikungunya", 2), "Mayaro"),
  age_group = c("< 40 years", ">= 40 years", "All ages (not stratified)"),
  resolved_by_14d      = c(pc(cs_y$acute),    pc(cs_o$acute),    pc(rbeta(N_MC, 1, 15))),
  resolved_14d_to_3m   = c(pc(cs_y$subacute), pc(cs_o$subacute), pc(rbeta(N_MC, 4, 12))),
  still_symptomatic_3m = c(pc(cs_y$chronic),  pc(cs_o$chronic),  pc(mv_chronic)),
  stringsAsFactors = FALSE)

chronic_notes <- data.frame(item = c(
  "Quantity", "Chikungunya derivation", "Mayaro derivation", "Marginal, not conditional",
  "Why medians are not summed", "Interval widths are not comparable",
  "Mayaro is prevalence, not survival", "Mayaro acute bin", "No pooled value"),
  detail = c(
  "Proportion of symptomatic cases still reporting arthralgia at approximately three months from symptom onset. Median (95% UI) over 200,000 Monte Carlo draws.",
  "One minus the cumulative probability of recovery by 90 days, i.e. the sum of the 6-month, 12-month and 30-month recovery rows, renormalised so the three phases sum to 1. The five O'Driscoll rows sum to 1.000 (<40) and 1.000 (>=40) using the published Beta hyperparameters, confirming they partition one cohort.",
  "Proportion with arthralgia at the 3-month visit, 11 of 16 patients. Drawn as the chronic component of a Dirichlet(1, 4, 11) so the three phases stay correlated and sum to 1 in every draw; the exact marginal is Beta(11, 5).",
  "Despite the 'within X after Y period' wording of the source tables, these are marginal shares of recovery time, not conditional recovery probabilities. Reading row 2 as conditional on not having recovered in row 1 gives a different number (26.7% rather than 23.9% for Mayaro).",
  "Components are summed per draw, not by adding medians. Medians are not additive in general; adding the three components' UI endpoints would assume perfect positive correlation and yields an interval about 1.6 times too wide, since the chikungunya Betas are drawn independently.",
  "The chikungunya intervals are narrow because O'Driscoll's Beta hyperparameters encode very large pooled sample sizes; they describe precision in the source meta-analysis, not transportability to Caldas Novas. The Mayaro interval reflects n = 16. The widths differ by roughly 17-fold and should not be read as relative confidence.",
  "Halsey reports arthralgia prevalence falling to 3/16 at day 20 and rising to 11/16 at 3 months, a biphasic pattern the authors treat as real. Recovery is therefore not absorbing and the 3-month figure is a POINT PREVALENCE, not the probability of never having recovered.",
  "Halsey's first follow-up visit was day 20, so the 14-day boundary was never observed. The acute share is 'no arthralgia at the acute visit' (1 of 16), which is why it is far below the chikungunya values.",
  "No age-pooled chikungunya value is given: the source stratifies at 40 years and a pooled figure would depend on the age distribution of symptomatic cases, making it a model output rather than a parameter."),
  stringsAsFactors = FALSE)

write_xlsx(list(still_symptomatic_3m = chronic_tbl,
                full_partition       = partition_tbl,
                notes                = chronic_notes),
           "outputs/combined_chronic_share.xlsx")
cat("\nSaved combined_chronic_share.xlsx (proportion still symptomatic at 3 months).\n")
print(chronic_tbl[, 1:3], row.names = FALSE)
cat("\nFull partition:\n"); print(partition_tbl, row.names = FALSE)
cat(sprintf("\nO'Driscoll rows sum to %.4f (<40) and %.4f (>=40) before renormalisation.\n",
            median(cs_y$raw_total), median(cs_o$raw_total)))
