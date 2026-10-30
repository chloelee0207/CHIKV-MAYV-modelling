# ============================================================
# MAYV_build_hybrid_envelope.R
# ------------------------------------------------------------
# Output: caldas_hybrid_season.rds & MAYV_hybrid_envelope.png (Fig. S4)
# ============================================================
suppressMessages(library(ggplot2))
dir.create("outputs", showWarnings = FALSE)

beta <- readRDS("outputs/caldas_beta_season.rds")
rain <- readRDS("outputs/caldas_rain_season.rds")
b_at <- function(o) beta[o - 23]
r_at <- function(o) rain[o - 23]

TRANS <- 63 # join at 2026-W10 because we want to constrain beta; beyond W10 beta goes up because susceptibles are depleting towards the end of the outbreak 
ords  <- 24:75
scale <- b_at(TRANS) / r_at(TRANS)
hybrid <- sapply(ords, function(o) if (o <= TRANS) b_at(o) else r_at(o) * scale)
hybrid <- hybrid / mean(hybrid)             # store mean-1

stopifnot(length(hybrid) == 52, abs(mean(hybrid) - 1) < 1e-9)
saveRDS(hybrid, "outputs/caldas_hybrid_season.rds")
cat(sprintf("Wrote caldas_hybrid_season.rds: %d wks, mean=%.3f, max=%.2f, min-after-peak=%.3f\n",
            length(hybrid), mean(hybrid), max(hybrid), min(hybrid[which.max(hybrid):length(hybrid)])))

wk_lab <- function(o) ifelse(o <= 53, o, o - 53)
LAB_BETA <- "CHIKV fitted weekly transmission rates \u03b2(t)"
LAB_HYB  <- "Hybrid seasonal envelope"
LAB_RAIN <- "Average rainfall in Caldas Novas"
ENV_COL  <- c("#F8766D", "#00BA38", "#619CFF")
names(ENV_COL) <- c(LAB_BETA, LAB_HYB, LAB_RAIN)

df <- rbind(
  data.frame(ord = 24:75, val = beta,   env = LAB_BETA),
  data.frame(ord = 24:75, val = rain,   env = LAB_RAIN),
  data.frame(ord = ords,  val = hybrid, env = LAB_HYB))
df$env <- factor(df$env, levels = names(ENV_COL))
p <- ggplot(df, aes(ord, val, colour = env)) + geom_line(linewidth = .7) +
  geom_vline(xintercept = TRANS, linetype = "dashed", colour = "grey50") +
  annotate("text", x = TRANS, y = Inf, label = "join\n2026-W10", vjust = 1.2, size = 2.5) +
  scale_x_continuous(breaks = c(24, 40, 50, 54, 63, 75),
                     labels = wk_lab(c(24, 40, 50, 54, 63, 75))) +
  scale_colour_manual(values = ENV_COL, name = NULL) +
  labs(x = "Week", y = "transmission envelope (mean-1)", colour = NULL) +
  theme_bw(9) + theme(legend.position = "bottom", legend.text = element_text(size = 7.2))
# supplementary figure (Fig. S4)
ggsave("outputs/MAYV_hybrid_envelope.png", p, width = 190 / 25.4, height = 3.6, dpi = 1000)
cat("Wrote MAYV_hybrid_envelope.png\n")
