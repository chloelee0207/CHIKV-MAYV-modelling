# ============================================================
# Build the MAYV transmission seasonality envelope from CHIRPS rainfall
# ============================================================
suppressMessages({library(dplyr); library(ggplot2)})
dir.create("outputs", showWarnings = FALSE)  

LAG_WEEKS <- 0   

# ------------------------------------------------------------
# 1. 52 epi weeks 2025-W24 -> 2026-W22
# ------------------------------------------------------------
win <- rbind(data.frame(Year = 2025, week = 24:53),
             data.frame(Year = 2026, week = 1:22))
win$week_index <- seq_len(nrow(win))
win$week_label <- sprintf("%d-W%02d", win$Year, win$week)
stopifnot(nrow(win) == 52)
w01 <- ifelse(win$Year == 2025, as.Date("2024-12-29"), as.Date("2026-01-04"))
win$date <- as.Date(w01, origin = "1970-01-01") + (win$week - 1) * 7 + 3

# ------------------------------------------------------------
# 2. Climatological rainfall grid
# ------------------------------------------------------------
r <- read.csv("data/rainfall_caldas.csv"); r$date <- as.Date(r$date)
r$md <- format(r$date, "%m-%d")
clim <- r %>% group_by(md) %>% summarise(rfh_avg = first(rfh_avg), .groups = "drop")
grid <- do.call(rbind, lapply(c(2025, 2026), function(y)
  data.frame(date = as.Date(paste0(y, "-", clim$md)) + 4, rain = clim$rfh_avg)))
grid <- grid[order(grid$date), ]

# ------------------------------------------------------------
# 3. Interpolate to weekly, apply lag, normalise to mean 1
# ------------------------------------------------------------
win$rain   <- approx(grid$date, grid$rain, xout = win$date - LAG_WEEKS * 7, rule = 2)$y
season     <- win$rain / mean(win$rain)          # mean-1 envelope
win$season <- season
stopifnot(abs(mean(season) - 1) < 1e-9, length(season) == 52)

cat(sprintf("Rainfall envelope: lag %d wk | peak %s | range [%.2f, %.2f] | mean %.3f\n",
            LAG_WEEKS, win$week_label[which.max(season)], min(season), max(season), mean(season)))

# ------------------------------------------------------------
# 4. Save
# ------------------------------------------------------------
saveRDS(season, "outputs/caldas_rain_season.rds")
write.csv(win[, c("week_index", "week_label", "rain", "season")],
          "outputs/caldas_rain_season.csv", row.names = FALSE)

year_break <- 30.5   # between 2025-W53 (idx 30) and 2026-W01 (idx 31)
tick_idx   <- c(7, 17, 27, 40, 50)
x_ticks    <- data.frame(week_index = tick_idx, week_label = win$week_label[tick_idx])
p <- ggplot(win, aes(week_index, season)) +
  annotate("rect", xmin = min(which(season >= 1)) - 0.5, xmax = max(which(season >= 1)) + 0.5,
           ymin = -Inf, ymax = Inf, fill = "#cfe6f2", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  geom_line(colour = "#2c7fb8", linewidth = 1) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week_label) +
  labs(x = "Epi week", y = "Transmission envelope (mean = 1)",
       title = "Caldas Novas MAYV seasonal envelope from CHIRPS climatological rainfall",
       subtitle = sprintf("rfh_avg, dekad->week, lag %d wk, normalised; shaded = above-mean (wet)", LAG_WEEKS)) +
  theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5),
                       plot.subtitle = element_text(hjust = 0.5, size = 9),
                       panel.grid.minor = element_blank())
ggsave("outputs/caldas_rain_season.png", p, width = 8, height = 4.2, dpi = 120)
cat("Saved caldas_rain_season.rds, caldas_rain_season.csv, caldas_rain_season.png\n")
print(win[, c("week_index", "week_label", "rain", "season")], row.names = FALSE)
