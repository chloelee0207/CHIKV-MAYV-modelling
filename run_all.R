# ============================================================
# run_all.R -- reproduce every result in the manuscript, in order.
#
# Run from the repository's top folder (open summer_project.Rproj, or setwd() there):
#     source("run_all.R")
# Inputs are read from data/ (see data/README.md); every table and figure is written to
# outputs/, and each step's console output to outputs/logs/. Each step runs in its own fresh R session, so no setting from one script can
# leak into the next. About 1.5 hours on a laptop, mostly the chikungunya fit (step 1)
# and the two analyses that re-fit it 1,000 times (steps 4 and 5).
# ============================================================
inputs <- c("weekly_case.xlsx", "population.xlsx", "disease_progression.xlsx",
            "model_calibration.xlsx", "rainfall_caldas.csv")
missing <- inputs[!file.exists(file.path("data", inputs))]
if (length(missing))
  stop("Missing input file(s) in data/: ", paste(missing, collapse = ", "),
       ". See data/README.md for where to download them.")

# Each step: a script, plus any setting it needs (R0_SCENARIO picks the Mayaro scenario).
steps <- list(
  list("CHIKV/CHIKV_lhs.R"),                    # fit to the 2025-26 outbreak, 1,000 draws
  list("CHIKV/CHIKV_engine.R"),                 # vaccination scenarios, draw by draw
  list("CHIKV/CHIKV_outputs.R"),                # chikungunya tables and figures
  list("CHIKV/CHIKV_owsa.R"),                   # one-way sensitivity; timing x coverage (Fig S6)
  list("CHIKV/CHIKV_rho_sensitivity.R"),        # reporting-rate sensitivity (Discussion)
  list("MAYV/rainfall_season.R"),               # climatological rainfall envelope
  list("MAYV/MAYV_build_hybrid_envelope.R"),    # Mayaro seasonal envelope (Fig S4)
  list("MAYV/MAYV_lhs.R",    'R0_SCENARIO <- "low"'),
  list("MAYV/MAYV_engine.R"),
  list("MAYV/MAYV_lhs.R",    'R0_SCENARIO <- "high"'),   # high last: the untagged results
  list("MAYV/MAYV_engine.R"),                             # files are the main scenario
  list("MAYV/MAYV_outputs.R"),                  # Mayaro tables and figures
  list("MAYV/MAYV_burden_outputs.R"),           # transmission-only Mayaro workbook
  list("MAYV/MAYV_owsa.R"),                     # one-way sensitivity and R0 response
  list("MAYV/MAYV_threshold_sensitivity.R"),    # transmission potential, seeding week
  list("common/combined_outputs.R")             # figures and tables combining both models
)

rscript <- file.path(R.home("bin"), "Rscript")
dir.create(file.path("outputs", "logs"), recursive = TRUE, showWarnings = FALSE)
t_all <- Sys.time()
for (i in seq_along(steps)) {
  s <- steps[[i]]
  expr <- paste0(if (length(s) > 1) paste0(s[[2]], "; ") else "", "source('", s[[1]], "')")
  cat(sprintf("[%2d/%d] %s%s ... ", i, length(steps), s[[1]],
              if (length(s) > 1) paste0(" (", s[[2]], ")") else ""))
  t0 <- Sys.time()
  log <- file.path("outputs", "logs", sprintf("%02d_%s.log", i, sub("[.]R$", "", basename(s[[1]]))))
  status <- system2(rscript, c("-e", shQuote(expr)), stdout = log, stderr = log)
  if (!identical(status, 0L))
    stop("Step ", i, " (", s[[1]], ") failed; see ", log, " for the error.")
  cat(sprintf("done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat(sprintf("\nAll %d steps finished in %.0f minutes. Results are in outputs/.\n",
            length(steps), as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
