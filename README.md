# Incremental health impact of chikungunya vaccine cross-protection against Mayaro virus disease in Brazil: a modeling study

The code fits an age-structured age-stratified SEIRV (susceptible-exposed-infectious-recovered-vaccinated) 
model to the 2025–26 chikungunya outbreak in municipality A, estimates the health impact of vaccination with IxchiqTM, and projects
the impact of the same vaccine against a hypothetical Mayaro virus outbreak in the same
municipality. Outcomes are infections, symptomatic cases, hospitalizations, deaths and
DALYs, over 52 epidemiological weeks (2025-W24 to 2026-W22, 8 June 2025 to 6 June 2026).

## Repository layout

```
CHIKV/    chikungunya model: fit, Monte Carlo engine, outputs, sensitivity analyses
MAYV/     Mayaro model: rainfall and seasonal envelope, Monte Carlo engine, outputs,
          sensitivity analyses
common/   code shared by both models (common.R) and the figures that combine them
          (combined_outputs.R)
data/     input data -- not in this repository; download it (see Data) and put it here
outputs/  created when the code runs; every table and figure is written here
          (and outputs/lhs_draws/, the sampled parameter values behind every estimate)
run_all.R runs the whole analysis in the right order
```

## Requirements

- R 4.1 or later (developed and tested with R 4.6.0).
- CRAN packages: `readxl`, `writexl`, `dplyr`, `tidyr`, `ggplot2`, `scales`, `patchwork`.
  Install them with
  ```r
  install.packages(c("readxl", "writexl", "dplyr", "tidyr", "ggplot2", "scales", "patchwork"))
  ```
- **macOS only:** [XQuartz](https://www.xquartz.org), which R needs to write the main-text
  figure as a PDF with embedded fonts. Windows and Linux builds of R need nothing extra.

## Data

The input files are deposited on Mendeley Data:
[doi:10.17632/zk98x2736j.1](https://doi.org/10.17632/zk98x2736j.1). Download them into
`data/`:

| File | Contents | Source |
|---|---|---|
| `weekly_case.xlsx` | Weekly reported chikungunya cases by age group, munucipality A (sheet `ca_combined`) | SINAN, Brazilian Ministry of Health |
| `population.xlsx` | Population by age group, munucipality A (sheet `prop_immune`) | IBGE 2022 census |
| `disease_progression.xlsx` | Severity, recovery, disability-weight and life-expectancy parameters | Kang et al. 2025 (Table S4) and sources cited in the supplementary materials |
| `model_calibration.xlsx` | Natural-history and reporting priors | Sources cited in the supplementary materials |
| `rainfall_caldas.csv` | Climatological rainfall for Caldas Novas | CHIRPS, via the Humanitarian Data Exchange (WFP, Brazil subnational rainfall) |

The Mendeley deposit also contains the sampled parameter values behind every published
estimate (the Latin hypercube draws). `run_all.R` writes the same files to `outputs/lhs_draws/`:

| File | Contents |
|---|---|
| `CHIKV_lhs_draws.csv` | Chikungunya transmission draws (one row per Latin hypercube draw) |
| `CHIKV_burden_draws.csv` | Chikungunya vaccine, severity and DALY draws, linked to the transmission draws by `transmission_row` |
| `MAYV_lhs_draws_low.csv`, `MAYV_lhs_draws_high.csv` | Mayaro transmission draws, one file per R0 scenario |
| `MAYV_burden_draws_low.csv`, `MAYV_burden_draws_high.csv` | Mayaro vaccine, severity and DALY draws, one file per R0 scenario |

## Running the analysis

Open the project in RStudio (`summer_project.Rproj`), or set R's working directory to the
repository's top folder, then run

```r
source("run_all.R")
```

This checks that the input files are in `data/`, then runs every script in order, each in
a fresh R session. All results are written to `outputs/`, and each step's console output to
`outputs/logs/`. It takes about 1.5 hours on a laptop; most of that is the chikungunya
model fit and the two analyses that re-fit it 1,000 times.

Follow the following order (each step reads what the previous one wrote):

| Step | Script | What it does |
|---|---|---|
| 1 | `CHIKV/CHIKV_lhs.R` | Fits the chikungunya model to the case data; 1,000 Latin hypercube draws |
| 2 | `CHIKV/CHIKV_engine.R` | Vaccination scenarios and health outcomes, draw by draw |
| 3 | `CHIKV/CHIKV_outputs.R` | Chikungunya tables and figures |
| 4 | `CHIKV/CHIKV_owsa.R` | One-way sensitivity analysis; campaign timing × coverage (Fig. S6) |
| 5 | `CHIKV/CHIKV_rho_sensitivity.R` | Reporting-rate sensitivity (Discussion) |
| 6 | `MAYV/rainfall_season.R` | Climatological rainfall envelope |
| 7 | `MAYV/MAYV_build_hybrid_envelope.R` | Seasonal transmission envelope for Mayaro (Fig. S4) |
| 8 | `MAYV/MAYV_lhs.R` | Mayaro outbreak draws, run once per R0 scenario (`R0_SCENARIO <- "low"`, then `"high"`) |
| 9 | `MAYV/MAYV_engine.R` | Mayaro vaccination scenario and health outcomes, after each step 8 |
| 10 | `MAYV/MAYV_outputs.R`, `MAYV/MAYV_burden_outputs.R` | Mayaro tables and figures |
| 11 | `MAYV/MAYV_owsa.R` | One-way sensitivity analysis |
| 12 | `MAYV/MAYV_threshold_sensitivity.R` | Transmission-potential and seeding-week sensitivity |
| 13 | `common/combined_outputs.R` | Figures and tables combining both models |

`MAYV/extract_caldas_rainfall.R` is a one-off that extracted `rainfall_caldas.csv` from the
full CHIRPS download; it is not part of the pipeline.

## Where each result comes from

Sheet names are in brackets. Unless stated otherwise, Mayaro results are from the high
transmission scenario (the second `MAYV/MAYV_lhs.R` run).

**Main text**

| Item | File in `outputs/` [sheet] | Script |
|---|---|---|
| Fig. 1 | Conceptual diagram, not produced by the code | – |
| Fig. 2 | `combined_master.pdf`; the percentages in panels B–E are in `CHIKV_residual_burden.xlsx` and `MAYV_residual_burden.xlsx` | `common/combined_outputs.R` (figure); `CHIKV/CHIKV_outputs.R`, `MAYV/MAYV_outputs.R` (workbooks) |
| Table 1 | `CHIKV_vacc_outputs.xlsx` and `MAYV_vacc_outputs.xlsx` [averted_per_100k_doses]; `CHIKV_nnv_outputs.xlsx` and `MAYV_nnv_outputs.xlsx` [nnv] | `CHIKV/CHIKV_outputs.R`, `MAYV/MAYV_outputs.R` |
| Results, % of burden averted | `CHIKV_residual_burden.xlsx`, `MAYV_residual_burden.xlsx` | `CHIKV/CHIKV_outputs.R`, `MAYV/MAYV_outputs.R` |
| Results, Mayaro mean Rn | `MAYV_threshold_sensitivity.xlsx` [transmission_potential] | `MAYV/MAYV_threshold_sensitivity.R` |
| Discussion, attack rates under the national and municipal reporting rates | `CHIKV_rho_sensitivity.xlsx` [summary] | `CHIKV/CHIKV_rho_sensitivity.R` |
| Discussion, chikungunya R0 | `lhs_draws/CHIKV_lhs_draws.csv` [`R0_peak` column] | `CHIKV/CHIKV_lhs.R` |

**Supplementary materials**

| Item | File in `outputs/` [sheet] | Script |
|---|---|---|
| Fig. S1 | Model schematic, not produced by the code | – |
| Table S1 | Inputs: `data/model_calibration.xlsx` and `data/disease_progression.xlsx` (chikungunya symptomatic fraction); the other priors (FOI, reporting rate, R0, Mayaro symptomatic fraction, vaccine parameters) are set in the scripts | `CHIKV/CHIKV_lhs.R`, `CHIKV/CHIKV_engine.R`, `MAYV/MAYV_lhs.R`, `MAYV/MAYV_engine.R` |
| Fig. S2 | Map, not produced by the code | – |
| Fig. S3 | `CHIKV_fit_cases_infections.png` | `CHIKV/CHIKV_outputs.R` |
| Fig. S4 | `MAYV_hybrid_envelope.png` | `MAYV/MAYV_build_hybrid_envelope.R` |
| Fig. S5 | Prognosis diagram, not produced by the code | – |
| Table S2 | Inputs: `data/disease_progression.xlsx`; Mayaro recovery shares (Halsey et al. 2013) are set in the script and summarised in `combined_chronic_share.xlsx` | `MAYV/MAYV_engine.R`, `common/combined_outputs.R` |
| Table S3 | `CHIKV_vacc_outputs.xlsx` and `MAYV_vacc_outputs.xlsx` [baseline_true_reported]; `CHIKV_daly_outputs.xlsx` and `MAYV_daly_outputs.xlsx` [daly_by_scenario] | `CHIKV/CHIKV_outputs.R`, `MAYV/MAYV_outputs.R` |
| Table S4 | `CHIKV_vacc_outputs.xlsx` [doses_wastage] | `CHIKV/CHIKV_outputs.R` |
| Fig. S6 | `CHIKV_timing_coverage.png`; values in `CHIKV_owsa.xlsx` [intended_vs_actual] | `CHIKV/CHIKV_owsa.R` |
| Table S5 | `MAYV_threshold_sensitivity.xlsx` [transmission_potential] (both scenarios) | `MAYV/MAYV_threshold_sensitivity.R` |
| Table S6 | `MAYV_owsa.xlsx` [r0_response_table] | `MAYV/MAYV_owsa.R` |
| Table S7 | `CHIKV_daly_outputs.xlsx` and `MAYV_daly_outputs.xlsx` [daly_by_scenario, yld_by_phase, daly_averted] | `CHIKV/CHIKV_outputs.R`, `MAYV/MAYV_outputs.R` |
| Fig. S7 | `CHIKV_MAYV_owsa.png`; values in `CHIKV_owsa.xlsx` and `MAYV_owsa.xlsx` [owsa] | `CHIKV/CHIKV_owsa.R`, `MAYV/MAYV_owsa.R` (analyses); `common/combined_outputs.R` (figure) |
| Fig. S8 | `MAYV_ve_coverage_heatmap.png`; values in `MAYV_ve_coverage_sweep.xlsx` | `MAYV/MAYV_owsa.R` |
| Table S8 | `MAYV_threshold_sensitivity.xlsx` [seeding_week] | `MAYV/MAYV_threshold_sensitivity.R` |

## Reproducibility

Every random step uses a fixed seed, so a re-run reproduces the published numbers exactly.

## Licence

The code is released under the MIT licence (see `LICENSE`). The data on Mendeley carry
their own licence.

## Citation

If you use this code, please cite the manuscript.

## Contact

Chloe Lee [chloe.lee@alumni.lshtm.ac.uk]
