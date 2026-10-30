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

| Manuscript item | File in `outputs/` | Script |
|---|---|---|
| Figure 1 | `combined_master.pdf` | `common/combined_outputs.R` |
| Table 2 | `CHIKV_vacc_outputs.xlsx`, `MAYV_vacc_outputs.xlsx`, `*_daly_outputs.xlsx` | `CHIKV/CHIKV_outputs.R`, `MAYV/MAYV_outputs.R` |
| Fig. S4 | `MAYV_hybrid_envelope.png` | `MAYV/MAYV_build_hybrid_envelope.R` |
| Fig. S6 | `CHIKV_timing_coverage.png` | `CHIKV/CHIKV_owsa.R` |
| [other supplementary items] | | |

## Reproducibility

Every random step uses a fixed seed, so a re-run reproduces the published numbers exactly.

## Licence

The code is released under the MIT licence (see `LICENSE`). The data on Mendeley carry
their own licence.

## Citation

If you use this code, please cite the manuscript.

## Contact

Chloe Lee [chloe.lee@alumni.lshtm.ac.uk]
