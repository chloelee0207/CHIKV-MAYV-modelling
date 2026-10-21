# Input data

The input files are not stored in this repository. Download them from Mendeley Data
([doi:10.17632/zk98x2736j.1](https://doi.org/10.17632/zk98x2736j.1)) and place them in
this folder:

| File | Contents |
|---|---|
| `weekly_case.xlsx` | Weekly reported chikungunya cases by age group, Caldas Novas (sheet `ca_combined`) |
| `population.xlsx` | Population by age group, Caldas Novas (sheet `prop_immune`) |
| `disease_progression.xlsx` | Severity, recovery, disability-weight and life-expectancy parameters |
| `model_calibration.xlsx` | Natural-history and reporting priors |
| `rainfall_caldas.csv` | Climatological rainfall for Caldas Novas (CHIRPS, via the Humanitarian Data Exchange) |

`run_all.R` checks that all five are here before it starts.
