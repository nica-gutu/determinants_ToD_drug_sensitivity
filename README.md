# Determinants of time-of-day drug sensitivity

Analysis and modeling code supporting the study **“A combined mathematical and experimental approach reveals the drivers of time-of-day drug sensitivity in human cells”** (*Communications Biology*, 2025).

**Publication:** https://doi.org/10.1038/s42003-025-07931-1

## Scientific question

Drug responses can vary with the time of day, but the observed temporal profile reflects several interacting biological and experimental factors. This project uses a combined mathematical and experimental framework to determine how **circadian properties, pharmacodynamic response characteristics, cellular growth, drug stability, and assay timing** shape time-of-day drug sensitivity in human cells.

The goal is to separate these contributions quantitatively and identify which properties determine the magnitude and timing of observed treatment-response rhythms.

## Repository contents

The repository is organized by the major determinants examined in the study.

### `Circadian_properties/`

Analysis and modeling of circadian parameters and their effect on time-of-day response profiles, including damped-oscillator models, fitting of luminescence time series, and exploration of circadian amplitude, period, phase, and damping.

### `Pharmacodynamics1/`

Analyses of how pharmacodynamic properties influence time-of-day drug sensitivity, including cytostatic and cytotoxic response models and relationships between drug-response parameters and temporal treatment effects.

### `Pharmacodynamics2/`

Analyses of additional pharmacodynamic and experimental factors, including drug stability, drug age, and survival responses.

### `Cell_growth_dynamics/`

Analyses of cell growth and density effects on measured treatment responses.

### `Evaluation_time/`

Analyses of how assay evaluation time influences inferred drug-response profiles.

### `Optimal_treatment_time/`

Model-based analyses of the treatment times associated with maximal predicted benefit under different circadian and pharmacodynamic conditions.

The scripts correspond to specific analyses from the study and are intended as research-analysis code rather than a general-purpose software package.

## Computational approaches

Methods represented in the repository include:

- mathematical modeling of circadian modulation and pharmacodynamic response;
- nonlinear curve fitting;
- damped-oscillator and dose-response models;
- correlation and statistical analyses;
- parameter sweeps and sensitivity analyses;
- simulation of time-of-day treatment-response profiles;
- quantitative analysis of cell-growth and assay-timing effects.

## Requirements

The scripts are written in Python. Core packages used across the repository include:

```text
numpy
pandas
matplotlib
seaborn
scipy
pyboat
```

Some analyses read Excel workbooks and may require a compatible pandas Excel engine.

## Data and reproducibility

The raw experimental datasets are **not included in the repository**. Scripts that use experimental data expect the corresponding input files (for example, files referenced under `Raw_Data/`) to be available locally.

Output paths and analysis parameters are defined within the individual scripts. To reproduce a specific analysis, provide the corresponding experimental input data, verify the file/path configuration, and run the script associated with the relevant project section.

## Citation

If you use this code or build on these analyses, please cite:

> Gutu, N., Ishikuma, H., Ector, C. et al. **A combined mathematical and experimental approach reveals the drivers of time-of-day drug sensitivity in human cells.** *Communications Biology* 8, 491 (2025). https://doi.org/10.1038/s42003-025-07931-1

## Contact

**Nica Gutu**  
Computational Biology / Data Science  
https://nica-gutu.github.io/website/
