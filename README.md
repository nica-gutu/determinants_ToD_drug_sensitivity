# Determinants of time-of-day drug sensitivity

Analysis and modeling code supporting:

**Gutu, N. et al. “A combined mathematical and experimental approach reveals the drivers of time-of-day drug sensitivity in human cells.” _Communications Biology_ 8, 491 (2025).**  
https://doi.org/10.1038/s42003-025-07931-1

## Overview

This project investigates why the measured effect of a drug can vary depending on **time of treatment** and which biological or experimental properties determine the magnitude and shape of this time-of-day (ToD) response.

Rather than attributing temporal drug sensitivity to the circadian clock alone, the study develops a general mathematical and experimental framework that separates the contributions of:

- circadian amplitude, period, phase, and damping;
- drug concentration;
- maximal drug effect and dose-response steepness;
- drug stability;
- cell-growth dynamics;
- assay duration and evaluation time.

The framework combines mathematical modeling, simulations, long-term live-cell imaging, and analysis of drug-response data across multiple human cell lines.

## Modeling framework

The computational model describes cell-population growth under pharmacological treatment while allowing the effective drug concentration to be modulated by a circadian signal.

Both **cytostatic** and **cytotoxic** drug effects are represented using dose-response relationships. Parameter sweeps are then used to determine how circadian and pharmacodynamic properties reshape the amplitude and timing of ToD sensitivity.

Circadian parameters are also extracted from experimental bioluminescence recordings by fitting damped oscillatory functions, while drug-response curves are parameterized using quantities such as **IC50, Emax, and the Hill coefficient**.

## Experimental context

The study integrates previously generated and new data from multiple human cell models, including breast cancer cell lines and U2OS cells.

Experimental analyses include drugs such as **Alpelisib, Paclitaxel, Alisertib, Adavosertib, Torin2, Cisplatin, Doxorubicin, and 5-FU**, depending on the analysis.

For long-term assay analyses, cells treated with Doxorubicin, Cisplatin, and 5-FU were followed by live imaging for up to 5 days, enabling direct assessment of how the apparent dose-response relationship changes with evaluation time.

## Main findings represented by the analyses

- Increasing circadian amplitude increases the magnitude of time-of-day drug-response differences, whereas circadian period and damping affect the temporal profile differently.
- The amplitude of the ToD response depends strongly on the location of the reference dose along the dose-response curve and is maximized near an intermediate, approximately half-maximal drug-effect region rather than simply increasing with dose.
- Drugs with stronger maximal effects and steeper dose-response relationships can produce larger time-of-day differences.
- Experimental drug-response trends across several drug/cell-line combinations are consistent with these model predictions.
- Cell-growth dynamics affect how drug response should be normalized and interpreted over time.
- **Assay duration is itself a major determinant of the measured drug response**: IC50, Emax, and Hill estimates change with evaluation time and approach more stable values at later measurements.
- Longer evaluation times can reveal larger ToD-response differences that may be underestimated in shorter assays.

The study therefore shows that observed time-of-day drug sensitivity emerges from the interaction of circadian, pharmacodynamic, growth, and experimental factors.

## Repository structure

### `Circadian_properties/`

Circadian-signal modeling and experimental parameter extraction:

- damped oscillator models;
- fitting of luminescence recordings;
- circadian amplitude, period, phase, and decay;
- simulation of circadian modulation of drug response.

### `Pharmacodynamics1/`

Effects of pharmacodynamic properties on ToD response:

- drug concentration;
- cytostatic and cytotoxic response models;
- maximal drug effect;
- dose-response sensitivity and correlation analyses.

### `Pharmacodynamics2/`

Additional drug properties, including:

- drug stability and preparation age;
- survival curves;
- time-dependent pharmacodynamic changes.

### `Cell_growth_dynamics/`

Analyses of cell growth, density, and normalization choices in temporal drug-response experiments.

### `Evaluation_time/`

Simulation and experimental analysis of how assay duration alters inferred IC50, Emax, Hill coefficient, and ToD-response magnitude.

### `Optimal_treatment_time/`

Model-based analyses of how circadian and pharmacodynamic parameters influence the treatment times associated with maximal predicted benefit.

## Computational methods

Methods represented in the repository include:

- nonlinear dynamical modeling of population growth;
- cytostatic and cytotoxic pharmacodynamic models;
- Hill-type dose-response modeling;
- nonlinear curve fitting;
- damped-oscillator fitting of circadian recordings;
- parameter sweeps and sensitivity analyses;
- Pearson/Spearman correlation analyses;
- simulation of time-of-day treatment-response profiles;
- quantitative analysis of assay-length and growth-rate effects.

## Requirements

Core Python packages used across the repository include:

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

The data generated and used for the published study are publicly available on Figshare:

**Dataset:** https://doi.org/10.6084/m9.figshare.28423709

The publication lists the primary code repository at:

https://github.com/Granada-Lab/drivers_ToD_drug_sensitivity

Some scripts in this personal repository retain the original research-directory structure and expect local input files such as those referenced under `Raw_Data/`. Paths must therefore be adjusted to reproduce individual analyses locally.

## Citation

> Gutu, N., Ishikuma, H., Ector, C. et al. **A combined mathematical and experimental approach reveals the drivers of time-of-day drug sensitivity in human cells.** _Communications Biology_ 8, 491 (2025). https://doi.org/10.1038/s42003-025-07931-1

## Contact

**Nica Gutu**  
Computational Biology / Data Science  
https://nica-gutu.github.io/website/
