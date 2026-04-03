# tb-erlang-progression-model

Mathematical model of tuberculosis(TB) progression comparing exponential and Gamma-distributed waiting-time assumptions, with applications to active case finding (ACF) interventions.

---

## Overview

This repository contains code and data to reproduce the analyses presented in our study on TB disease progression. We compare conventional exponential (memoryless) models with multi-stage Gamma-distributed formulations to better capture biologically realistic distributions of asymptomatic disease duration.

The framework is used to evaluate how structural assumptions about disease progression influence long-term epidemiological projections and the estimated impact of ACF interventions.

---

## Key Features

* Implementation of both **exponential** and **Gamma-distributed (multi-stage)** TB progression models
* Bayesian calibration to epidemiological targets (India TB incidence)
* Simulation of **active case finding (ACF)** interventions
* Model comparison using **WAIC** and **LOO-CV**
* Reproducible pipeline for model fitting and forward projection

---

## Repository Structure
```
.
├── data/
│   ├── who_data/
│   └── nti_cohort_data/
│
├── src/
│   ├── calibration/
│   ├── evaluation/
│   ├── models/
│   └── visualization/
│
├── results/
│   ├── figures/
│   ├── tables/
│   └── logs/
│       ├── acf_sensitivity/       
│       └── calibration/
│           ├── India/
│           └── nti/
│
└── README.md
```
---

## Data Sources

Epidemiological data used for calibration were obtained from the following sources:

* **World Health Organization (WHO) Global Tuberculosis Database**
  Publicly available at:
  https://worldhealthorg.shinyapps.io/tb_profiles/

* **NTI cohort data**
  Derived from published cohort studies of tuberculosis progression and/or publicly available datasets.
  Where applicable, data were obtained from publicly accessible repositories ([e.g., GitHub](https://github.com/ERC-TBornotTB/Inf-Dis)) and processed for analysis.

All data sources are publicly available or derived from previously published studies.

---

## Methods Summary

We construct a compartmental transmission model of TB with an explicit asymptomatic stage. Two alternative formulations are considered:

* **Exponential model**: assumes a memoryless progression process
* **Erlang model**: represents progression as a sequence of (k) stages with identical transition rates

The Gamma-distributed formulation preserves the same mean duration as the exponential model while modifying the shape of the dwell-time distribution, reducing the probability of unrealistically short durations.

Model parameters are estimated using Bayesian calibration with an adaptive Metropolis algorithm. Model comparison is conducted using WAIC and leave-one-out cross-validation (LOO-CV).

---

## Reproducing Results (MATLAB)

All analyses were conducted in MATLAB.

### 1. Simple model calibration

Run:

```matlab
bayes_calibration_simple.m
bayes_calibration_multi.m
```
## Model Components

The modeling framework is organized into modular MATLAB functions:

### ODE Models

* `simpleTBmodel.m`
  Implements the exponential (memoryless) TB progression model.

* `erlangODE.m`
  Implements the Gamma-distributed (multi-stage) progression model with (k) subclinical compartments.

---


### Likelihood Functions

* `logPosterior.m`
  Log-posterior for the Gamma-distributed model.

* `logPosterior_simple.m`
  Log-posterior for the simple(exponential) model.

---

### Diagnostics

* `modelResidual.m`
  Computes residuals between model predictions and observed data.

---

### Notes

All functions are located in the `src/` directory and are automatically called by the main calibration scripts. No manual function linking is required.

### 2. Model comparison

Compute WAIC and LOO-CV:

```matlab
model_evaluation_simple.m
model_evaluation.m
```

### 3. Complex (Extended) model calibration

Run:

```matlab
complex_calibration.m
```

### 4. Forward projection and figures

Run:

```matlab
UI.m
IndianPrediction_2trajectory.m
IndianPrediction_data_full.m
IndianPrediction_efsimple.m
```

---

## Code Availability

All MATLAB code used for model calibration, evaluation, and simulation is provided in this repository.

---


## License

This project is licensed under the MIT License.

---

## Contact

Zhichao Zhou
