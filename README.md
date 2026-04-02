# tb-erlang-progression-model

Mathematical model of tuberculosis progression comparing exponential and Erlang waiting-time assumptions, with applications to active case finding (ACF) interventions.

---

## Overview

This repository contains code and data to reproduce the analyses presented in our study on tuberculosis (TB) disease progression. We compare conventional exponential (memoryless) models with multi-stage Erlang formulations to better capture biologically realistic distributions of asymptomatic disease duration.

The framework is used to evaluate how structural assumptions about disease progression influence long-term epidemiological projections and the estimated impact of active case finding (ACF) interventions.

---

## Key Features

* Implementation of both **exponential** and **Erlang (multi-stage)** TB progression models
* Bayesian calibration to epidemiological targets (India TB incidence)
* Simulation of **active case finding (ACF)** interventions
* Model comparison using **WAIC** and **LOO-CV**
* Reproducible pipeline for model fitting and forward projection

---

## Repository Structure

```
.
├── data/                  # Input epidemiological data
├── src/                   # MATLAB model and calibration code
├── evaluation/            # Model comparison
├── results/               # Output figures and posterior samples
├── README.md
```

---

## Data Sources

Epidemiological data used for calibration were obtained from:

* World Health Organization (WHO) Global Tuberculosis Database

Publicly available at:
https://worldhealthorg.shinyapps.io/tb_profiles/

---

## Methods Summary

We construct a compartmental transmission model of TB with an explicit asymptomatic stage. Two alternative formulations are considered:

* **Exponential model**: assumes a memoryless progression process
* **Erlang model**: represents progression as a sequence of (k) stages with identical transition rates

The Erlang formulation preserves the same mean duration as the exponential model while modifying the shape of the dwell-time distribution, reducing the probability of unrealistically short durations.

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
  Implements the Erlang (multi-stage) progression model with (k) subclinical compartments.

---


### Likelihood Functions

* `logPosterior.m`
  Log-posterior for the Erlang model.

* `logPosterior_simple.m`
  Log-posterior for the exponential model.

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

### 3. Complex (Erlang) model calibration

Run:

```matlab
complex_calibration.m
```

### 4. Forward projection and figures

Run:

```matlab
testbaseUI.m
```

---

## Code Availability

All MATLAB code used for model calibration, evaluation, and simulation is provided in this repository.

---

## Citation

If you use this code or build upon this work, please cite:

*Author(s). Title. Journal (Year).*

---

## License

This project is licensed under the MIT License.

---

## Contact

Zhichao Zhou
