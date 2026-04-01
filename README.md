# tb-erlang-progression-model
Mathematical model of tuberculosis progression comparing exponential and Erlang waiting-time assumptions, with applications to active case finding (ACF) interventions.
# TB Erlang Progression Model

## Overview

This repository contains code and data to reproduce the analyses presented in our study on tuberculosis (TB) disease progression. We compare conventional exponential (memoryless) models with multi-stage Erlang formulations to better capture biologically realistic distributions of asymptomatic disease duration.

The framework is used to evaluate how structural assumptions about disease progression influence long-term epidemiological projections and the estimated impact of active case finding (ACF) interventions.

---

## Key Features

* Implementation of both **exponential** and **Erlang (multi-stage)** TB progression models
* Bayesian calibration to epidemiological targets (India TB incidence)
* Simulation of **active case finding (ACF)** interventions
* Quantification of differences in projected incidence and cases averted
* Reproducible pipeline for model fitting and forward projection

---

## Repository Structure

```
.
├── data/                  # Input epidemiological data
├── src/                   # Model implementation and calibration
├── scripts/               # Scripts to run analyses
├── results/               # Output figures and posterior samples
├── README.md
└── requirements.txt
```

---

## Data Sources

Epidemiological data used for calibration were obtained from:

* World Health Organization (WHO) Global Tuberculosis Database

These data are publicly available and can be accessed at:
[https://worldhealthorg.shinyapps.io/tb_profiles/](https://worldhealthorg.shinyapps.io/tb_profiles/?_inputs_&tab=%22charts%22&lan=%22EN%22&iso3=%22IND%22&entity_type=%22country%22)

---

## Methods Summary

We construct a compartmental transmission model of TB with an explicit asymptomatic stage. Two alternative formulations are considered:

* **Exponential model**: assumes a memoryless progression process
* **Erlang model**: represents progression as a sequence of (k) stages with identical transition rates

The Erlang formulation preserves the same mean duration as the exponential model while modifying the shape of the dwell-time distribution, reducing the probability of unrealistically short durations.

Model parameters are estimated using Bayesian calibration with an adaptive Metropolis algorithm. Annual TB incidence targets are matched using likelihood-based inference.

---

## Reproducing Results

### 1. Install dependencies

```
pip install -r requirements.txt
```

### 2. Run calibration

```
python src/calibration.py
```

### 3. Run forward simulation

```
python src/acf_simulation.py
```

### 4. Generate figures

```
python scripts/plot_results.py
```

---

## Code Availability

All code used for model calibration, simulation, and analysis is provided in this repository.

---

## Citation

If you use this code or build upon this work, please cite:

*Author(s). Title. Journal (Year).*

---

## License

This project is licensed under the MIT License.

---

## Contact

For questions, please contact:
Zhichao Zhou
