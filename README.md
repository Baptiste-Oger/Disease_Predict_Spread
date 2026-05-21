# Disease Prediction and Spread

This repository provides the R code associated with a statistical framework designed to predict local disease severity from imperfect indicators and reconstruct large-scale epidemiological dynamics over space and time.  

The framework combines:

- a stacked hurdle model based on multiple random forests for local-scale disease severity prediction;
- a semi-parametric spatio-temporal model used to reconstruct disease spread dynamics from either predictions or field observations.

In this repository, the methodology is applied to sugar beet yellows disease in France using satellite remote sensing data and field observations collected between 2019 and 2023.

### Data availability

The data are available on Zenodo at the following link:  

https://doi.org/10.5281/zenodo.20311758

Download the `Data.zip` archive and extract the `Data/` folder at the root of the R project.

### Project configuration

The `Config.R` script is used to configure file paths and project directories. It relies on the `here` package.

### Main scripts

The main scripts are located in the `Script/` folder.

- `Main_disease_predict_and_spread.R` is the main entry point of the framework. It runs:
  - the stacked hurdle model for local disease severity prediction from remote sensing indicators;
  - the semi-parametric spatio-temporal model for reconstructing disease spread dynamics over space and time.

The script relies on a collection of functions available in the `Script/Source/` subfolder.

### Repository structure

```text
.
├── Config.R
├── Data/
└── Script/
    ├── Main_disease_predict_and_spread.R
    └── Source/
```