
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spaco (Rcpp Implementation of SPACO)

<!-- badges: start -->

[![R-CMD-check](https://github.com/kexinwang3/spaco/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kexinwang3/spaco/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/kexinwang3/spaco/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/kexinwang3/spaco/actions/workflows/test-coverage.yaml)
[![lint](https://github.com/kexinwang3/spaco/actions/workflows/lint.yaml/badge.svg)](https://github.com/kexinwang3/spaco/actions/workflows/lint.yaml)
<!-- badges: end -->

## Description

`spaco` is an R package for modelling sparsely observed multivariate
longitudinal data using SPACO model (Guan, 2022). It provides tools for
initializing and estimating model parameters.

This package contains 3 main functions:

- `rank_selection`: performs rank selection via cross-validation.
- `train_prepare`: prepares user-level input data for model training.
- `train_spaco`: trains SPACO model on the prepared input data.

There are 4 other functions for case studies:

- `impact_data_wrangling`: transforms and maps raw data form of IMPACT
  into the desired format.
- `immune_data_wrangling`: transforms and maps raw data form of IMMUNE
  into the desired format.
- `feature_predict`: evaluates the predictive power of SPACO model.
- `feature_plot`: plots the observed versus estimated values for the
  selected feature.

The methodology background underlying the `spaco` package can be found
in the original paper. A Python implementation of the model is also
available on [GitHub](https://github.com/).

- [Reference](https://arxiv.org/abs/2104.05184): Guan, Leying. “Smooth
  and probabilistic PARAFAC model with auxiliary covariates.” arXiv
  preprint arXiv:2104.05184 (2022).
- [Python Implementation](https://github.com/LeyingGuan/SPACO): Guan,
  Leying. “SPACO (Python Implementation of SPACO).” GitHub repository:
  LeyingGuan/SPACO.

## Installation

A Conda environment is required to support Python modules in R
interface.

``` r
install.packages("reticulate")
reticulate::install_miniconda(force = TRUE)
reticulate::conda_install(envname = "spaco_env",
                          packages = c("numpy", "tensorly", "scikit-learn"))
```

You can install the development version of spaco from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("kexinwang3/spaco")
```

## Example

This is an example which shows you how to train SPACO model on `IMPACT`
dataset:

- Set up the Conda environment and load the library

``` r
reticulate::use_condaenv("spaco_env", required = TRUE)
library(spaco)
```

- Prepare `IMPACT` dataset

``` r
data("impact_imputed")
data("impact_missing")
impact <- impact_data_wrangling(impact_missing, impact_imputed)
```

- Perform rank selection

``` r
ranks <- c(2:10)
rank <- rank_selection(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
                       Z = impact$Z, ranks = ranks, early_stop = TRUE)
```

- Train SPACO model

``` r
spaco_object <- train_prepare(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
                              Z = impact$Z, K = rank, mean_removal = FALSE)
spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
                            tol = 1e-4, trace = TRUE)
```

- Evaluate predictive power

``` r
spaco_object <- feature_predict(spaco_object)
```

- Plot observed versus estimated values for selected features

``` r
feature_plot(spaco_object, "TcellsofLivecells", impact$columns_feature)
```

<img src="man/figures/README-impact-plot-1.png" width="100%" />

``` r
feature_plot(spaco_object, "TotalNeutrophilsofLivecells", impact$columns_feature)
```

<img src="man/figures/README-impact-plot-2.png" width="100%" />

``` r
feature_plot(spaco_object, "HLA.DR.ofTotalMono", impact$columns_feature)
```

<img src="man/figures/README-impact-plot-3.png" width="100%" />

``` r
feature_plot(spaco_object, "IL6", impact$columns_feature)
```

<img src="man/figures/README-impact-plot-4.png" width="100%" />
