# Run tests for wrapper.R

# rank_selection() -------------------------------------------------------------

test_that(
  "The rank_selection() returns a likelihood matrix when early_stop is true.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    ranks <- c(2:10)
    neglik <- suppressWarnings(rank_selection(X = impact$X, OBS = impact$OBS,
                                              T1 = impact$T1, Z = impact$Z,
                                              ranks = ranks, early_stop = TRUE,
                                              trace = FALSE))
    expect_true(inherits(neglik, "matrix"))
  }
)

test_that(
  "The rank_selection() expects errors when ranks contain integers
  less than or equal to 1.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    ranks <- c(1:10)
    expect_error(rank_selection(X = impact$X, OBS = impact$OBS,
                                T1 = impact$T1, Z = impact$Z,
                                ranks = ranks, early_stop = TRUE,
                                trace = FALSE))
  }
)



# train_prepare() --------------------------------------------------------------

test_that(
  "The train_prepare() expects warinings when the range of dfs for modeling the
  trajectories is not specified.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    expect_warning(train_prepare(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
                                 Z = impact$Z, K = 4, mean_removal = FALSE))
  }
)

test_that(
  "The train_prepare() runs without errors when run_prepare and run_initial
  are false.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    expect_no_error(train_prepare(X = impact$X, OBS = impact$OBS,
                                  T1 = impact$T1, Z = impact$Z, K = 4,
                                  lambda1_dfmax = 33, lambda1_dfmin = 1,
                                  mean_removal = FALSE, run_prepare = FALSE,
                                  run_initial = FALSE))
  }
)



# train_spaco() --------------------------------------------------------------

test_that(
  "The train_spaco() prints out the progress when trace is true.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    spaco_object <- train_prepare(X = impact$X, OBS = impact$OBS,
                                  T1 = impact$T1, Z = impact$Z, K = 4,
                                  lambda1_dfmax = 33, lambda1_dfmin = 1,
                                  mean_removal = FALSE)
    expect_message(train_spaco(spaco_object, max_iter = 30, min_iter = 1,
                               tol = 1e-4, trace = FALSE))
  }
)
