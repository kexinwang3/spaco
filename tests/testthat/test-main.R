# Run tests for main.R

# V_init() ---------------------------------------------------------------------

test_that(
  "The V_init() runs without errors when initialization type is random.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    expect_no_error(train_prepare(X = impact$X, OBS = impact$OBS,
                                  T1 = impact$T1, Z = impact$Z, K = 4,
                                  initial_type = "random", mean_removal = FALSE,
                                  lambda1_dfmax = 33, lambda1_dfmin = 1))
  }
)



# Phi_init() -------------------------------------------------------------------

test_that(
  "The Phi_init() runs without errors when initialization type is random.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    expect_no_error(train_prepare(X = impact$X, OBS = impact$OBS,
                                  T1 = impact$T1, Z = impact$Z, K = 4,
                                  initial_type = "random", mean_removal = FALSE,
                                  lambda1_dfmax = 33, lambda1_dfmin = 1))
  }
)



# beta_update() ----------------------------------------------------------------

test_that(
  "The beta_update() runs without errors when update_lasso_penalty is false.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    expect_no_error(train_prepare(X = impact$X, OBS = impact$OBS,
                                  T1 = impact$T1, Z = impact$Z, K = 4,
                                  mean_removal = FALSE,
                                  update_lasso_penalty = FALSE,
                                  lambda1_dfmax = 33, lambda1_dfmin = 1))
  }
)



# prepare() --------------------------------------------------------------------

test_that(
  "The prepare() runs without errors when mean_removal is true.", {
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    expect_no_error(train_prepare(X = impact$X, OBS = impact$OBS,
                                  T1 = impact$T1, Z = impact$Z, K = 4,
                                  mean_removal = TRUE,
                                  lambda1_dfmax = 33, lambda1_dfmin = 1))
  }
)
