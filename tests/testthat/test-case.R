# Run tests for case.R

# impact_data_wrangling() ------------------------------------------------------

test_that(
  "The impact_data_wrangling() returns a list.", {
    skip_if_no_modules()
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    expect_true(inherits(impact, "list"))
    }
  )



# immune_data_wrangling() ------------------------------------------------------

test_that(
  "The immune_data_wrangling() returns a list.", {
    skip_if_no_modules()
    data("immune_original")
    immune <- immune_data_wrangling(immune_original)
    expect_true(inherits(immune, "list"))
  }
)



# feature_predict() ------------------------------------------------------------

test_that(
  "The feature_predict() returns a list.", {
    skip_if_no_modules()
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    spaco_object <- train_prepare(X = impact$X, OBS = impact$OBS,
                                  T1 = impact$T1, Z = impact$Z, K = 4,
                                  lambda1_dfmax = 33, lambda1_dfmin = 1,
                                  mean_removal = FALSE)
    spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
                                tol = 1e-4, trace = FALSE)
    spaco_object <- feature_predict(spaco_object)
    expect_true(inherits(spaco_object, "list"))
  }
)



# feature_plot() ---------------------------------------------------------------

test_that(
  "The feature_plot() runs without errors.", {
    skip_if_no_modules()
    data("impact_imputed")
    data("impact_missing")
    impact <- impact_data_wrangling(impact_missing, impact_imputed)
    spaco_object <- train_prepare(X = impact$X, OBS = impact$OBS,
                                  T1 = impact$T1, Z = impact$Z, K = 4,
                                  lambda1_dfmax = 33, lambda1_dfmin = 1,
                                  mean_removal = FALSE)
    spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
                                tol = 1e-4, trace = FALSE)
    spaco_object <- feature_predict(spaco_object)
    expect_no_error(feature_plot(spaco_object, "TcellsofLivecells",
                                 impact$columns_feature))
  }
)
