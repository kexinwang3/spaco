data("impact_imputed")
data("impact_missing")
impact <- impact_data_wrangling(impact_missing, impact_imputed)

# rank selection
ranks <- c(2:10)
rank <- rank_selection(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
                         Z = impact$Z, ranks = ranks, early_stop = TRUE)

# model training
spaco_object <- train_prepare(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
                              Z = impact$Z, K = rank, mean_removal = FALSE)
spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
                            tol = 1e-4, trace = TRUE)
spaco_object <- impact_predict(spaco_object)
impact_plot(spaco_object, "TcellsofLivecells",
                        impact$imputed_pt, impact$filtered_feature_idx)
