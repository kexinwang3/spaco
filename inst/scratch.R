# IMPACT

# data preparation
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
spaco_object <- feature_predict(spaco_object)
feature_plot(spaco_object, "TcellsofLivecells", impact$columns_feature)
feature_plot(spaco_object, "TotalNeutrophilsofLivecells",
             impact$columns_feature)
feature_plot(spaco_object, "HLA.DR.ofTotalMono", impact$columns_feature)
feature_plot(spaco_object, "IL6", impact$columns_feature)



# IMMUNE

# data preparation
data("immune_original")
immune <- immune_data_wrangling(immune_original)

# rank selection
ranks <- c(2:10)
rank <- rank_selection(X = immune$X, OBS = immune$OBS, T1 = immune$T1,
                       Z = immune$Z, ranks = ranks, early_stop = TRUE)
# model training
spaco_object <- train_prepare(X = immune$X, OBS = immune$OBS, T1 = immune$T1,
                              Z = immune$Z, K = rank, mean_removal = FALSE)
spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
                            tol = 1e-4, trace = TRUE)
spaco_object <- feature_predict(spaco_object)
feature_plot(spaco_object, "IgG..LY", immune$columns_feature)
feature_plot(spaco_object, "CD19._CD20..LY.1", immune$columns_feature)
feature_plot(spaco_object, "CD4._CD8..CD3.", immune$columns_feature)
feature_plot(spaco_object, "LY.All_CD45", immune$columns_feature)
