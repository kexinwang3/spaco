# Wrapper for run SPACO starting from the input data

#' @title Rank Selection via Cross-validation
#' @description This function performs rank selection as suggested in SupCP
#' (Lock and Li, 2018). Model parameters are estimated on the training sets
#' and evaluated on the test sets. The rank is chosen to be the value that
#' achieves the lowest log-likelihood of test sets.
#' @param X An `I` (number of subjects) × `T` (number of times) × `J` (number
#' of features) array containing time-series data.
#' @param OBS An `I` (number of subjects) × `T` (number of times) × `J` (number
#' of features) array indicating whether an observation is available in array
#' `X`.
#' @param T1 A length `T` (number of times) vector containing measured time.
#' @param Z An `I` (number of subjects) × `q` (number of covariates) matrix
#' containing auxiliary covariates.
#' @param ranks An array of ranks containing integers greater than 1 from
#' which to select potential ranks.
#' @param early_stop Whether to stop the process when the cross-validated
#' marginal log-likelihood is no longer decreasing. Default is TRUE.
#' @param max_iter Maximum number of iterations in the model training
#' process. Default is 30.
#' @param cv_iter Maximum number of iterations in the cross-validated
#' training process. Default is 5.
#' @param nfolds Number of folds. Default is 5.
#' @param extra_std Coefficient of standard deviation when calculating the
#' standard deviation of log-likelihood. Default is 0.
#' @param random_state Control the randomness of each fold. When shuffling the
#' data before splitting into batches, it affects the ordering of the indices.
#' Default is 2022.
#' @param trace Whether to print the model training process. Default is TRUE.
#' @return An integer indicating selected rank.
#' @importFrom stats rnorm var
#' @examples
#' data("impact_imputed")
#' data("impact_missing")
#' impact <- impact_data_wrangling(impact_missing, impact_imputed)
#' ranks <- c(2:10)
#' rank <- rank_selection(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
#'                        Z = impact$Z, ranks = ranks, early_stop = TRUE)
#' @export
rank_selection <- function(X, OBS, T1, Z, ranks, early_stop = TRUE,
                           max_iter = 30, cv_iter = 5, nfolds = 5,
                           extra_std = 0, random_state = 2022, trace = TRUE) {
  if (!is.integer(ranks) || sum(ranks < 2)) {
    stop("Please specify a list of ranks containing integers greater than 1.")
  }
  neglik <- matrix(NA, nrow = dim(X)[1], ncol = length(ranks))
  means <- array(NA, dim = length(ranks))
  means_std <- array(NA, dim = length(ranks))
  for (k in seq(1, length(ranks))) {
    rank <- ranks[k]
    message(paste0("Rank ", rank, ": "))
    self <- train_prepare(X = X, OBS = OBS, T1 = T1, Z = Z, K = rank,
                          mean_removal = FALSE)
    self <- train_spaco(self, max_iter = max_iter, min_iter = 1, tol = 1e-4,
                        trace = trace)
    ids_split <- cut_foldid(n = dim(X)[1], nfolds = 5,
                            random_state = random_state)
    self <- cross_validation_train(self, ids_split$train_ids,
                                   ids_split$test_ids, max_iter = cv_iter,
                                   min_iter = 1, tol = 1e-3, trace = trace)
    self <- log_likelihood(self)
    neglik[, k] <- self$cross_likloss
    means[k] <- mean(neglik[, k])
    means_std[k] <- means[k] + extra_std *
      sqrt(var(neglik[, k]) / length(neglik[, k]) * (length(neglik[, k]) - 1)) /
      sqrt(dim(X)[1])
    if (early_stop) {
      if (k > 1) {
        if (means_std[k] > means[k - 1]) {
          break
        }
      }
    }
  }
  means <- colMeans(neglik)
  means_std <- means + sqrt(colMeans(sweep(neglik, 2, colMeans(neglik))^2)) /
    sqrt(dim(X)[1]) * 0.5
  means <- means[!is.na(means)]
  means_std <- means_std[!is.na(means_std)]
  idx_min <- which.min(means)
  rank <- as.integer(ranks[idx_min])
  return(rank)
}



train_prepare <- function(X, OBS, T1, Z, K = NULL,
                          run_prepare = TRUE, run_initial = TRUE,
                          random_number = 10, random_state = 0,
                          mean_removal = TRUE, homo_noise = TRUE,
                          fit_intercept = TRUE, initial_type = "pca",
                          update_sigma_noise = TRUE, update_sigma_mu = TRUE,
                          update_smooth_penalty = TRUE,
                          update_lasso_penalty = TRUE,
                          lambda1_dfmin = NULL, lambda1_dfmax = NULL,
                          Omega = NULL, h = NULL, nlambda1 = 10, nlambda2 = 100,
                          lambda1 = NULL, lambda2 = NULL, nlambda0 = 10,
                          foldid = NULL, betafolds = 10, kappa = 1e-2,
                          orthogonal = 0, ridge_traj = 1e-4, lam1criterion = 1,
                          lasso_max_iter = 1e4, lasso_tol = 1e-8, eps = 1e-10) {
  self <- list("X" = X, "OBS" = OBS, "T1" = T1, "Z" = Z, "K" = K,
               "random_number" = random_number, "random_state" = random_state,
               "homo_noise" = homo_noise, "fit_intercept" = fit_intercept,
               "initial_type" = initial_type, "foldid" = foldid,
               "kappa" = kappa, "Omega" = Omega, "h" = h,
               "update_smooth_penalty" = update_smooth_penalty,
               "update_lasso_penalty" = update_lasso_penalty,
               "update_sigma_noise" = update_sigma_noise,
               "update_sigma_mu" = update_sigma_mu,
               "lambda1_dfmin" = lambda1_dfmin, "lambda1_dfmax" = lambda1_dfmax,
               "nlambda1" = nlambda1, "nlambda2" = nlambda2,
               "nlambda0" = nlambda0, "eps" = eps,
               "lambda1" = lambda1, "lambda2" = lambda2,
               "orthogonal" = orthogonal, "ridge_traj" = ridge_traj,
               "lam1criterion" = lam1criterion, "betafolds" = betafolds,
               "lasso_max_iter" = lasso_max_iter, "lasso_tol" = lasso_tol)
  self$num_subjects <- dim(self$X)[1]
  self$num_times <- dim(self$X)[2]
  self$num_features <- dim(self$X)[3]
  if (run_prepare) {
    message("Start data preparation: ")
    data_prepare <- prepare(X, OBS, T1, mean_removal = mean_removal,
                            nlam = nlambda0, lams = NULL, kappa = kappa)
    Omega <- data_prepare$Omega
    h <- data_prepare$h
    R <- data_prepare$R
    message("Data preparation done.")
  }
  if (run_initial) {
    message("Start parameter initialization: ")
    self <- initialization(self, X = R, OBS = OBS, T1 = T1,
                           V_norm = 1, Phi_norm = sqrt(h),
                           eps = 1 / sqrt(self$num_features * self$num_times))
    set.seed(self$random_state)
    self <- V_init(self)
    self <- Phi_init(self)
    self <- UG_init(self)
    self <- G_PARAFAC(self)
    self <- Var_init(self)
    message("Parameter initialization done.")
  }
  if (is.null(self$Omega)) {
    self$Omega <- Omega
  }
  if (is.null(self$h)) {
    self$h <- h
  }
  if (!is.null(self$update_smooth_penalty)) {
    if (is.null(self$lambda1_dfmax) || is.null(self$lambda1_dfmin)) {
      warning(paste0("Did not specify the range of dfs (lambda1_dfmin or ",
                     "lambda1_dfmax) for modeling the trajectories."))
    }
    self$lambda1_dfmax <- self$num_times - 2
    self$lambda1_dfmin <- 1
  }
  self <- check_init(self)
  return(self)
}



train_spaco <- function(self, max_iter = 100, min_iter = 1, tol = 1e-6,
                        update_sigma_mu = TRUE, update_sigma_noise = TRUE,
                        update_lasso_penalty = TRUE,
                        update_smooth_penalty = TRUE,
                        trace = FALSE, reorder = TRUE) {
  error <- tol * 2
  iter <- 0
  self$update_lasso_penalty <- update_lasso_penalty
  self$update_smooth_penalty <- update_smooth_penalty
  self <- posteriorU(self)
  message("Start parameter update: ")
  while (iter < max_iter && error > tol) {
    if (trace) {
      print(paste0("Iteration ", iter, ": ", error))
    }
    self <- Phi_update(self)
    Phierror <- self$Phierror
    self <- V_update(self)
    Verror <- self$Verror
    self <- posteriorU(self)
    if (!is.null(self$Z)) {
      self <- beta_update(self)
    }
    if (update_sigma_mu) {
      self <- sigma_mu_update(self)
    }
    if (update_sigma_noise) {
      self <- sigma_noise_update(self)
    }
    self <- reordering(self)
    error <- Phierror + Verror
    iter <- iter + 1
  }
  message("Parameter update done.")
  return(self)
}
