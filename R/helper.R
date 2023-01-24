# Helpers for initialization and estimation

tnsr_unfold <- function(tnsr) {
  # unfold tensor in the subject(I) dimension
  tnsr_mode1 <- tnsr[, , 1]
  for (i in 2:dim(tnsr)[3]) {
    tnsr_mode1 <- cbind(tnsr_mode1, tnsr[, , i])
  }

  # unfold tensor in the time(T) dimension
  tnsr_mode2 <- t(tnsr[, , 1])
  for (i in 2:dim(tnsr)[3]) {
    tnsr_mode2 <- cbind(tnsr_mode2, t(tnsr[, , i]))
  }

  # unfold tensor in the feature(J) dimension
  tnsr_mode3 <- t(tnsr[, 1, ])
  for (i in 2:dim(tnsr)[2]) {
    tnsr_mode3 <- cbind(tnsr_mode3, t(tnsr[, i, ]))
  }
  multi_return <- function() {
    return_list <- list("tnsr_mode1" = tnsr_mode1, "tnsr_mode2" = tnsr_mode2,
                        "tnsr_mode3" = tnsr_mode3)
    return(return_list)
  }
  mode <- multi_return()
  return(mode)
}



#' @importFrom stats lm
#' @importFrom pracma Reshape
smooth_FPCA <- function(W_total, h, num_times, num_subjects) {
  W_smooth <- matrix(0, nrow = num_times, ncol = num_times)
  location <- location_search(num_times, num_subjects, W_total)
  y <- location$y
  y <- array(y)
  row_location <- location$row_location
  row_location <- array(row_location) + 1
  column_location <- location$column_location
  column_location <- array(column_location) + 1
  idx1 <- list()
  for (i in seq(1, length(row_location))) {
    a <- unlist(row_location[i])
    b <- unlist(column_location[i])
    if (a != b) {
      idx1 <- append(idx1, i)
    }
  }

  # off-diagonal
  idx1 <- unlist(idx1)
  row_location10 <- row_location[idx1]
  column_location10 <- column_location[idx1]
  y10 <- y[idx1]
  y10 <- array(unlist(y10))
  for (s in seq(1, num_times)) {
    for (t in seq(1, num_times)) {
      if (s < t) {
        x10 <- t(rbind((unlist(row_location10) - s),
                       (unlist(column_location10) - t)))
        weights1 <- exp((- x10[, 1]^2 - x10[, 2]^2) / (2 * h))
        W_smooth[s, t] <- lm(y10 ~ x10, weights = weights1)$coef[1]
        W_smooth[t, s] <- lm(y10 ~ x10, weights = weights1)$coef[1]
      }
    }
  }

  # diagonal
  idx2 <- list()
  for (i in seq(1, length(row_location))) {
    a <- unlist(row_location[i])
    b <- unlist(column_location[i])
    if (a == b) {
      idx2 <- append(idx2, i)
    }
  }
  idx2 <- unlist(idx2)
  row_location20 <- row_location[idx2]
  y20 <- y[idx2]
  y20 <- array(unlist(y20))
  for (t in seq(1, num_times)) {
    x20 <- unlist(row_location20) - t
    weights2 <- exp(- x20^2 / (2 * h))
    x20 <- Reshape(x20, length(x20), 1)
    W_smooth[t, t] <- lm(y20 ~ x20, weights = weights2)$coef[1]
  }

  return(W_smooth)
}



G_PARAFAC <- function(self) {
  # perform PARAFAC on G
  set.seed(self$random_state)
  G <- self$G

  # svd for matrix initialization
  Gparafac_svd <- parafac_svd(G, self$K, self$random_state)
  V <- matrix(0, nrow = dim(self$V0)[1], ncol = dim(self$V0)[2])
  Phi <- matrix(0, nrow = dim(self$Phi0)[1], ncol = dim(self$Phi0)[2])
  VPhi <- matrix(0, nrow = dim(self$V0)[1] * dim(self$Phi0)[1],
                 ncol = dim(self$Phi0)[2])
  B_svd <- Gparafac_svd[[2]]
  C_svd <- Gparafac_svd[[3]]
  VPhi_svd <- kronecker_VPhi(self$K, self$V0, self$Phi0,
                             V, Phi, VPhi, B_svd, C_svd)
  # set V, Phi, VPhi
  V <- VPhi_svd$V
  Phi <- VPhi_svd$Phi
  VPhi <- VPhi_svd$VPhi

  xhat <- matrix(0, nrow = dim(self$U0)[1],
                 ncol = dim(self$V0)[1] * dim(self$Phi0)[1])
  rescale_U_svd <- rescale_U(self$num_subjects, self$num_times,
                             self$num_features, self$K, self$Xmodes$XI,
                             self$Omodes$OI, VPhi, self$eps)
  xhat <- rescale_U_svd$xhat
  error_svd <- mean((self$Xmodes$XI[self$Omodes$OI == 1] -
                       xhat[self$Omodes$OI == 1])^2)

  if (self$random_number > 0) {
    errors_random <- array(0, self$random_number)
    seeds <- array(0, self$random_number)
    for (i in seq(1, self$random_number)) {
      seeds[i] <- as.integer(self$random_number + self$random_state)
    }
    for (j in seq(1, self$random_number)) {
      set.seed(seeds[j])
      # random for matrix initialization
      Gparafac_random <- parafac_random(G, self$K, seeds[j])
      V <- matrix(0, nrow = dim(self$V0)[1], ncol = dim(self$V0)[2])
      Phi <- matrix(0, nrow = dim(self$Phi0)[1], ncol = dim(self$Phi0)[2])
      VPhi <- matrix(0, nrow = dim(self$V0)[1] * dim(self$Phi0)[1],
                     ncol = dim(self$Phi0)[2])
      B_random <- Gparafac_random[[2]]
      C_random <- Gparafac_random[[3]]
      VPhi_random <- kronecker_VPhi(self$K, self$V0, self$Phi0,
                                    V, Phi, VPhi, B_random, C_random)
      V <- VPhi_random$V
      Phi <- VPhi_random$Phi
      VPhi <- VPhi_random$VPhi
      xhat <- matrix(0, nrow = dim(self$U0)[1],
                     ncol = dim(self$V0)[1] * dim(self$Phi0)[1])
      rescale_U_random <- rescale_U(self$num_subjects, self$num_times,
                                    self$num_features, self$K, self$Xmodes$XI,
                                    self$Omodes$OI, VPhi, self$eps)
      xhat <- rescale_U_random$xhat
      errors_random[j] <- mean((self$Xmodes$XI[self$Omodes$OI == 1] -
                                  xhat[self$Omodes$OI == 1])^2)
    }
    error_random <- min(errors_random)
    idx <- which(errors_random == error_random)[1]
    Gparafac_random <- parafac_random(G, self$K, seeds[idx])
  } else {
    error_random <- Inf
    Gparafac_random <- Gparafac_svd
  }

  if (error_random < error_svd) {
    Gparafac <- Gparafac_random
  } else {
    Gparafac <- Gparafac_svd
  }

  # stack as the columns of A, B, C
  A <- Gparafac[[1]]
  B <- Gparafac[[2]]
  C <- Gparafac[[3]]

  # set U, Phi, V
  self$U <- self$U0 %*% A
  self$Phi <- self$Phi0 %*% B
  self$V <- self$V0 %*% C

  # rescale Phi, V
  self$VPhi <- matrix(0, nrow = dim(self$V0)[1] * dim(self$Phi0)[1],
                      ncol = self$K)
  rescale_VPhi_list <- rescale_VPhi(self$K, self$Phi_norm, self$V_norm, self$V,
                                    self$Phi, self$VPhi)
  self$V <- rescale_VPhi_list$V
  self$Phi <- rescale_VPhi_list$Phi
  self$VPhi <- rescale_VPhi_list$VPhi

  # rescale U
  self$xhat <- matrix(0, nrow = dim(self$U0)[1],
                      ncol = dim(self$V0)[1] * dim(self$Phi0)[1])
  self$U <- matrix(0, self$num_subjects, self$K)
  rescale_U_list <- rescale_U(self$num_subjects, self$num_times,
                              self$num_features, self$K, self$Xmodes$XI,
                              self$Omodes$OI, self$VPhi, self$eps)
  self$U <- rescale_U_list$U
  self$xhat <- rescale_U_list$xhat
  return(self)
}



#' @importFrom glmnet glmnet
#' @importFrom stats coef
my_glmnet <- function(x, y, nlambda, lambda, penalty_factor) {
  dim(y) <- c(length(y), 1)
  fit <- glmnet(x, y, family = c("gaussian"), nlambda = nlambda,
                intercept = FALSE, penalty.factor = penalty_factor)
  coefficients0 <- coef(fit, s = lambda)
  summary.coef <- summary(coefficients0)
  coefficients0 <- array(summary.coef)
  coefficients <- array(0, dim(x)[2])
  if (is.na(coefficients0[[1]][3])) {
    coefficients <- array(0, dim(x)[2])
  } else {
    for (i in seq(1, dim(coefficients0)[1])) {
      if (as.integer(coefficients0[[i]][1]) - 1 > 0) {
        coefficients[as.integer(coefficients0[[i]][1]) - 1] <-
          coefficients0[[i]][3]
      } else {
        coefficients[dim(x)[2] - abs(as.integer(coefficients0[[i]][1]) - 1)] <-
          coefficients0[[i]][3]
      }
    }
  }
  return(coefficients)
}



#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
my_cv_glmnet <- function(x, y, nlambda, penalty_factor, nfolds, foldid) {
  dim(y) <- c(length(y), 1)
  cvfit <- cv.glmnet(x, y, nlambda = nlambda, intercept = FALSE,
                     nfolds = nfolds, foldid = foldid)
  cvm <- cvfit$cvm
  lambdas <- cvfit$lambda
  cvsd <- cvfit$cvsd
  lambda.min <- cvfit$lambda.min
  cv_coefficients0 <- coef(cvfit, s = "lambda.min")
  summary.cvcoef <- summary(cv_coefficients0)
  cv_coefficients0 <- array(summary.cvcoef)
  cv_coefficients <- array(0, dim(x)[2])
  if (is.na(cv_coefficients0[[1]][3])) {
    cv_coefficients <- array(0, dim(x)[2])
  } else {
    for (i in seq(1, dim(cv_coefficients0)[1])) {
      if (as.integer(cv_coefficients0[[i]][1]) - 1 > 0) {
        cv_coefficients[as.integer(cv_coefficients0[[i]][1]) - 1] <-
          cv_coefficients0[[i]][3]
      } else {
        cv_coefficients[dim(x)[2] -
                          abs(as.integer(cv_coefficients0[[i]][1]) - 1)] <-
          cv_coefficients0[[i]][3]
      }
    }
  }
  multi_return <- function() {
    return_list <- list("coefficients" = cv_coefficients, "lambdas" = lambdas,
                        "cvm" = cvm, "cvsd" = cvsd, "lambda.min" = lambda.min)
    return(return_list)
  }
  cv <- multi_return()
  return(cv)
}



subject_cv_mean <- function(X, Psi, O1, s0, Omega, lams) {
  # remove the mean time curve if needed
  nlam <- length(lams)
  cv <- array(0, dim = c(dim(X)[1], dim(X)[3], nlam))
  cv <- cv_creator(nlam, dim(X)[1], dim(X)[2], dim(X)[3], lams, s0 - 1, O1, Psi,
                   Omega, X, cv)
  cv <- colMeans(colMeans(cv))
  idx <- which(cv == min(cv))
  lam <- lams[idx]
  mean_curve <- mean_curve_update(dim(X)[1], dim(X)[2], dim(X)[3],
                                  lam, Psi, O1, Omega, X)
  B0 <- mean_curve$B0
  multi_return <- function() {
    return_list <- list("cv" = cv, "B0" = B0)
    return(return_list)
  }
  cv_mean <- multi_return()
  return(cv_mean)
}



parafac_svd <- function(G, K, random_state) {
  K <- as.integer(K)
  random_state <- as.integer(random_state)
  numpy$random$sample(random_state)
  # svd for matrix initialization
  Gparafac <- tensorly$decomposition$parafac(tensor = G, rank = K,
                                             n_iter_max = as.integer(10000),
                                             init = "svd",
                                             random_state = random_state)[1]
  return(Gparafac)
}



parafac_random <- function(G, K, random_state) {
  K <- as.integer(K)
  random_state <- as.integer(random_state)
  numpy$random$sample(random_state)
  # random for matrix initialization
  Gparafac <- tensorly$decomposition$parafac(tensor = G, rank = K,
                                             n_iter_max = as.integer(100),
                                             init = "random",
                                             random_state = random_state)[1]
  return(Gparafac)
}



cut_foldid <- function(n, nfolds, random_state) {
  n <- as.integer(n)
  nfolds <- as.integer(nfolds)
  random_state <- as.integer(random_state)
  k_fold <- sklearn$model_selection$KFold(n_splits = nfolds, shuffle = TRUE,
                                          random_state = random_state)
  subject_ids <- seq(1, n)
  split_id_obj <- k_fold$split(X = subject_ids)
  train_ids <- list()
  test_ids <- list()
  index <- py$list(py$enumerate(split_id_obj))
  for (i in seq(1, nfolds)) {
    train_ids[[i]] <- index[[i]][[2]][[1]] + 1
    test_ids[[i]] <- index[[i]][[2]][[2]] + 1
  }
  multi_return <- function() {
    data_split <- list("train_ids" = train_ids, "test_ids" = test_ids)
    return(data_split)
  }
  data_split <- multi_return()
  return(data_split)
}



#' @importFrom testthat skip
#' @importFrom reticulate py_module_available
skip_if_no_modules <- function() {
  # skip tests if we don't have the required modules
  have_numpy <- reticulate::py_module_available("numpy")
  if (!have_numpy) {
    skip("numpy not available for testing")
  }
  have_tensorly <- py_module_available("tensorly")
  if (!have_tensorly) {
    skip("tensorly not available for testing")
  }
  have_sklearn <- py_module_available("sklearn")
  if (!have_sklearn) {
    skip("sklearn not available for testing")
  }
}
