# Initialization and estimation of model parameters

#' @importFrom stats rnorm
V_init <- function(self) {
  set.seed(self$random_state)
  K0 <- self$K
  if (K0 > dim(self$X)[3]) {
    K0 <- dim(self$X)[3]
  }
  if (self$initial_type == "pca") {
    sums <- colSums(1.0 - self$Omodes$OJ)
    idx <- which(sums == 0)
    tmp <- self$Xmodes$XJ[, idx]
    scaled_tmp <- sweep(tmp, 2, colMeans(tmp))
    svd_tmp <- svd(scaled_tmp)
    tmp <- matrix(unlist(svd_tmp[2]), nrow = dim(scaled_tmp)[1])
    self$V0 <- tmp[, 1:K0]
  } else {
    # random initialization
    self$V0 <- matrix(rnorm(dim(self$X)[3] * self$K), dim(self$X)[3], self$K)
    for (i in 1:self$K) {
      self$V0[, i] <- self$V0[, i] / sqrt(sum(self$V0[, i]^2))
    }
  }
  return(self)
}



#' @importFrom stats rnorm
Phi_init <- function(self) {
  set.seed(self$random_state)
  K0 <- self$K
  if (K0 > dim(self$X)[2]) {
    K0 <- dim(self$X)[2]
  }
  if (self$initial_type == "pca") {
    # projection along the kth column of V0
    Y <- array(0, dim = c(self$num_subjects, self$num_times, K0))
    Y <- V_projection(K0, self$num_times, self$num_subjects, self$X, self$V0)

    # empirical estimate of the total product matrix
    W_total <- array(0, dim = c(self$num_times, self$num_times,
                                self$num_subjects))
    W_total <- empirical_W_total(self$num_times, self$num_subjects,
                                 K0, self$OBS, Y)

    # number of non-zero pairs in Wcount
    count <- 0
    count <- pairs_count(self$num_times, self$num_subjects, W_total)

    if (is.null(self$h)) {
      h <- self$num_times / sqrt(count)
    } else {
      h <- self$h
    }
    W_smooth <- smooth_FPCA(W_total, h, self$num_times, self$num_subjects)
    # return the top K left singular vectors
    svd_W_smooth <- svd(W_smooth)
    tmp <- matrix(unlist(svd_W_smooth[2]), nrow = dim(W_smooth)[1],
                  ncol = dim(W_smooth)[2])
    self$Phi0 <- tmp[, 1:K0]
  } else {
    # random initialization
    self$Phi0 <- matrix(rnorm(dim(self$X)[2] * K0), dim(self$X)[2], K0)
    for (i in 1:K0) {
      self$Phi0[, i] <- self$Phi0[, i] / sqrt(sum(self$Phi0[, i]^2))
    }
  }
  return(self)
}



UG_init <- function(self) {
  # kronecker product of V0 & Phi0
  set.seed(self$random_state)
  V0Phi0 <- kronecker(self$V0, self$Phi0)

  # set the default value of eps
  if (is.null(self$eps)) {
    self$eps <- 1 / sqrt(self$num_times * self$num_features)
  }

  # run a ridge-penalized regression of rows of Xi on VPhi
  Utmp <- matrix(0, self$num_subjects, self$K * self$K)
  Utmp <- ridge_penalized(self$num_subjects, self$num_times, self$num_features,
                 self$K, self$Xmodes$XI, self$Omodes$OI, V0Phi0, self$eps)

  # return the top K left singular vectors
  svd_tmp1 <- svd(Utmp)
  tmp1 <- matrix(unlist(svd_tmp1[2]), nrow = dim(Utmp)[1], ncol = dim(Utmp)[2])
  self$U0 <- tmp1[, 1:self$K]

  # estimate G from the regression coefficients
  self$G <- array(0, dim = c(self$K, self$K, self$K))
  self$G <- rearrange_G(self$K, Utmp, self$U0, self$G)
  return(self)
}



Var_init <- function(self) {
  # initialize variances for uk
  self$sigma_mu <- colMeans(self$U ^ 2)
  self$sigma_mu <- array(self$sigma_mu)
  idx <- which(self$Omodes$OI == 1, arr.ind = TRUE) - 1
  Osum <- sum(self$Omodes$OI)
  if (self$homo_noise) {
    # initialize variances for Xj
    self$sigma_noise <- sigma_init(self$homo_noise, self$num_features, Osum,
                                   self$xhat, self$Xmodes$XI, idx)
    self$sigma_noise <- array(self$sigma_noise)
  }
  return(self)
}



#' @importFrom stats rnorm
check_init <- function(self) {
  self$intercepts <- array(0, self$K)
  self$vec <- matrix(0, nrow = self$num_subjects, ncol = self$K)
  self$mat <- array(0, dim = c(self$num_subjects, self$K, self$K))
  if (is.null(self$h)) {
    self$h <- self$num_times
  }
  if (is.null(self$V)) {
    self$V <- matrix(rnorm(self$num_features * self$K),
                     self$num_features, self$K)
  }
  for (i in 1:self$K) {
    self$V[, i] <- self$V[, i] / sqrt(sum(self$V[, i]^2))
  }
  if (is.null(self$Phi)) {
    self$Phi <- matrix(rnorm(self$num_times * self$K), self$num_times, self$K)
  }
  for (i in 1:self$K) {
    self$Phi[, i] <- self$Phi[, i] / sqrt(sum(self$Phi[, i]^2) / self$h)
  }
  if (!is.null(self$Z)) {
    self$beta <- matrix(0, nrow = dim(self$Z)[2], ncol = self$K)
  }
  if (is.null(self$mu)) {
    self$mu <- matrix(0, nrow = self$num_subjects, ncol = self$K)
  }
  if (is.null(self$cov)) {
    self$cov <- array(0, dim = c(self$num_subjects, self$K, self$K))
  }
  if (is.null(self$sigma_mu)) {
    self$sigma_mu <- array(1, dim = self$K)
  }
  if (is.null(self$lambda1)) {
    self$lambda1 <- array(1 * 1e-2, dim = self$K)
  }
  if (is.null(self$lambda2)) {
    self$lambda2 <- array(1 * 1e-2, dim = self$K)
  }
  if (is.null(self$sigma_noise)) {
    self$sigma_noise <- array(1, dim = self$num_features)
  }
  if (is.null(self$orthogonal)) {
    self$orthogonal <- 2 * sqrt(self$num_subjects) *
      log(self$num_subjects * self$num_times * self$num_subjects)
  }
  return(self)
}



posteriorU <- function(self) {
  # kronecker product of V & Phi
  VPhi <- matrix(0, self$num_times * self$num_features, self$K)
  VPhi <- posterior_kronecker_VPhi(self$K, self$V, self$Phi, VPhi)

  s <- array(0, self$num_times * self$num_features)
  s <- sigma_transform(self$num_times, self$num_features, self$sigma_noise, s)

  # posterior mean and covariance of U
  self$vec <- matrix(0, nrow = self$num_subjects, ncol = self$K)
  self$mat <- array(0, dim = c(self$num_subjects, self$K, self$K))
  posterior <- posterior_distribution(self$num_subjects, self$num_times,
                                      self$num_features, self$K,
                                      self$fit_intercept, self$intercepts,
                                      self$Xmodes$XI, self$Omodes$OI, VPhi, s,
                                      self$sigma_mu, self$Z, self$beta)
  self$vec <- posterior$vec
  self$mat <- posterior$mat
  self$mu <- posterior$mu
  self$cov <- posterior$cov
  return(self)
}



V_update <- function(self) {
  # kronecker product of Phi & mu
  A0 <- matrix(0, self$num_subjects * self$num_times, self$K)
  A0 <- kronecker_Phimu(self$K, self$num_subjects, self$num_times,
                        self$Phi, self$mu, A0)

  # hj & Mj
  h <- matrix(0, self$K, self$num_features)
  M <- array(0, dim = c(self$K, self$K, self$num_features))
  Mj <- Mj_creator(self$num_subjects, self$num_times, self$num_features,
                   self$K, self$sigma_noise, self$Xmodes$XJ, self$Omodes$OJ,
                   A0, self$Phi, self$cov, self$OBS)
  h <- Mj$h
  M <- Mj$M

  # update V
  Vtilde <- self$V
  self$V <- V_solver(self$K, self$num_features, self$V, h, M)
  self$Verror <- sum((Vtilde - self$V)^2)
  return(self)
}



Phi_update <- function(self) {
  # kronecker product of V & mu
  A0 <- matrix(0, self$num_subjects * self$num_features, self$K)
  A0 <- kronecker_Vmu(self$K, self$num_subjects, self$num_features,
                      self$V, self$mu, A0)

  s2 <- array(0, self$num_features * self$num_subjects)
  s2 <- s2_transform(self$num_subjects, self$num_features, self$sigma_noise, s2)

  # ht & Mt
  h <- matrix(0, self$K, self$num_times)
  M <- array(0, dim = c(self$K, self$K, self$num_times))
  Mt <- Mt_creator(self$num_subjects, self$num_times, self$num_features,
                   self$K, s2, self$sigma_noise, self$Xmodes$XT, self$Omodes$OT,
                   A0, self$V, self$cov, self$OBS)
  h <- Mt$h
  M <- Mt$M

  # update Phi
  Phitilde <- self$Phi
  Phi_solver_list <- Phi_solver(self$K, self$num_times, self$nlambda1,
                                self$lambda1_dfmin, self$lambda1_dfmax,
                                self$ridge_traj, self$lam1_criterion,
                                self$h, self$update_smooth_penalty,
                                self$lambda1, self$Phi, h, self$Omega, M)
  self$Phi <- Phi_solver_list$Phi
  self$lambda1 <- Phi_solver_list$lambda1
  self$Phierror <- sum((Phitilde - self$Phi)^2) / (self$h)^2
  return(self)
}



beta_fit <- function(Z, mat, vec, beta, intercepts, lambda2,
                     fit_intercept = TRUE, update_lasso_penalty = TRUE,
                     nlambda2 = 100, beta_folds = 5, foldid = NULL,
                     max_iter = 1, tol = 0.01) {
  Z0 <- Z
  beta0 <- beta
  beta1 <- beta
  intercepts0 <- intercepts
  intercepts1 <- intercepts
  K <- dim(vec)[2]
  error <- 2 * tol
  iter <- 0
  while (iter < max_iter && error > tol) {
    if (fit_intercept) {
      a <- array(1, dim(Z)[1])
      Z0 <- cbind(Z, a)
    }
    for (k in 1:K) {
      b <- beta0[, -k]
      tmp0 <- Z %*% b
      if (fit_intercept) {
        tmp1 <- tmp0 + matrix(rep(intercepts0[-k], each = nrow(tmp0)),
                              nrow = dim(tmp0)[1], ncol = dim(tmp0)[2])
      }
      mutrans <- array(0, dim(Z)[1])
      Ztrans <- matrix(0, nrow = dim(Z0)[1], ncol = dim(Z0)[2])
      trans_list <- muZ_transform(dim(Z)[1], k, K, mutrans, Ztrans,
                                  vec, Z0, tmp1, mat)
      mutrans <- trans_list$mutrans
      Ztrans <- trans_list$Ztrans
      penalty_factor <- array(1, dim(Ztrans)[2])
      if (fit_intercept) {
        penalty_factor[dim(Ztrans)[2]] <- 0
      }
      if (update_lasso_penalty) {
        cv <- my_cv_glmnet(Ztrans, mutrans, nlambda2, penalty_factor,
                           beta_folds, foldid)
        lambda2[k] <- cv$lambda.min
        coefficients <- cv$coefficients
      } else {
        coefficients <- my_glmnet(Ztrans, mutrans, nlambda2, lambda2[k],
                                  penalty_factor)
      }
      coefficients <- array(unlist(coefficients))
      p <- length(coefficients)
      if (fit_intercept) {
        intercepts1[k] <- coefficients[p]
      }
      beta1[, k] <- coefficients[1:(p - 1)]
    }
    iter <- iter + 1
    error <- sum((beta1 - beta0)^2)
    beta0 <- beta1
    intercepts0 <- intercepts1
  }
  beta_error <- sum((beta - beta0)^2)
  beta <- beta1
  intercepts <- intercepts1
  beta_list <- list("beta" = beta, "intercepts" = intercepts,
                    "lambda2" = lambda2, "beta_error" = beta_error)
  return(beta_list)
}



beta_update <- function(self) {
  if (self$update_lasso_penalty) {
    beta_list <- beta_fit(self$Z, self$mat, self$vec, self$beta,
                          self$intercepts, self$lambda2, self$fit_intercept,
                          self$update_lasso_penalty, self$nlambda2,
                          self$beta_folds, self$foldid,
                          max_iter = 1, tol = 0.01)
    self$beta <- beta_list$beta
    self$intercepts <- beta_list$intercepts
    self$lambda2 <- beta_list$lambda2
  } else {
    beta_list <- beta_fit(self$Z, self$mat, self$vec, self$beta,
                          self$intercepts, self$lambda2, self$fit_intercept,
                          self$update_lasso_penalty, self$nlambda2,
                          self$beta_folds, self$foldid,
                          max_iter = 1, tol = 0.01)
    self$beta <- beta_list$beta
    self$intercepts <- beta_list$intercepts
    self$lambda2 <- beta_list$lambda2
  }
  return(self)
}



sigma_noise_update <- function(self) {
  # kronecker product of Phi & mu
  A0 <- matrix(0, self$num_subjects * self$num_times, self$K)
  A0 <- kronecker_muPhi(self$K, self$num_subjects, self$num_times,
                        self$Phi, self$mu, A0)

  idx <- which(self$Omodes$OJ[1, ] == 1) - 1
  idxlength <- length(idx)
  A1 <- matrix(0, self$K, self$K)
  O1 <- self$OBS[, , 1]
  A1 <- sigma_helper(self$num_subjects, self$num_times, self$K, idxlength, idx,
                     A0, A1, self$cov, O1, self$Phi)

  # update sigma_noise (variances for Xj)
  sigma <- array(0, self$num_features)
  sigma <- sigma_solver(self$num_subjects, self$num_times, self$num_features,
                        idxlength, self$Xmodes$XJ, A0, A1, self$V, idx)
  sigma <- array(sigma)

  if (self$homo_noise) {
    self$sigma_noise <- sigma
  }
  return(self)
}



sigma_mu_update <- function(self) {
  fit <- matrix(0, self$num_subjects, self$K)
  fit <- sigmaF_helper(self$num_subjects, self$K, self$fit_intercept,
                       self$intercepts, self$Z, self$beta)

  # update sigma_mu (variances for uk)
  sigmaF <- array(0, self$K)
  sigmaF <- sigmaF_solver(self$num_subjects, self$K, fit, self$mu, self$cov)
  sigmaF <- array(sigmaF)
  self$sigma_mu <- sigmaF / self$num_subjects
  return(self)
}



prepare <- function(X, OBS, T1, mean_removal = TRUE, nlam = NULL, lams = NULL,
                    kappa = 1e-2) {
  # prepare input data for initialization
  h <- length(T1)
  basis <- diag(h)
  left <- matrix(0, nrow = h - 1, ncol = h)
  left <- left_creator(h, T1, left)
  Omega <- t(left) %*% left
  Omega <- Omega + kappa * basis
  R <- X
  B0 <- NULL
  cv <- NULL
  dfs <- NULL
  if (mean_removal) {
    if (is.null(nlam) && is.null(lams)) {
      warning("No penalty supplied!")
    }
    O1 <- OBS[, , 1]
    ts_list <- ts_creator(dim(X)[1], dim(X)[2], T1, O1)
    ts <- array(ts_list$ts)
    s0 <- array(ts_list$s0)
    Psi <- basis[ts, ]
    A <- t(Psi) %*% Psi
    if (is.null(lams)) {
      dfs <- 1 + (1 - (array(1:nlam) - 1) / (nlam - 1)) *
        (length(ts) / log(length(ts)) - 1)
      lams <- array(0, nlam)
      lams <- penalty_search(dim(X)[2], nlam, 1000, lams, A, Omega, dfs)
      lams <- array(lams)
    } else {
      nlam <- length(lams)
    }
    cv_mean <- subject_cv_mean(X, Psi, O1, s0, Omega, lams)
    cv <- cv_mean$cv
    B0 <- cv_mean$B0
    R <- R_creator(dim(X)[1], dim(X)[2], dim(X)[3], basis, B0, R)
  }
  data_prepare <- list("Omega" = Omega, "h" = h, "R" = R, "B0" = B0,
                       "cv" = cv, "lams" = lams, "dfs" = dfs)
  return(data_prepare)
}



initialization <- function(self, X, OBS, T1, V_norm, Phi_norm, eps) {
  self$X <- X
  self$OBS <- OBS
  self$T1 <- T1
  Xmodes <- tnsr_unfold(X)
  Omodes <- tnsr_unfold(OBS)
  XI <- Xmodes$tnsr_mode1
  XT <- Xmodes$tnsr_mode2
  XJ <- Xmodes$tnsr_mode3
  OI <- Omodes$tnsr_mode1
  OT <- Omodes$tnsr_mode2
  OJ <- Omodes$tnsr_mode3
  self$Xmodes <- list("XI" = XI, "XT" = XT, "XJ" = XJ)
  self$Omodes <- list("OI" = OI, "OT" = OT, "OJ" = OJ)
  self$sigma_mu <- array(1, self$K)
  self$sigma_noise <- array(1, dim(self$X)[3])
  self$V_norm <- V_norm
  self$Phi_norm <- Phi_norm
  self$eps <- eps
  return(self)
}



reordering <- function(self) {
  mags <- sqrt(colMeans(sweep(self$mu, 2, colMeans(self$mu))^2))
  idx <- sort(-mags, index.return = TRUE)$ix
  self$mu <- self$mu[, idx]
  self$cov <- self$cov[, , idx][, idx, ]
  self$V <- self$V[, idx]
  self$Phi <- self$Phi[, idx]
  self$sigma_mu <- self$sigma_mu[idx]
  self$vec <- self$vec[, idx]
  self$mat <- self$mat[, , idx][, idx, ]
  self$lambda1 <- self$lambda1[idx]
  self$lambda2 <- self$lambda2[idx]
  if (!is.null(self$beta)) {
    self$beta <- self$beta[, idx]
    self$intercepts <- self$intercepts[idx]
  }
  return(self)
}



cross_validation_train <- function(self, train_ids, test_ids, max_iter = 10,
                                   min_iter = 1, tol = 1e-3, trace = TRUE) {
  s <- array(0, self$num_times * self$num_features)
  s <- sigma_transform(self$num_times, self$num_features, self$sigma_noise, s)
  self$train_ids <- train_ids
  self$test_ids <- test_ids
  self$nfolds <- length(train_ids)
  self$cross_V <- array(0, dim = c(self$num_features, self$K, self$nfolds))
  self$cross_Phi <- array(0, dim = c(self$num_times, self$K, self$nfolds))
  self$cross_sigma_mu <- matrix(0, nrow = self$K, ncol = self$nfolds)
  self$cross_mu <- self$mu
  self$cross_cov <- self$cov
  self$cross_vec <- self$vec
  self$cross_mat <- self$mat
  if (!is.null(self$Z)) {
    self$cross_beta <- array(0, dim = c(dim(self$Z)[2], self$K, self$nfolds))
    self$cross_intercept <- matrix(0, nrow = self$K, ncol = self$nfolds)
  } else {
    self$cross_beta <- NULL
    self$cross_intercept <- NULL
  }
  for (fold in seq(1, self$nfolds)) {
    message(paste0("Fold ", fold, ": "))
    # initialize parameters from full fit
    cv <- self
    cv$X <- cv$X[self$train_ids[[fold]], , ]
    cv$OBS <- cv$OBS[self$train_ids[[fold]], , ]
    cv$mu <- cv$mu[self$train_ids[[fold]], ]
    cv$cov <- cv$cov[self$train_ids[[fold]], , ]
    if (!is.null(self$Z)) {
      cv$Z <- cv$Z[self$train_ids[[fold]], ]
    }
    cv$num_subjects <- length(self$train_ids[[fold]])
    Xmodes <- tnsr_unfold(cv$X)
    Omodes <- tnsr_unfold(cv$OBS)
    cv$Xmodes$XI <- Xmodes$tnsr_mode1
    cv$Xmodes$XT <- Xmodes$tnsr_mode2
    cv$Xmodes$XJ <- Xmodes$tnsr_mode3
    cv$Omodes$OI <- Omodes$tnsr_mode1
    cv$Omodes$OT <- Omodes$tnsr_mode2
    cv$Omodes$OJ <- Omodes$tnsr_mode3
    cv$vec <- matrix(0, nrow = cv$num_subjects, cv$K)
    cv$mat <- array(0, dim = c(cv$num_subjects, cv$K, cv$K))
    cv <- check_init(cv)
    cv <- train_spaco(cv, max_iter = max_iter, min_iter = min_iter, tol = tol,
                      update_sigma_mu = TRUE, update_sigma_noise = FALSE,
                      update_lasso_penalty = FALSE,
                      update_smooth_penalty = TRUE,
                      trace = trace, reorder = FALSE)
    self$cross_V[, , fold] <- cv$V
    self$cross_Phi[, , fold] <- cv$Phi
    self$cross_sigma_mu[, fold] <- cv$sigma_mu
    if (!is.null(self$Z)) {
      self$cross_beta[, , fold] <- cv$beta
      if (!is.null(self$fit_intercept)) {
        self$cross_intercept[, fold] <- cv$intercepts
      }
    }
    x_test <- self$Xmodes$XI[self$test_ids[[fold]], ]
    o_test <- self$Omodes$OI[self$test_ids[[fold]], ]
    if (!is.null(self$Z)) {
      z_test <- self$Z[self$test_ids[[fold]], ]
    } else {
      z_test <- NULL
    }

    # find the most correlated VPhi dimension and flip the sign
    VPhi <- matrix(0, nrow = self$num_times * self$num_features, ncol = self$K)
    VPhi <- posterior_kronecker_VPhi(self$K, cv$V, cv$Phi, VPhi)
    posterior <- posterior_distribution(dim(x_test)[1], self$num_times,
                                        self$num_features, self$K,
                                        cv$fit_intercept, cv$intercepts,
                                        x_test, o_test, VPhi, s, cv$sigma_mu,
                                        z_test, cv$beta)
    self$cross_mu[self$test_ids[[fold]], ] <- posterior$mu
    self$cross_cov[self$test_ids[[fold]], , ] <- posterior$cov
    self$cross_vec[self$test_ids[[fold]], ] <- posterior$vec
    self$cross_mat[self$test_ids[[fold]], , ] <- posterior$mat
    if (!is.null(self$Z)) {
      self$cross_beta[, , fold] <- cv$beta
      self$cross_intercept[, fold] <- cv$intercepts
    }
  }
  return(self)
}



log_likelihood <- function(self) {
  log_lik <- array(0, self$num_subjects)
  s <- array(0, self$num_times * self$num_features)
  s <- sigma_transform(self$num_times, self$num_features, self$sigma_noise, s)
  for (k in seq(1, length(self$test_ids))) {
    Phi <- self$cross_Phi[, , k]
    V <- self$cross_V[, , k]
    sigma_mu <- self$cross_sigma_mu[, k]
    cov <- self$cross_cov
    beta <- self$cross_beta[, , k]
    mu <- self$Z %*% beta
    VPhi <- matrix(0, nrow = self$num_times * self$num_features, ncol = self$K)
    VPhi <- posterior_kronecker_VPhi(self$K, V, Phi, VPhi)
    for (i in self$test_ids[[k]]) {
      x1 <- self$Xmodes$XI[i, ]
      o1 <- self$Omodes$OI[i, ]
      o1 <- which(o1 == 1)
      VPhi0 <- VPhi[o1, ]
      x0 <- x1[o1]
      s0 <- s[o1]
      xhat <- array(0, length(x0))
      det_sigma_mu <- 1
      for (j in 1:self$K) {
        det_sigma_mu <- det_sigma_mu * sigma_mu[j]
        xhat <- xhat + VPhi0[, j] * mu[i, j]
      }
      fi <- x0 - xhat
      log_lik[i] <- sum(fi^2 / s0)
      tmp <- (fi / s0) %*% VPhi0
      log_lik[i] <- log_lik[i] - sum(tmp %*% cov[i, , ] * tmp) +
        sum(log(sigma_mu)) + sum(log(s0)) - log(det(cov[i, , ]))
    }
  }
  self$cross_likloss <- log_lik
  return(self)
}
