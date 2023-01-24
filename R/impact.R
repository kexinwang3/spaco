# Case study: impact

#' @importFrom stats qnorm
impact_data_wrangling <- function(impact_missing, impact_imputed) {
  missing <- impact_missing
  imputed <- impact_imputed

  # sex
  missing$sex[missing$sex == "M"] <- 1
  missing$sex[missing$sex == "F"] <- 0
  imputed$sex[imputed$sex == "M"] <- 1
  imputed$sex[imputed$sex == "F"] <- 0

  # filter out features with more than 20% missingness
  missing_pt <- missing[substr(missing$ID, 1, 2) == "Pt", ]
  missing_X_pt <- missing_pt[, 25:159]
  filtered_feature_idx <- which(colMeans(is.na(missing_X_pt)) <= 0.2)

  imputed_pt <- imputed[substr(imputed$ID, 1, 2) == "Pt", ]
  Z_pt <- imputed_pt[, 1:24]
  X_pt <- imputed_pt[, 25:159]

  X_pt <- X_pt[, filtered_feature_idx]
  n <- dim(X_pt)[1]
  a <- (array(0:(n - 1)) + 0.5) / n
  normalization <- qnorm(a)
  colnames(X_pt) <- NULL
  X_pt <- matrix(unlist(X_pt), nrow = dim(X_pt)[1])
  for (i in seq(1, dim(X_pt)[2])) {
    X_pt[, i] <- normalization[numpy$argsort(numpy$argsort(X_pt[, i])) + 1]
  }
  ids <- substr(Z_pt$ID, 1, 5)
  ids_unique <- unique(ids)
  DFSO <- Z_pt$DFSO
  T0 <- sort(unique(DFSO))
  I <- length(ids_unique)
  T1 <- length(T0)
  J <- dim(X_pt)[2]


  Zname <- c("COVIDRISK_1", "COVIDRISK_2", "COVIDRISK_3", "COVIDRISK_4",
             "COVIDRISK_5", "Age", "sex", "BMI")
  Yname <- c("ICU", "Clinicalscore", "LengthofStay", "Coagulopathy")
  Z1 <- as.matrix(sapply(Z_pt[, Zname], as.numeric))
  Y1 <- as.matrix(sapply(Z_pt[, Yname], as.numeric))
  X0 <- array(NA, dim = c(I, T1, J)) # 98,35,111
  Y0 <- array(NA, dim = c(I, T1, length(Yname))) # 98,35,4
  Z0 <- array(NA, dim = c(I, T1, length(Zname))) # 98,35,8
  for (i in 1:I) {
    idx_ids <- which(ids == ids_unique[i])
    ti <- DFSO[idx_ids]
    for (t in seq(1, length(ti))) {
      idx_t <- which(T0 == ti[t])
      X0[i, idx_t, ] <- X_pt[idx_ids[t], ]
      Y0[i, idx_t, ] <- Y1[idx_ids[t], ]
      Z0[i, idx_t, ] <- Z1[idx_ids[t], ]
    }
  }

  OBS <- array(0, dim = c(I, T1, J))
  OBS[!is.na(X0)] <- 1
  Z <- matrix(NA, nrow = I, ncol = length(Zname))
  Y <- matrix(NA, nrow = I, ncol = length(Yname))
  for (i in seq(1, I)) {
    for (j in seq(1, length(Yname))) {
      Y[i, j] <- suppressWarnings(max(Y0[i, , j], na.rm = TRUE))
    }
    for (j in seq(1, dim(Z)[2])) {
      idx <- which(!is.na(Z0[i, , j]))
      if (length(idx) == 0) {
        Z[i, j] <- 0
      } else {
        Z[i, j] <- Z0[i, idx[1], j]
      }
    }
  }
  # standardize Z
  for (j in seq(1, length(Zname))) {
    Z[, j] <- (Z[, j] - mean(Z[, j])) /
      sqrt(colMeans(sweep(Z, 2, colMeans(Z))^2))[j]
  }

  Y[is.infinite(Y)] <- NA
  Z[is.na(Z)] <- 0
  X <- X0
  T1 <- sqrt(seq(1, length(T0)) - 1)

  multi_return <- function() {
    return_list <- list("X" = X, "OBS" = OBS, "Z" = Z, "T1" = T1, "Y" = Y,
                        "imputed_pt" = imputed_pt,
                        "filtered_feature_idx" = filtered_feature_idx)
    return(return_list)
  }
  impact <- multi_return()
  return(impact)
}


#' @importFrom pracma Reshape
impact_predict <- function(spaco_object) {
  muPhi <- matrix(0, nrow = spaco_object$num_times * spaco_object$num_subjects,
                  ncol = spaco_object$K)
  for (k in 1:spaco_object$K) {
    muPhi[, k] <- Reshape(kronecker(Reshape(spaco_object$Phi[, k],
                                            spaco_object$num_times, 1),
                                    Reshape(spaco_object$mu[, k],
                                            spaco_object$num_subjects, 1)),
                         1, spaco_object$num_times * spaco_object$num_subjects)
  }
  spaco_object$Xhat <- array(0, dim = c(dim(spaco_object$X)[1],
                                        dim(spaco_object$X)[2],
                                        dim(spaco_object$X)[3]))
  for (j in seq(1, dim(spaco_object$X)[3])) {
    tmp <- muPhi %*% spaco_object$V[j, ]
    tmp.reshape <- array(tmp)
    tmp.reshape <- matrix(tmp.reshape, nrow = spaco_object$num_times,
                          ncol = spaco_object$num_subjects, byrow = TRUE)
    spaco_object$Xhat[, , j] <- t(tmp.reshape)
  }
  return(spaco_object)
}


#' @importFrom graphics par lines points
impact_plot <- function(spaco_object, feature_name,
                        imputed_pt, filtered_feature_idx) {
  par(pty = "s")
  columns_feature <- colnames(imputed_pt[, 25:159])[filtered_feature_idx]
  j0 <- which(columns_feature == feature_name)
  x0 <- Reshape(spaco_object$X[1, , j0], 1, spaco_object$num_times)
  e0 <- Reshape(spaco_object$Xhat[1, , j0], 1, spaco_object$num_times)
  o0 <- spaco_object$OBS[1, , 1]
  idx0 <- which(o0 == 1)
  x0 <- x0[idx0]
  e0 <- e0[idx0]

  if (length(idx0) > 1) {
    s1 <- sort(e0, index.return = TRUE)$ix
    s2 <- sort(s1, index.return = TRUE)$ix
    plot(e0, x0, type = "l", lwd = 2, col = 4,
         xlab = "Estimated", ylab = "Observed", xlim = c(-2, 2),
         ylim = c(-3, 3), cex.lab = 1, cex.main = 1.2, main = feature_name)
  } else {
    plot(e0, x0, pch = 20, col = 4, xlim = c(-1, 1), ylim = c(-2, 2),
              main = feature_name)
  }

  for (i in seq(2, spaco_object$num_subjects)) {
    x1 <- Reshape(spaco_object$X[i, , j0], 1, spaco_object$num_times)
    e1 <- Reshape(spaco_object$Xhat[i, , j0], 1, spaco_object$num_times)
    o1 <- spaco_object$OBS[i, , 1]
    idx <- which(o1 == 1)
    x1 <- x1[idx]
    e1 <- e1[idx]
    if (length(idx) > 1) {
      s1 <- sort(e1, index.return = TRUE)$ix
      s2 <- sort(s1, index.return = TRUE)$ix
      lines(e1[s2], x1[s2], lwd = 2, col = 4)
    } else {
      points(e1, x1, pch = 20, col = 4)
    }
  }
}
