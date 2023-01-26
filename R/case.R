# Case study: impact

#' @title Data Wrangling for IMPACT Dataset
#' @description This function transforms and maps raw data form of IMPACT into
#' the desired format. Features with more than 20% missingness are excluded.
#' Several covariates are selected including `COVIDRISK_1`, `COVIDRISK_2`,
#' `COVIDRISK_3`, `COVIDRISK_4`, `COVIDRISK_5`, `Age`, `sex`, and `BMI`.
#' @param impact_missing IMPACT dataset with missing values.
#' @param impact_imputed IMPACT dataset after imputing missing values.
#' @return A list `impact` containing the following elements:
#' \item{X}{A 98 (number of subjects) × 35 (number of times) × 111 (number
#' of features) array containing time-series data.}
#' \item{OBS}{A 98 (number of subjects) × 35 (number of times) × 111 (number
#' of features) array indicating whether an observation is available in array
#' `X`.}
#' \item{T1}{A length 35 (number of times) vector containing measured time.}
#' \item{Z}{A 98 (number of subjects) × 8 (number of covariates) matrix
#' containing auxiliary covariates.}
#' \item{Y}{A 98 (number of subjects) × 4 (number of responses) matrix
#' containing responses.}
#' \item{columns_feature}{Name of all features.}
#' @importFrom stats qnorm
#' @examples
#' # IMPACT
#' data("impact_imputed")
#' data("impact_missing")
#' impact <- impact_data_wrangling(impact_missing, impact_imputed)
#' spaco_object <- train_prepare(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
#'                               Z = impact$Z, K = 4, mean_removal = FALSE)
#' spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
#'                             tol = 1e-4, trace = TRUE)
#' @export
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
  columns_feature <- colnames(imputed_pt[, 25:159])[filtered_feature_idx]

  multi_return <- function() {
    return_list <- list("X" = X, "OBS" = OBS, "T1" = T1, "Z" = Z,
                        "columns_feature" = columns_feature)
    return(return_list)
  }
  impact <- multi_return()
  return(impact)
}



#' @title Data Wrangling for IMMUNE Dataset
#' @description This function transforms and maps raw data form of IMMUNE into
#' the desired format. Features with missing values are excluded. Several
#' covariates are selected including `age`, `sex`, `obesity`, `hospitalized`,
#' `severity_group`, `intubated`, `alive`, `tocilizumab`, `heme`, and `bmt`.
#' @param immune_original Original IMMUNE dataset.
#' @return A list `immune` containing the following elements:
#' \item{X}{A 36 (number of subjects) × 40 (number of times) × 133 (number
#' of features) array containing time-series data.}
#' \item{OBS}{A 36 (number of subjects) × 40 (number of times) × 133 (number
#' of features) array indicating whether an observation is available in array
#' `X`.}
#' \item{T1}{A length 40 (number of times) vector containing measured time.}
#' \item{Z}{A 36 (number of subjects) × 10 (number of covariates) matrix
#' containing auxiliary covariates.}
#' \item{columns_feature}{Name of all features.}
#' @importFrom stats qnorm
#' @examples
#' # IMMUNE
#' data("immune_original")
#' immune <- immune_data_wrangling(immune_original)
#' spaco_object <- train_prepare(X = immune$X, OBS = immune$OBS, T1 = immune$T1,
#'                               Z = immune$Z, K = 5, mean_removal = FALSE)
#' spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
#'                             tol = 1e-4, trace = TRUE)
#' @export
immune_data_wrangling <- function(immune_original) {
  columns <- colnames(immune_original)
  columns_covariate <- colnames(immune_original)[1:26]
  columns_feature <- colnames(immune_original)[27:170]
  covariate <- immune_original[, columns_covariate]
  feature <- immune_original[, columns_feature]

  for (i in seq(1, dim(feature)[1])) {
    for (j in seq(1, dim(feature)[2])) {
      if (grepl("%", feature[i, j], fixed = TRUE)) {
        feature[i, j] <- as.numeric(sub("%", "", feature[i, j])) / 100
      }
    }
  }

  chars <- sapply(feature, is.character)
  feature[, chars] <- as.data.frame(apply(feature[, chars], 2, as.numeric))
  covariate <- covariate[-106, ]
  feature <- feature[-106, ]
  # remove columns with na
  feature <- t(na.omit(t(feature)))
  columns_feature <- colnames(feature)
  severity <- covariate$severity_group
  mild_idx <- which(severity == "mild")
  severe_idx <- which(severity == "severe")
  idx <- c(mild_idx, severe_idx)
  feature0 <- feature[idx, ]
  covariate0 <- covariate[idx, ]
  n <- dim(feature0)[1]
  a <- (array(0:(n - 1)) + 0.5) / n
  normalization <- qnorm(a)
  for (i in seq(1, dim(feature0)[2])) {
    feature0[, i] <- normalization[numpy$argsort(numpy$argsort(feature0[, i]))
                                   + 1]
  }

  T0 <- covariate0$time.days.from.symptoms.start
  T0 <- as.double(T0)
  patient_code <- covariate0$patient_code
  patient_code_uni <- unique(patient_code)
  T0_uni <- unique(T0)
  X <- array(0, dim = c(length(patient_code_uni), length(T0_uni),
                        dim(feature0)[2])) # 36,40,133
  X[X == 0] <- NA

  patient_code_uni <- sort.default(patient_code_uni)
  T0_uni <- sort.default(T0_uni)
  Zname <- c("age", "sex", "race", "obesity", "bmi", "hospitalized",
             "severity_group", "intubated", "alive", "tocilizumab",
             "heme", "bmt")
  Z0 <- matrix(0, nrow = length(patient_code_uni), ncol = length(Zname))
  for (i in seq(1, length(patient_code_uni))) {
    idx <- which(covariate0$patient_code == patient_code_uni[i])
    f <- feature0[idx, ]
    c <- covariate0[idx, ]
    for (j in seq(1, length(Zname))) {
      z <- which(columns_covariate == Zname[j])
      Z0[i, j] <- c[1, z]
    }
    for (t in seq(1, length(idx))) {
      tmp <- as.double(c[t, 15])
      t_idx <- which(T0_uni == tmp)
      if (is.null(dim(f))) {
        X[i, t_idx, ] <- as.double(f)
      } else {
        X[i, t_idx, ] <- as.double(f[t, ])
      }
    }
  }

  OBS <- array(0, dim = c(dim(X)[1], dim(X)[2], dim(X)[3]))
  OBS[!is.na(X)] <- 1
  T1 <- sqrt(c(0:(length(T0_uni) - 1))) # 40

  Zdf <- as.data.frame(Z0)
  colnames(Zdf) <- Zname
  Zdf$age <- as.double(Zdf$age)
  Zdf$bmi <- as.double(Zdf$bmi)
  Zdf$severity_group[Zdf$severity_group == "severe"] <- 1
  Zdf$severity_group[Zdf$severity_group == "mild"] <- 0
  Zdf$intubated[Zdf$intubated == "intubated"] <- 1
  Zdf$intubated[Zdf$intubated == "not intubated"] <- 0
  Zdf$tocilizumab[Zdf$tocilizumab == "yes"] <- 1
  Zdf$tocilizumab[Zdf$tocilizumab == "no"] <- 0
  Zdf$sex[Zdf$sex == "f"] <- 1
  Zdf$sex[Zdf$sex == "F"] <- 1
  Zdf$sex[Zdf$sex == "m"] <- 0
  Zdf$sex[Zdf$sex == "M"] <- 0
  Zdf$alive[Zdf$alive == "dead"] <- 1
  Zdf$alive[Zdf$alive == "alive"] <- 0
  Zdf$hospitalized[Zdf$hospitalized == "yes"] <- 1
  Zdf$hospitalized[Zdf$hospitalized == "no"] <- 0
  Zdf$heme[Zdf$heme == "yes"] <- 1
  Zdf$heme[Zdf$heme == "no"] <- 0
  Zdf$bmt[Zdf$bmt == "yes"] <- 1
  Zdf$bmt[Zdf$bmt == "no"] <- 0
  Zdf$obesity[Zdf$obesity == "obese"] <- 2
  Zdf$obesity[Zdf$obesity == "overweight"] <- 1
  Zdf$obesity[Zdf$obesity == "overwheight"] <- 1
  Zdf$obesity[Zdf$obesity == "nonobese"] <- 0
  Zdf$obesity[Zdf$obesity == ""] <- 0
  Zname <- c("age", "sex", "obesity", "hospitalized", "severity_group",
             "intubated", "alive", "tocilizumab", "heme", "bmt")
  Z <- Zdf[, Zname] # 36,10
  Z <- as.matrix(Z)
  class(Z) <- "numeric"
  for (i in seq(1, dim(Z)[2])) {
    Zn <- length(Z[, i])
    Z[, i] <- (Z[, i] - mean(Z[, i])) / sqrt((var(Z[, i]) * (Zn - 1) / Zn))
  }

  multi_return <- function() {
    return_list <- list("X" = X, "OBS" = OBS, "T1" = T1, "Z" = Z,
                        "columns_feature" = columns_feature)
    return(return_list)
  }
  immune <- multi_return()
  return(immune)
}



#' @title Prediction for IMPACT Dataset
#' @description This function evaluates the predictive power of SPACO model
#' trained on IMPACT dataset. The observation for each subject is predicted by
#' \deqn{\hat{x}_{itj} =  \sum_{k=1}^K \hat{u}_{ik} \phi_{kt} v_{jk}}
#' @param spaco_object A list containing the results of model training. It is
#' assumed to include elements `num_subjects`, `num_times`, `K`, `X`, `V`,
#' `Phi`, and `mu`.
#' @return A list containing the predictive values.
#' @importFrom pracma Reshape
#' @examples
#' # IMPACT
#' data("impact_imputed")
#' data("impact_missing")
#' impact <- impact_data_wrangling(impact_missing, impact_imputed)
#' spaco_object <- train_prepare(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
#'                               Z = impact$Z, K = 4, mean_removal = FALSE)
#' spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
#'                             tol = 1e-4, trace = TRUE)
#' spaco_object <- feature_predict(spaco_object)
#'
#' # IMMUNE
#' data("immune_original")
#' immune <- immune_data_wrangling(immune_original)
#' spaco_object <- train_prepare(X = immune$X, OBS = immune$OBS, T1 = immune$T1,
#'                               Z = immune$Z, K = 5, mean_removal = FALSE)
#' spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
#'                             tol = 1e-4, trace = TRUE)
#' spaco_object <- feature_predict(spaco_object)
#' @export
feature_predict <- function(spaco_object) {
  muPhi <- matrix(0, nrow = spaco_object$num_times * spaco_object$num_subjects,
                  ncol = spaco_object$K)
  for (k in seq(1, spaco_object$K)) {
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



#' @title Plot for IMPACT Dataset
#' @description This function plots the observed versus estimated values for
#' the selected feature.
#' @param spaco_object A list containing the results of model training. It is
#' assumed to include elements `num_subjects`, `num_times`, `X`, `OBS`, and
#' `Xhat`.
#' @param selected_feature Name of the selected feature.
#' @param columns_feature Name of all features.
#' @return A plot of the observed versus estimated values.
#' @importFrom graphics par lines points
#' @examples
#' # IMPACT
#' data("impact_imputed")
#' data("impact_missing")
#' impact <- impact_data_wrangling(impact_missing, impact_imputed)
#' spaco_object <- train_prepare(X = impact$X, OBS = impact$OBS, T1 = impact$T1,
#'                               Z = impact$Z, K = 4, mean_removal = FALSE)
#' spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
#'                             tol = 1e-4, trace = TRUE)
#' spaco_object <- feature_predict(spaco_object)
#' feature_plot(spaco_object, "TcellsofLivecells", impact$columns_feature)
#' feature_plot(spaco_object, "TotalNeutrophilsofLivecells",
#'             impact$columns_feature)
#' feature_plot(spaco_object, "HLA.DR.ofTotalMono", impact$columns_feature)
#' feature_plot(spaco_object, "IL6", impact$columns_feature)
#'
#' # IMMUNE
#' data("immune_original")
#' immune <- immune_data_wrangling(immune_original)
#' spaco_object <- train_prepare(X = immune$X, OBS = immune$OBS, T1 = immune$T1,
#'                               Z = immune$Z, K = 5, mean_removal = FALSE)
#' spaco_object <- train_spaco(spaco_object, max_iter = 30, min_iter = 1,
#'                             tol = 1e-4, trace = TRUE)
#' spaco_object <- feature_predict(spaco_object)
#' feature_plot(spaco_object, "IgG..LY", immune$columns_feature)
#' feature_plot(spaco_object, "CD19._CD20..LY.1", immune$columns_feature)
#' feature_plot(spaco_object, "CD4._CD8..CD3.", immune$columns_feature)
#' feature_plot(spaco_object, "LY.All_CD45", immune$columns_feature)
#' @export
feature_plot <- function(spaco_object, selected_feature, columns_feature) {
  par(pty = "s")
  j0 <- which(columns_feature == selected_feature)
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
         ylim = c(-3, 3), cex.lab = 1, cex.main = 1.2, main = selected_feature)
  } else {
    plot(e0, x0, pch = 20, col = 4, xlim = c(-1, 1), ylim = c(-2, 2),
              main = selected_feature)
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
