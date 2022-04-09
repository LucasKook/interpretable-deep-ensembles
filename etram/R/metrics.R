#' Calculate accuracy
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_acc(cdf = cdf, y_true = y_true)
#' @export
get_acc <- function(cdf, y_true) {
  cdf <- as.matrix(cdf)
  pdf <- t(apply(cbind(0, cdf), 1, diff))
  predicted <- apply(pdf, 1,  which.max)
  target <- apply(y_true, 1, which.max)
  ret <- mean(predicted == target)
  return(ret)
}

#' Calculate mean accuracy of ensemble members
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_acc(lys_cdf = lys_cdf, y_true = y_true)
#' @export
get_avg_acc <- function(lys_cdf, y_true, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  weighted.mean(unlist(lapply(lys_cdf, get_acc, y_true = y_true)), w = weights)
}

#' Calculate binary accuracy
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.1, 0.8, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_binacc(cdf = cdf, y_true = y_true, cutoff = 2)
#' @export
get_binacc <- function(cdf, y_true, cutoff = 3) {
  cdf <- as.matrix(cdf)
  pdfbin <- cdf[, cutoff]
  predicted <- ifelse(pdfbin >= 0.5, 1, 0)
  target <- apply(y_true[, 1L:cutoff, drop = FALSE], 1, sum)
  ret <- mean(predicted == target)
  return(ret)
}

#' Calculate mean binary accuracy of ensemble members
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_binacc(lys_cdf = lys_cdf, y_true = y_true, cutoff = 2)
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @export
get_avg_binacc <- function(lys_cdf, y_true, cutoff = 3, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  weighted.mean(unlist(lapply(lys_cdf, get_binacc, y_true = y_true, cutoff = cutoff)), w = weights)
}

#' Calculate ranked probability score
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_rps(cdf = cdf, y_true = y_true)
#' @export
get_rps <- function(cdf, y_true) {
  cdf <- as.matrix(cdf)
  K <- ncol(y_true)
  cdf_true <- t(apply(y_true, 1, cumsum))
  briers <- 1/(K - 1) * apply((cdf_true - cdf)^2, 1, sum)
  ret <- mean(briers)
  return(ret)
}

#' Calculate mean ranked probability score of ensemble members
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_rps(lys_cdf = lys_cdf, y_true = y_true)
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @export
get_avg_rps <- function(lys_cdf, y_true, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  weighted.mean(unlist(lapply(lys_cdf, get_rps, y_true = y_true)), w = weights)
}

#' Calculate negative log-likelihood
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_nll(cdf = cdf, y_true = y_true)
#' @export
get_nll <- function(cdf, y_true) {
  cdf <- as.matrix(cdf)
  pdf <- t(apply(cbind(0, cdf), 1, diff))
  li <- apply(y_true * pdf, 1, sum)
  lli <- log(li)
  nll <- -mean(lli)
  return(nll)
}

#' Calculate mean negative log-likelihood of ensemble members
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_nll(lys_cdf = lys_cdf, y_true = y_true)
#' @export
get_avg_nll <- function(lys_cdf, y_true, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  weighted.mean(unlist(lapply(lys_cdf, get_nll, y_true = y_true)), w = weights)
}

#' Calculate binary negative log-likelihood
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_binnll(cdf = cdf, y_true = y_true, cutoff = 2)
#' @export
get_binnll <- function(cdf, y_true, cutoff = 3) {
  cdf <- as.matrix(cdf)
  pdfbin <- cdf[, cutoff]
  ytruebin <- apply(y_true[, 1L:cutoff, drop = FALSE], 1, sum)
  binli <- ytruebin * pdfbin + (1 - ytruebin) * (1 - pdfbin)
  binnll <- -mean(log(binli))
  return(binnll)
}

#' Calculate mean binary negative log-likelihood of ensemble members
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_binnll(lys_cdf = lys_cdf, y_true = y_true, cutoff = 2)
#' @export
get_avg_binnll <- function(lys_cdf, y_true, cutoff = 3, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  weighted.mean(unlist(lapply(lys_cdf, get_binnll, y_true = y_true, cutoff = cutoff)), w = weights)
}

#' Calculate area under the ROC curve
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_auc(cdf = cdf, y_true = y_true, cutoff = 2)
#' @export
get_auc <- function(cdf, y_true, cutoff = 3) {
  cdf <- as.matrix(cdf)
  pdfbin <- cdf[, cutoff]
  ytruebin <- apply(y_true[, 1L:cutoff, drop = FALSE], 1, sum)
  auc1 <- pROC::auc(ytruebin, pdfbin, levels = c(0, 1), direction = "<")
  auc2 <- pROC::auc(ytruebin, pdfbin, levels = c(1, 0), direction = "<")
  ret <- max(auc1, auc2)
  return(ret)
}

#' Calculate mean area under the ROC curve of ensemble members
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_auc(lys_cdf = lys_cdf, y_true = y_true, cutoff = 2)
#' @export
get_avg_auc <- function(lys_cdf, y_true, cutoff = 3, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  weighted.mean(unlist(lapply(lys_cdf, get_auc, y_true = y_true, cutoff = cutoff)), w = weights)
}

#' Calculate quadratic weighted kappa
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.3, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_qwk(cdf = cdf, y_true = y_true, p = 2)
#' @export
get_qwk <- function(cdf, y_true, p = 2) {
  cdf <- as.matrix(cdf)
  K <- ncol(y_true)
  weights <- ontram:::weight_scheme(K, p)
  pdf <- t(apply(cbind(0, cdf), 1, diff))
  yt <- factor(apply(y_true, 1, which.max), levels = 1:K) # otherwise unused levels will be skipped in table()
  pt <- factor(apply(pdf, 1, which.max), levels = 1:K)
  cmat <- table(yt, pt)
  observed_margin <- apply(cmat, 1, sum)
  predicted_margin <- apply(cmat, 2, sum)
  expected_cmat <- (matrix(observed_margin, ncol = 1) %*%
                      predicted_margin) / nrow(y_true)

  # (sum(weights * cmat) - sum(weights * expected_cmat)) /
  #   (1 - sum(weights * expected_cmat)) # different parameterization requires 1 - weights

  1 - sum(weights * cmat) / sum(weights * expected_cmat)
}

#' Calculate mean quadratic weighted kappa of ensemble members
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_qwk(lys_cdf = lys_cdf, y_true = y_true, p = 2)
#' @export
get_avg_qwk <- function(lys_cdf, y_true, p = 2, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  weighted.mean(unlist(lapply(lys_cdf, get_qwk, y_true = y_true, p = p)), w = weights)
}

#' Calculate calibration intercept and slope per class (K-1)
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_cal_perclass(cdf = cdf, y_true = y_true)
#' @export
get_cal_perclass <- function(cdf, y_true) {
  cdf <- as.matrix(cdf)
  K <- ncol(y_true)
  pdf <- t(apply(cbind(0, cdf), 1, diff))
  pdf[pdf == 0] <- 1e-20
  cint <- numeric(K-1)
  cslope <- numeric(K-1)
  for (cl in seq_len(K-1)) {
    yt <- apply(y_true, 1, function(x) sum(x[1:cl]))
    yp <- apply(pdf, 1, function(x) sum(x[1:cl]))
    yp[yp == 1] <- 0.999
    logits <- qlogis(yp)
    df <- as.data.frame(cbind(y = yt, l = logits))
    m_cint <- glm(I(y == 1) ~ offset(l), data = df, family = "binomial")
    cint[cl] <- coef(m_cint)[1]
    m_cslope <- glm(I(y == 1) ~ l, data = df, family = "binomial")
    cslope[cl] <- coef(m_cslope)[2]
  }
  ret <- list(cint = cint,
              cslope = cslope)
  return(ret)
}

#' Calculate mean calibration intercept and slope across classes
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_cal(cdf = cdf, y_true = y_true)
#' @export
get_cal <- function(cdf, y_true) {
  cdf <- as.matrix(cdf)
  cint <- get_cal_perclass(cdf, y_true)$cint
  cslope <- get_cal_perclass(cdf, y_true)$cslope
  ret <- list(cint = mean(cint),
              cslope = mean(cslope))
  return(ret)
}

#' Calculate mean calibration intercept and slope across classes and ensemble members
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_cal(lys_cdf = lys_cdf, y_true = y_true)
#' @export
get_avg_cal <- function(lys_cdf, y_true, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  cint <- sapply(lapply(lys_cdf, get_cal, y_true = y_true), "[[", 1)
  cslope <- sapply(lapply(lys_cdf, get_cal, y_true = y_true), "[[", 2)
  ret <- list(cint = weighted.mean(cint, w = weights),
              cslope = weighted.mean(cslope, w = weights))
  return(ret)
}

#' Calculate binary brier score
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.3, 0.4, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_brier(cdf = cdf, y_true = y_true, cutoff = 2)
#' @export
get_brier <- function(cdf, y_true, cutoff = 3) {
  cdf <- as.matrix(cdf)
  pdfbin <- cdf[, cutoff]
  ytruebin <- apply(y_true[, 1L:cutoff, drop = FALSE], 1, sum)
  ret <- mean((ytruebin - pdfbin)^2)
  return(ret)
}

#' Calculate mean brier score across ensemble members
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @examples
#' cdf1 <- data.frame(matrix(c(0.2, 0.6, 1,
#'                             0.5, 0.7, 1,
#'                             0.3, 0.4, 1),
#'                           nrow = 3, byrow = T))
#' cdf2 <- data.frame(matrix(c(0.1, 0.8, 1,
#'                             0.3, 0.4, 1,
#'                             0.2, 0.6, 1),
#'                           nrow = 3, byrow = T))
#' lys_cdf <- list(cdf1, cdf2)
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               0, 1, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_avg_brier(lys_cdf = lys_cdf, y_true = y_true, cutoff = 2)
#' @export
get_avg_brier <- function(lys_cdf, y_true, cutoff = 3, weights = rep(1, length(lys_cdf))) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  weighted.mean(unlist(lapply(lys_cdf, get_brier, y_true = y_true, cutoff = cutoff)), w = weights)
}

