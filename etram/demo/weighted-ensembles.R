
library(etram)
library(tram)

set.seed(74823)

# Functions ---------------------------------------------------------------

dgp <- function(x, b, p_baseline) {
  K <- length(p_baseline) + 1
  v <- log(p_baseline / (1 - p_baseline))
  p <- matrix(0, ncol = K, nrow = nrow(x))
  for (k in seq_len(K)) {
    if (k == 1) {
      p[, k] <- plogis(v[k] - as.matrix(x) %*% b)
    } else if (k == K) {
      p[,k] <- 1 - plogis(v[k-1] - as.matrix(x) %*% b)
    } else {
      p[,k] <- plogis(v[k] - as.matrix(x) %*% b) - plogis(v[k-1] - as.matrix(x) %*% b)
    }
  }
  y <- apply(p, 1, function(p) sample(K, 1, prob = p))
  ret <- data.frame(cbind(y, x))
  ret$y <- ordered(y, levels = 1:K)
  return(ret)
}

# Sim ---------------------------------------------------------------------

K <- 4
n <- 200
b <- log(c(3, 2.5))
x <- data.frame(x1 = rnorm(n, 0, 1), x2 = rnorm(n, 0, 1))
p <- seq(0.2, 0.8, length.out = K - 1)
dat_true <- dgp(x, b, p)
y_true <- data.frame(as.matrix(model.matrix(~ 0 + y, data = dat_true)))

# sample from dgp distribution
d1 <- dgp(x, b, p)
d2 <- dgp(x, b, p)
d3 <- dgp(x, b, p)

# fit models to get coefficients
m1 <- Polr(y ~ x1 + x2, data = d1)
m2 <- Polr(y ~ x1 + x2, data = d2)
m3 <- Polr(y ~ x1 + x2, data = d3)

# predict to get cdf
cdf1 <- t(predict(m1, newdata = d1[, -1L], type = "distribution"))
cdf2 <- t(predict(m1, newdata = d2[, -1L], type = "distribution"))
cdf3 <- t(predict(m1, newdata = d3[, -1L], type = "distribution"))

# mix predictions of one model
pdf3 <- t(apply(cbind(0, cdf3), 1, diff))
pdf3_mix <- pdf3[, c(4, 3, 2, 1)]
cdf3_mix <- t(apply(pdf3_mix, 1, cumsum))

lys_cdf <- list(cdf1, cdf2, cdf3)
lys_cdf_mix <- list(cdf1, cdf2, cdf3_mix)

# Results -----------------------------------------------------------------

# All CDFs were generated from the same data generating distribution
# --> should result in equal weights
get_w(lys_cdf_val = lys_cdf, y_true_val = y_true, "linear", "nll")
get_w(lys_cdf_val = lys_cdf, y_true_val = y_true, "log-linear", "nll")
get_w(lys_cdf_val = lys_cdf, y_true_val = y_true, "trafo", "nll")
get_w(lys_cdf_val = lys_cdf, y_true_val = y_true, "linear", "rps")
get_w(lys_cdf_val = lys_cdf, y_true_val = y_true, "log-linear", "rps")
get_w(lys_cdf_val = lys_cdf, y_true_val = y_true, "trafo", "rps")

# CDF 3 is mixed --> should get lower weight than the others
get_w(lys_cdf_val = lys_cdf_mix, y_true_val = y_true, "linear", "nll")
get_w(lys_cdf_val = lys_cdf_mix, y_true_val = y_true, "log-linear", "nll")
get_w(lys_cdf_val = lys_cdf_mix, y_true_val = y_true, "trafo", "nll")
get_w(lys_cdf_val = lys_cdf_mix, y_true_val = y_true, "linear", "rps")
get_w(lys_cdf_val = lys_cdf_mix, y_true_val = y_true, "log-linear", "rps")
w <- get_w(lys_cdf_val = lys_cdf_mix, y_true_val = y_true, "trafo", "rps")

get_ensemble(lys_cdf = lys_cdf, type = "trafo", weights = w)
