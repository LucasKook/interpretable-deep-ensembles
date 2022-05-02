
#' Prepare input to fit \code{k_ontram} or \code{k_ontram_ci} model
#' @examples
#' nn <- function(output_shape = NULL, mbl = FALSE, input_shape = dim(im)[2:4], ...) {
#'   m <- keras_model_sequential() %>%
#'     layer_conv_2d(filters = 16,
#'                   input_shape = input_shape, kernel_size = c(3, 3), activation = "relu") %>%
#'     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#'     layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
#'     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#'     layer_flatten() %>%
#'     layer_dense(units = 50, activation = "relu") %>%
#'     layer_dense(output_shape, activation = "linear", use_bias = FALSE) %>%
#'     {if (mbl){
#'       ontram::layer_trafo_intercept()(.)
#'     } else .
#'     }
#'   return(m)
#' }
#'
#' mnist <- dataset_mnist()
#' c(c(x_train, y_train), c(x_test, y_test)) %<-% mnist
#' x_train <- array_reshape(x_train, c(60000, 28, 28, 1))
#' im <- x_train / 255
#' im_train <- im[1:50,,,, drop = FALSE]
#' im_val <- im[51:70,,,, drop = FALSE]
#' im_test <- im[71:90,,,, drop = FALSE]
#'
#' inp <- get_input(mod = "sics", im_train = im_train, im_val = im_val, im_test = im_test)
#' inp_train <- inp$inp_train
#' inp_val <- inp$inp_val
#' inp_test <- inp$inp_test
#' @export
get_input <- function(mod = c("silscs", "sics", "cils", "ci"),
                      x_train = NULL, x_val = NULL, x_test = NULL,
                      im_train, im_val, im_test) {
  mod <- match.arg(mod)
  one_train <- matrix(1, nrow = nrow(im_train))
  one_val <- matrix(1, nrow = nrow(im_val))
  one_test <- matrix(1, nrow = nrow(im_test))

  if (mod == "silscs") {
    inp_train <- list(one_train, x_train, im_train)
    inp_val <- list(one_val, x_val, im_val)
    inp_test <- list(one_test, x_test, im_test)
  } else if (mod == "sics") {
    inp_train <- list(one_train, im_train)
    inp_val <- list(one_val, im_val)
    inp_test <- list(one_test, im_test)
  } else if (mod == "cils") {
    inp_train <- list(im_train, x_train)
    inp_val <- list(im_val, x_val)
    inp_test <- list(im_test, x_test)
  } else if (mod == "ci") {
    inp_train <- list(im_train, one_train)
    inp_val <- list(im_val, one_val)
    inp_test <- list(im_test, one_test)
  }
  ret <- list(inp_train = inp_train,
              inp_val = inp_val,
              inp_test = inp_test)
  return(ret)
}

#' Set up a \code{k_ontram} or \code{k_ontram_ci} model
#' @param nn function to function to build neural network used for modeling complex intercept or complex shift terms.
#' Arguments \code{input_shape}, \code{output_shape} and \code{mbl} are required.
#' @examples
#' nn <- function(output_shape = NULL, mbl = FALSE, input_shape = dim(im)[2:4], ...) {
#'   m <- keras_model_sequential() %>%
#'     layer_conv_2d(filters = 16,
#'                   input_shape = input_shape, kernel_size = c(3, 3), activation = "relu") %>%
#'     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#'     layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
#'     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#'     layer_flatten() %>%
#'     layer_dense(units = 50, activation = "relu") %>%
#'     layer_dense(output_shape, activation = "linear", use_bias = FALSE) %>%
#'     {if (mbl){
#'       ontram::layer_trafo_intercept()(.)
#'     } else .
#'     }
#'   return(m)
#' }
#'
#' get_model(mod = "sics", y_dim = 10L, x_dim = NULL, nn = nn, input_shape = c(28, 28, 1))
#' @export
get_model <- function(mod = c("silscs", "sics", "cils", "ci"),
                      y_dim, x_dim = NULL,
                      nn, input_shape) {
  mod <- match.arg(mod)
  if (mod == "silscs") {
    mod_si <- k_mod_baseline(y_dim)
    mod_ls <- mod_shift(x_dim)
    mod_cs <- nn(output_shape = 1L, mbl = FALSE, input_shape = input_shape)
    m <- k_ontram(mod_si, list(mod_ls, mod_cs))
  } else if (mod == "sics") {
    mod_si <- k_mod_baseline(y_dim)
    mod_cs <- nn(output_shape = 1L, mbl = FALSE, input_shape = input_shape)
    m <- k_ontram(mod_si, mod_cs)
  } else if (mod == "cils") {
    mod_ci <- nn(output_shape = y_dim - 1L, mbl = TRUE, input_shape = input_shape)
    mod_ls <- mod_shift(x_dim)
    m <- k_ontram(mod_ci, mod_ls, complex_intercept = TRUE)
  } else if (mod == "ci") {
    mod_ci <- nn(output_shape = y_dim - 1, mbl = TRUE, input_shape = input_shape)
    m <- k_ontram(mod_ci, complex_intercept = TRUE)
  }
  return(m)
}

#' Initialize weights for simple or complex intercept and linear shift terms
#' @examples
#' set.seed(2022)
#' nn <- function(output_shape = NULL, mbl = FALSE, input_shape = dim(im)[2:4], ...) {
#'   m <- keras_model_sequential() %>%
#'     layer_conv_2d(filters = 16,
#'                   input_shape = input_shape, kernel_size = c(3, 3), activation = "relu") %>%
#'     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#'     layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
#'     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#'     layer_flatten() %>%
#'     layer_dense(units = 50, activation = "relu") %>%
#'     layer_dense(output_shape, activation = "linear", use_bias = FALSE) %>%
#'     {if (mbl){
#'       ontram::layer_trafo_intercept()(.)
#'     } else .
#'     }
#'   return(m)
#' }
#'
#' mnist <- dataset_mnist()
#' c(c(x_train, y_train), c(x_test, y_test)) %<-% mnist
#'
#' m <- get_model(mod = "silscs", y_dim = 10L, x_dim = 2L, nn = nn, input_shape = c(28, 28, 1))
#'
#' x <- model.frame(data.frame(x1 = rnorm(100), x2 = rnorm(100)))
#' y <- to_categorical(y_train[1:100])
#'
#' get_weights(m$mod_baseline)
#' get_weights(m$list_of_shift_models[[0]])
#' m <- warm_mod(m, mod = "silscs", x = x, y = y, binary = FALSE)
#' get_weights(m$mod_baseline)
#' get_weights(m$list_of_shift_models[[0]])
#' @export
warm_mod <- function(m, mod = c("silscs", "sics", "cils", "ci"), x = NULL, y, binary = FALSE) {
  mod <- match.arg(mod)
  K <- ncol(y)
  y <- apply(y, 1, which.max)
  df <- as.data.frame(cbind(y, x))
  if (!binary) {
    df$y <- factor(df$y, ordered = TRUE)
    if (mod %in% c("silscs", "cils")) {
      mp <- Polr(y ~ ., data = df)
      betas <- coef(mp)
    } else if (mod %in% c("sics", "ci")) {
      mp <- Polr(y ~ 1, data = df)
    }
    thetas <- coef(mp, with_baseline = TRUE)[1:(K - 1)]
  } else if (binary) {
    df$y <- factor(df$y, levels = c(1, 2))
    if (mod %in% c("silscs", "cils")) {
      mg <- glm(y ~ ., family = "binomial", data = df)
      betas <- coef(mg)[2:length(coef(mg))]
    } else if (mod %in% c("sics", "ci")) {
      mg <- glm(y ~ 1, family = "binomial", data = df)
    }
    thetas <- -coef(mg)[1L]
  }

  if (mod %in% c("silscs", "cils")) {
    mw <- warmstart(m, thetas = thetas, betas = betas, which = "all") # also ci from cils
  } else if (mod %in% c("sics", "ci")) {
    mw <- warmstart(m, thetas = thetas, which = "baseline only") # add ci
  }
  return(mw)
}

#' Ordered indices of CDFs based on NLL or RPS
#' @examples
#' cdf1 <- data.frame(matrix(c(0.1, 0.2, 1,
#'                             0.3, 0.7, 1,
#'                             0.5, 0.7, 1),
#'                           nrow = 3, byrow = TRUE))
#' cdf2 <- data.frame(matrix(c(0.2, 0.3, 1,
#'                             0.5, 0.8, 1,
#'                             0.2, 0.9, 1),
#'                           nrow = 3, byrow = TRUE))
#' cdf3 <- data.frame(matrix(c(0.7, 0.8, 1,
#'                             0.6, 0.8, 1,
#'                             0.2, 0.8, 1),
#'                           nrow = 3, byrow = TRUE))
#' lys_cdf <- list(cdf1, cdf2, cdf3)
#' y_true <- data.frame(matrix(c(0, 0, 1,
#'                               1, 0, 0,
#'                               0, 1, 0),
#'                            nrow = 3, byrow = TRUE))
#'
#' get_order(lys_cdf_val = lys_cdf, y_true = y_true, order_metric = "rps")
#' get_order(lys_cdf_val = lys_cdf, y_true = y_true, order_metric = "nll")
#' @export
get_order <- function(lys_cdf_val, y_true_val, order_metric = c("nll", "rps")) {
  order_metric <- match.arg(order_metric)
  lys_cdf_val <- lapply(lys_cdf_val, as.matrix)
  if (order_metric == "rps") {
    indiv <- lapply(lys_cdf_val, get_rps, y_true = y_true_val)
  } else if (order_metric == "nll") {
    indiv <- lapply(lys_cdf_val, get_nll, y_true = y_true_val)
  }
  indiv <- do.call("rbind", indiv)
  switch(
    order_metric,
    "rps" = order(indiv, decreasing = FALSE),
    "nll" = order(indiv, decreasing = FALSE),
  )
}

#' Ensemble
#' @examples
#' cdf1 <- data.frame(matrix(c(0.1, 0.2, 1,
#'                             0.3, 0.7, 1,
#'                             0.5, 0.7, 1),
#'                           nrow = 3, byrow = TRUE))
#' cdf2 <- data.frame(matrix(c(0.2, 0.3, 1,
#'                             0.5, 0.8, 1,
#'                             0.2, 0.9, 1),
#'                           nrow = 3, byrow = TRUE))
#' cdf3 <- data.frame(matrix(c(0.7, 0.8, 1,
#'                             0.6, 0.8, 1,
#'                             0.2, 0.8, 1),
#'                           nrow = 3, byrow = TRUE))
#' lys_cdf <- list(cdf1, cdf2, cdf3)
#' y_true <- data.frame(matrix(c(0, 0, 1,
#'                               1, 0, 0,
#'                               0, 1, 0),
#'                            nrow = 3, byrow = TRUE))
#'
#' # equal weights
#' get_ensemble(lys_cdf = lys_cdf, type = "linear")
#'
#' # unequal weights
#' w <- get_w(lys_cdf_val = lys_cdf, y_true_val = y_true, type = "linear", optim_metric = "rps")
#' get_ensemble(lys_cdf = lys_cdf, type = "linear", weights = w)
#' @export
get_ensemble <- function(lys_cdf,
                         type = c("linear", "log-linear", "trafo"),
                         weights = rep(1, length(lys_cdf))) {
  type <- match.arg(type)
  lys_cdf <- lapply(lys_cdf, as.matrix)
  w <- weights
  switch(
    type,
    "linear" = apply(simplify2array(lys_cdf), 1:2,
                     function(x) weighted.mean(x = x, w = w)),
    "log-linear" = exp(apply(simplify2array(lapply(lys_cdf, log)), 1:2,
                             function(x) weighted.mean(x = x, w = w))),
    "trafo" = plogis(apply(simplify2array(lapply(lys_cdf, qlogis)), 1:2,
                           function(x) weighted.mean(x = x, w = w)))
  )
}

#' Get weights per ensemble member
#' @examples
#' cdf1 <- data.frame(matrix(c(0.1, 0.2, 1,
#'                             0.3, 0.7, 1,
#'                             0.5, 0.7, 1),
#'                           nrow = 3, byrow = TRUE))
#' cdf2 <- data.frame(matrix(c(0.2, 0.3, 1,
#'                             0.5, 0.8, 1,
#'                             0.2, 0.9, 1),
#'                           nrow = 3, byrow = TRUE))
#' cdf3 <- data.frame(matrix(c(0.7, 0.8, 1,
#'                             0.6, 0.8, 1,
#'                             0.2, 0.8, 1),
#'                           nrow = 3, byrow = TRUE))
#' lys_cdf <- list(cdf1, cdf2, cdf3)
#' y_true <- data.frame(matrix(c(0, 0, 1,
#'                               1, 0, 0,
#'                               0, 1, 0),
#'                            nrow = 3, byrow = TRUE))
#'
#' get_w(lys_cdf_val = lys_cdf, y_true_val = y_true, type = "trafo", optim_metric = "nll")
#' @export
get_w <- function(lys_cdf_val,
                  y_true_val,
                  type = c("linear", "log-linear", "trafo"),
                  optim_metric = c("nll", "rps")) {
  lys_cdf_val <- lapply(lys_cdf_val, as.matrix)
  nens <- length(lys_cdf_val)
  start <- 1/nens
  w <- optim(par = rep(start, nens), fn = .opt, method = "L-BFGS-B",
             lower = rep(0, nens), upper = rep(1, nens),
             lys_cdf_val = lys_cdf_val, y_true_val = y_true_val, type = type,
             optim_metric = optim_metric)$par
  ret <- w / sum(w)
  return(ret)
}

.opt <- function(w,
                 lys_cdf_val,
                 y_true_val,
                 type = c("linear", "log-linear", "trafo"),
                 optim_metric = c("nll", "rps")) {
  optim_metric <- match.arg(optim_metric)
  type <- match.arg(type)
  w <- w / sum(w) # ensures that weights sum up to 1
  weighted_ens <- get_ensemble(lys_cdf = lys_cdf_val, type = type, weights = w)
  if (optim_metric == "nll") {
    ret <- get_nll(cdf = weighted_ens, y_true = y_true_val)
    if (!is.finite(ret)) {
      ret <- 10^6 # optim needs finite values
    }
  } else if (optim_metric == "rps") {
    ret <- get_rps(cdf = weighted_ens, y_true = y_true_val)
    if (!is.finite(ret)) {
      ret <- 1 # optim needs finite values
    }
  }
  return(ret)
}

#' Get predicted and observed probabilities
#' @examples
#' cdf <- data.frame(matrix(c(0.2, 0.6, 1,
#'                            0.5, 0.7, 1,
#'                            0.1, 0.8, 1),
#'                          nrow = 3, byrow = T))
#' y_true <- data.frame(matrix(c(0, 1, 0,
#'                               1, 0, 0,
#'                               0, 0, 1),
#'                             nrow = 3, byrow = T))
#' get_bins(cdf = cdf, y_true = y_true, bins = 10)
#' @export
get_bins <- function(cdf, y_true, bins = 10) {
  cdf <- as.matrix(cdf)
  K <- ncol(y_true)
  pdf <- t(apply(cbind(0, cdf), 1, diff))
  ret <- data.frame(n = numeric(),
                    pred = numeric(),
                    obs = numeric(),
                    se = numeric(),
                    uci = numeric(),
                    lci = numeric(),
                    class = factor())
  for (cl in 1:(K-1)) {
    yt <- apply(y_true, 1, function(x) sum(x[1:cl]))
    yp <- apply(pdf, 1, function(x) sum(x[1:cl]))
    df <- data.frame(pred = yp, true = yt)
    tmp <- df %>% mutate(bin = cut(pred, breaks = seq(0, 1, by = 1/bins), labels = FALSE)) %>%
      group_by(bin) %>%
      summarise(n = n(),
                pred = mean(pred),
                obs = mean(true),
                se = sqrt((obs * (1 - obs)) / n),
                uci = obs + qnorm(0.975) * se,
                lci = obs - qnorm(0.975) * se,
                class = cl)
    ret <- rbind(ret, tmp)
  }
  return(ret)
}

#' Data generator
#' @description  Data generator that feeds data real-time to model.
#' @export
generator <- function(mod = c("silscs", "sics", "cils", "ci"),
                      im, x = NULL, y,
                      batch_size = 32,
                      shuffle = TRUE) {
  mod <- match.arg(mod)
  n <- nrow(im)
  # start iterator
  i <- 1

  function() {
    if (shuffle) {
      ridx <- sample(seq_len(n), size = batch_size)
    } else {
      if (i + batch_size >= n) {
        i <<- 1
      }
      ridx <- c(i:min(i + batch_size - 1, n))
      i <<- i + length(ridx)
    }
    if (!is.null(x)) {
      xsub <- x[ridx, , drop = FALSE]
    } else {
      xsub <- NULL
    }
    imsub <- ontram:::.batch_subset(im, ridx, dim = dim(im))

    xret <- .inp(mod = mod, im = imsub, x = xsub)
    yret <- y[ridx, , drop = FALSE]
    ret <- list(xret, yret)
    return(ret)
  }
}

.inp <- function(mod = c("silscs", "sics", "cils", "ci"),
                 im, x = NULL) {
  mod <- match.arg(mod)
  n <- nrow(im)
  one <- matrix(1, nrow = n)

  if (mod == "silscs") {
    ret <- list(one, x, im)
  } else if (mod == "sics") {
    ret <- list(one, im)
  } else if (mod == "cils") {
    ret <- list(im, x)
  } else if (mod == "ci") {
    ret <- list(im, one)
  }
  return(ret)
}

#' Bootstrap confidence intervals across all splits
#' @description  Calculates bootstrap median and confidence intervals across all splits
#' for binary (NLL, Brier score, 1-AUC, 1-accuracy, calibration-in-the-large, calibration slope)
#' or ordinal (NLL, RPS, 1-QWK, 1-accuracy, calibration-in-the-large, calibration slope) metrics.
#' If desired the confidence interval is additionally calculated for the test error relative to a
#' reference model (performance treated as fixed).
#' @param lys_cdf_all list of lists. Each sublist contains the CDFs of all ensemble members
#' or the ensemble CDF per split.
#' @param y_true_all list of all observed responses (one-hot encoded).
#' @param met_ref optional. Test performance of the reference model (e.g. simple linear shift model)
#' as a data frame. Must contain the variables \code{spl} and \code{metric} with the following levels:
#' \code{c('nll', 'brier', 'eauc', 'eacc', 'cint', 'cslope')} or
#' \code{c('nll', 'rps', 'eqwk', 'eacc', 'cint', 'cslope')}.
#' @param R number of bootstrap samples.
#' @param weights list of optimized weights. Each element contains the weights for each
#' ensemble member as numeric vector.
get_bs <- function(lys_cdf_all, y_true_all, met_ref = NULL,
                   R = 1000,
                   weights = rep(list(rep(1, length(lys_cdf_all[[1]]))),
                                 length(lys_cdf_all)), binary = FALSE) {
  K <- ncol(lys_cdf_all[[1]][[1]])
  spl <- length(lys_cdf_all)
  mem <- length(lys_cdf_all[[1]])
  if (!is.null(met_ref)) {
    ref <- TRUE
    met_ref$metric <- paste0("w", met_ref$metric) # same names as boot res
    met_ref$spl <- factor(met_ref$spl) # otw error in left_join
  } else {
    ref <- FALSE
  }
  # revert one-hot encoding
  yt_all <- lapply(y_true_all, function(spl) apply(spl, 1, which.max))
  # convert to df, create vars spl, mem, weights
  for (s in seq_len(spl)) {
    for (m in seq_len(mem)) {
      lys_cdf_all[[s]][[m]] <- data.frame(lys_cdf_all[[s]][[m]])
      lys_cdf_all[[s]][[m]]$yt <- yt_all[[s]]
      lys_cdf_all[[s]][[m]]$mem <- factor(m)
      lys_cdf_all[[s]][[m]]$w <- weights[[s]][m]
    }
  }
  lys_cdf_all <- lapply(lys_cdf_all, function(spl) do.call("rbind", spl))
  for (s in seq_len(spl)) {
    lys_cdf_all[[s]]$spl <- factor(s)
  }
  # df with first k cols cdf, y_true, mem, spl
  dd <- do.call("rbind", lys_cdf_all)
  res <- boot(dd, statistic = .statFun, R = R, K = K, binary = binary, met_ref = met_ref)
  nmet <- ncol(res$t)/spl
  ret <- list()
  for (m in seq_len(nmet)) {
    avg_acr_spl <- rowMeans(res$t[, ((m - 1) * spl + 1) : ((m - 1) * spl + spl)])
    ret <- c(ret, list(data.frame(t(quantile(avg_acr_spl, c(0.025, 0.5, 0.975), na.rm = TRUE)))))
  }
  ret <- do.call("rbind", ret)
  if (binary) {
    if (!ref) {
      ret$metric <- c("nll", "brier", "eauc", "eacc", "cint", "cslope")
    } else {
      ret$metric <- c("nll", "brier", "eauc", "eacc", "cint", "cslope",
                      "dnll", "dbrier", "deauc", "deacc")
    }
  } else {
    if (!ref) {
      ret$metric <- c("nll", "rps", "eqwk", "eacc", "cint", "cslope")
    } else {
      ret$metric <- c("nll", "rps", "eqwk", "eacc", "cint", "cslope",
                      "dnll", "drps", "deqwk", "deacc")
    }
  }
  colnames(ret) <- c("lwr", "med", "upr", "metric")
  return(ret)
}

.onehot <- function(y, K) {
  rnks <- as.numeric(as.factor(rank(y)))
  ret <- matrix(0, nrow = length(rnks), ncol = K)
  for (r in seq_len(nrow(ret))) {
    ret[r, rnks[r]] <- ret[r, rnks[r]] + 1
  }
  return(ret)
}

.get_nll_bs <- function(d, K) {
  # cols 1:K cdf, K + 1: y_true, cdf ref
  cdf <- d[, 1:K]
  y_true <- data.frame(.onehot(d[, K+1], K = K))
  get_nll(cdf = cdf, y_true = y_true)
}

.get_rps_bs <- function(d, K) {
  cdf <- d[, 1:K]
  y_true <- data.frame(.onehot(d[, K+1], K = K))
  get_rps(cdf = cdf, y_true = y_true)
}

.get_eacc_bs <- function(d, K) {
  cdf <- d[, 1:K]
  y_true <- data.frame(.onehot(d[, K+1], K = K))
  1 - get_acc(cdf = cdf, y_true = y_true)
}

.get_eauc_bs <- function(d, K) {
  cdf <- d[, 1:K]
  y_true <- data.frame(.onehot(d[, K+1], K = K))
  1 - get_auc(cdf = cdf, y_true = y_true)
}

.get_eqwk_bs <- function(d, K, p = 2) {
  cdf <- d[, 1:K]
  y_true <- data.frame(.onehot(d[, K+1], K = K))
  1 - get_qwk(cdf = cdf, y_true = y_true, p = 2)
}

.get_cint_bs <- function(d, K) {
  cdf <- d[, 1:K]
  y_true <- data.frame(.onehot(d[, K+1], K = K))
  get_cal(cdf = cdf, y_true = y_true)[[1]]
}

.get_cslope_bs <- function(d, K) {
  cdf <- d[, 1:K]
  y_true <- data.frame(.onehot(d[, K+1], K = K))
  get_cal(cdf = cdf, y_true = y_true)[[2]]
}

.get_brier_bs <- function(d, K) {
  cdf <- d[, 1:K]
  y_true <- data.frame(.onehot(d[, K+1], K = K))
  get_brier(cdf = cdf, y_true = y_true)
}


.statFun <- function(d, i, K, binary, met_ref) {
  if (!is.null(met_ref)) ref <- TRUE else ref <- FALSE
  d %>%
    slice(i) %>%
    group_by(mem, spl) %>%
    {if (binary) {
      {.} %>% dplyr::summarise(nll = .get_nll_bs(cur_data(), K = K),
                               brier = .get_brier_bs(cur_data(), K = K),
                               eauc = .get_eauc_bs(cur_data(), K = K),
                               eacc = .get_eacc_bs(cur_data(), K = K),
                               cint = .get_cint_bs(cur_data(), K = K),
                               cslope = .get_cslope_bs(cur_data(), K = K),
                               w = unique(w))
    } else {
      {.} %>% dplyr::summarise(nll = .get_nll_bs(cur_data(), K = K),
                               rps = .get_rps_bs(cur_data(), K = K),
                               eqwk = .get_eqwk_bs(cur_data(), K = K),
                               eacc = .get_eacc_bs(cur_data(), K = K),
                               cint = .get_cint_bs(cur_data(), K = K),
                               cslope = .get_cslope_bs(cur_data(), K = K),
                               w = unique(w))
    }} %>%
    group_by(spl) %>%
    {if (binary) {
      {.} %>% dplyr::summarise(wnll = weighted.mean(nll, w),
                               wbrier = weighted.mean(brier, w),
                               weauc = weighted.mean(eauc, w),
                               weacc = weighted.mean(eacc, w),
                               wcint = weighted.mean(cint, w),
                               wcslope = weighted.mean(cslope, w))
    } else {
      {.} %>% dplyr::summarise(wnll = weighted.mean(nll, w),
                               wrps = weighted.mean(rps, w),
                               weqwk = weighted.mean(eqwk, w),
                               weacc = weighted.mean(eacc, w),
                               wcint = weighted.mean(cint, w),
                               wcslope = weighted.mean(cslope, w))
    }} %>%
    {if (ref) {
      {.} %>%
        gather("metric", "val", -spl) %>%
        left_join(met_ref, by = c("metric", "spl")) %>%
        mutate(diff = val.x - val.y,
               dmetric = paste0("d", metric)) %>%
        select(metric, val.x, diff) %>%
        gather("key", "val", -metric) %>%
        mutate(metric = case_when(key == "val.x" ~ metric,
                                  key == "diff" ~ paste0("d", metric))) %>%
        select(metric, val) %>%
        pivot_wider(names_from = metric, values_from = val) %>% unnest()
    } else (.)} %>%
    {if (binary) {
      if (!ref) {
        {.} %>% dplyr::select(wnll, wbrier, weauc, weacc, wcint, wcslope) %>%  as.matrix(ncol = 1)
      } else {
        {.} %>% dplyr::select(wnll, wbrier, weauc, weacc, wcint, wcslope,
                              dwnll, dwbrier, dweauc, dweacc) %>%  as.matrix(ncol = 1)
      }
    } else {
      if (!ref) {
        {.} %>% dplyr::select(wnll, wrps, weqwk, weacc, wcint, wcslope) %>% as.matrix(ncol = 1)
      } else {
        {.} %>%  dplyr::select(wnll, wrps, weqwk, weacc, wcint, wcslope,
                               dwnll, dwrps, dweqwk, dweacc) %>% as.matrix(ncol = 1)
      }
    }}
}
