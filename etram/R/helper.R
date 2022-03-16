
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

#' Initialize weights for simple intercept and linear shift term
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
warm_mod <- function(m, mod = c("silscs", "sics", "cils"), x = NULL, y, binary = FALSE) {
  mod <- match.arg(mod)
  K <- ncol(y)
  y <- apply(y, 1, which.max)
  df <- as.data.frame(cbind(y, x))
  if (!binary) {
    df$y <- factor(df$y, ordered = TRUE)
    if (mod %in% c("silscs", "cils")) {
      mp <- Polr(y ~ ., data = df)
      betas <- coef(mp)
    } else if (mod == "sics") {
      mp <- Polr(y ~ 1, data = df)
    }
    thetas <- coef(mp, with_baseline = TRUE)[1:(K - 1)]
  } else if (binary) {
    df$y <- factor(df$y, levels = c(1, 2))
    if (mod %in% c("silscs", "cils")) {
      mg <- glm(y ~ ., family = "binomial", data = df)
      betas <- coef(mg)[2:length(coef(mg))]
    } else if (mod == "sics") {
      mg <- glm(y ~ 1, family = "binomial", data = df)
    }
    thetas <- coef(mg)[1L]
  }

  if (mod == "silscs") {
    mw <- warmstart(m, thetas = thetas, betas = betas, which = "all")
  } else if (mod == "sics") {
    mw <- warmstart(m, thetas = thetas, which = "baseline only")
  } else if (mod == "cils") {
    mw <- warmstart(m, betas = betas, which = "shift only")
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
#'
#' get_ensemble(lys_cdf = lys_cdf, type = "linear")
#' get_ensemble(lys_cdf = lys_cdf, type = "log-linear")
#' get_ensemble(lys_cdf = lys_cdf, type = "trafo")
#' @export
get_ensemble <- function(lys_cdf, type = c("linear", "log-linear", "trafo")) {
  type <- match.arg(type)
  lys_cdf <- lapply(lys_cdf, as.matrix)
  switch(
    type,
    "linear" = apply(simplify2array(lys_cdf), 1:2, mean),
    "log-linear" = exp(apply(simplify2array(lapply(lys_cdf, log)), 1:2, mean)),
    "trafo" = plogis(apply(simplify2array(lapply(lys_cdf, qlogis)), 1:2, mean))
  )
}
