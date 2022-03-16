
context("helper")

set.seed(6742)

nn <- function(output_shape = NULL, mbl = FALSE, input_shape = dim(im)[2:4], ...) {
  m <- keras_model_sequential() %>%
    layer_conv_2d(filters = 16,
                  input_shape = input_shape, kernel_size = c(3, 3), activation = "relu") %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_flatten() %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(output_shape, activation = "linear", use_bias = FALSE) %>%
    {if (mbl){
      ontram::layer_trafo_intercept()(.)
    } else .
    }
  return(m)
}

mnist <- dataset_mnist()
c(c(x_train, y_train), c(x_test, y_test)) %<-% mnist
x_train <- array_reshape(x_train, c(60000, 28, 28, 1))
im <- x_train / 255
im_train <- im[1:50,,,, drop = FALSE]
im_val <- im[51:70,,,, drop = FALSE]
im_test <- im[71:90,,,, drop = FALSE]
x_train <- data.matrix(rnorm(50))
x_val <- data.matrix(rnorm(20))
x_test <- data.matrix(rnorm(20))
y_train <- y_train[1:50]
y_train_bin <- rbinom(50, 1, 0.6)


test_that("get_input works", {
  mods <- c("silscs", "sics", "cils", "ci")
  lapply(mods, function(mod) {
    inp <- get_input(mod = mod, x_train = x_train, x_val = x_val, x_test = x_test,
                     im_train = im_train, im_val = im_val, im_test = im_test)
    inp_train <- inp$inp_train
    inp_val <- inp$inp_val
    inp_test <- inp$inp_test

    if (mod == "silscs") {
      expect_true(length(inp_train) == 3 &
                  length(inp_val) == 3 &
                  length(inp_test) == 3)
    } else if (mod %in% c("sics", "cils")) {
      expect_true(length(inp_train) == 2 &
                    length(inp_val) == 2 &
                    length(inp_test) == 2)
    }
    if (mod %in% c("silscs", "sics")) {
      expect_true(all(inp_train[[1]] == 1) &
                  all(inp_val[[1]] == 1) &
                  all(inp_test[[1]] == 1))
    } else if (mod == "ci") {
      expect_true(all(inp_train[[2]] == 1) &
                  all(inp_val[[2]] == 1) &
                  all(inp_test[[2]] == 1) &
                  length(inp_train) == 2 &
                  length(inp_val) == 2 &
                  length(inp_test) == 2)
    }
  })
})

test_that("get_model works", {
  mods <- c("silscs", "sics", "cils", "ci")
  lapply(mods, function(mod) {
    m <- get_model(mod, 10L, x_dim = 1L, nn = nn, input_shape = c(28L, 28L, 1L))
    expect_true("k_ontram" %in% class(m))
    if (mod %in% c("cils", "ci")) {
      expect_true("k_ontram_ci" %in% class(m))
    }
  })
})


test_that("warm_mod works for oridinal responses", {
  mods <- c("silscs", "sics", "cils")
  expect_error(warm_mod("ci", x = x_train, binary = FALSE))
  lapply(mods, function(mod) {
    m <- get_model(mod, 10L, x_dim = 1L, nn = nn, input_shape = c(28L, 28L, 1L))
    m <- warm_mod(m, mod, x = x_train, y = to_categorical(y_train), binary = FALSE)
    df <- data.frame(y = y_train, x = x_train)
    df$y <- factor(df$y, ordered = TRUE)
    theta <- coef(Polr(y ~ x, data = df), with_baseline = TRUE)[1:9]
    theta_uncond <- coef(Polr(y ~ 1, data = df), with_baseline = TRUE)[1:9]
    beta <- coef(Polr(y ~ x, data = df))
    if (mod == "silscs") {
      expect_equal(as.numeric(get_weights(m$list_of_shift_models[[0]])[[1]]), as.numeric(beta))
      expect_equal(as.numeric(get_weights(m$mod_baseline)[[1]]), as.numeric(ontram:::.to_gamma(theta)),
                   tolerance = 1e-5)
    } else if (mod == "sics") {
      expect_equal(as.numeric(get_weights(m$mod_baseline)[[1]]), as.numeric(ontram:::.to_gamma(theta_uncond)),
                   tolerance = 1e-5)
    } else if (mod == "cils") {
      expect_equal(as.numeric(get_weights(m$list_of_shift_models)[[1]]), as.numeric(beta))
    }
  })
})

test_that("warm_mod works for binary responses", {
  mods <- c("silscs", "sics", "cils")
  expect_error(warm_mod("ci", x = x_train, binary = TRUE))
  lapply(mods, function(mod) {
    m <- get_model(mod, 2L, x_dim = 1L, nn = nn, input_shape = c(28L, 28L, 1L))
    m <- warm_mod(m, mod, x = x_train, y = to_categorical(y_train_bin), binary = TRUE)
    df <- data.frame(y = y_train_bin, x = x_train)
    df$y <- factor(df$y, levels = c(0, 1))
    theta <- coef(glm(y ~ x, data = df, family = "binomial"))[1]
    theta_uncond <- coef(glm(y ~ 1, data = df, family = "binomial"))[1]
    beta <- coef(glm(y ~ x, data = df, family = "binomial"))[2]
    if (mod == "silscs") {
      expect_equal(as.numeric(get_weights(m$list_of_shift_models[[0]])[[1]]), as.numeric(beta))
      expect_equal(as.numeric(get_weights(m$mod_baseline)[[1]]), as.numeric(ontram:::.to_gamma(theta)),
                   tolerance = 1e-5)
    } else if (mod == "sics") {
      expect_equal(as.numeric(get_weights(m$mod_baseline)[[1]]), as.numeric(ontram:::.to_gamma(theta_uncond)),
                   tolerance = 1e-5)
    } else if (mod == "cils") {
      expect_equal(as.numeric(get_weights(m$list_of_shift_models)[[1]]), as.numeric(beta))
    }

  })
})

cdf1 <- data.frame(matrix(c(0.1, 0.2, 1,
                            0.3, 0.7, 1,
                            0.5, 0.7, 1),
                          nrow = 3, byrow = TRUE))
cdf2 <- data.frame(matrix(c(0.2, 0.3, 1,
                            0.5, 0.8, 1,
                            0.2, 0.9, 1),
                          nrow = 3, byrow = TRUE))
cdf3 <- data.frame(matrix(c(0.7, 0.8, 1,
                            0.6, 0.8, 1,
                            0.2, 0.8, 1),
                          nrow = 3, byrow = TRUE))
y_true <- data.frame(matrix(c(0, 0, 1,
                              1, 0, 0,
                              0, 1, 0),
                            nrow = 3, byrow = TRUE))
lys_cdf <- list(cdf1, cdf2, cdf3)

test_that("get_ensemble works", {
  ttype <- c("linear", "log-linear", "trafo")
  lapply(ttype, function(t) {
    cdf <- get_ensemble(lys_cdf, type = t)
    expect_true(is.matrix(cdf))
    if (t == "linear") {
      expect_equal(as.numeric(cdf[1, 1]), mean(c(0.1, 0.2, 0.7)))
    } else if (t == "log-linear") {
      expect_equal(as.numeric(cdf[1, 1]), exp(mean(log(c(0.1, 0.2, 0.7)))))
    } else if (t == "trafo") {
      expect_equal(as.numeric(cdf[1, 1]), plogis(mean(qlogis(c(0.1, 0.2, 0.7)))))
    }
  })
})

test_that("get_order works", {
  met <- c("nll", "rps")
  lapply(met, function(m) {
    o <- get_order(lys_cdf, y_true, order_metric = m)
    if (m == "nll") {
      expect_equal(o, c(2, 3, 1))
    } else {
      expect_equal(o, c(2, 1, 3))
    }
  })
})
