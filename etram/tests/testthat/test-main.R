
context("Transformation ensembles")

set.seed(3627)

nn <- function(output_shape = NULL, mbl = FALSE, input_shape = c(28, 28, 1), ...) {
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
im <- im[1:20,,,, drop = FALSE]
x1 <- rnorm(20)
x2 <- rnorm(20)
y <- y_train[1:20]
y_bin <- rbinom(20, 1, 0.6)
tab_dat <- data.frame(y = y,
                      y_bin = y_bin,
                      x1 = x1,
                      x2 = x2)
tab_dat$y <- factor(tab_dat$y, ordered = TRUE)
tab_dat$y_bin <- factor(tab_dat$y_bin, ordered = TRUE)

get_fml <- function(type = c("ordinal", "binary"),
                    mod = c("silscs", "sics", "cils", "ci")) {
  type <- match.arg(type)
  mod <- match.arg(mod)

  if (type == "ordinal") {
    if (mod %in% c("silscs", "cils")) {
      fml <- y ~ x1 + x2
    } else {
      fml <- y ~ 1
    }
  } else if (type == "binary") {
    if (mod %in% c("silscs", "cils")) {
      fml <- y_bin ~ x1 + x2
    } else {
      fml <- y_bin ~ 1
    }
  }
  return(fml)
}

test_that("ensemble works", {
  types <- c("ordinal", "binary")
  mods <- c("silscs", "sics", "cils", "ci")
  losses <- c("nll", "rps")
  lapply(types, function(ttype) {
    lapply(mods, function(tmod) {
      fml <- get_fml(type = ttype, mod = tmod)
      ridx <- save_ridx(nsplts = 2, prptest = 0.1, prpval = 0.1,
                        tab_dat = tab_dat, fml = fml, out_dir = "./", fname = "ttest")
      ridx <- get_ridx(in_dir = "./", fname = "ttest")
      expect_false(any(ridx[ridx$spl == 1 & ridx$type == "test", "idx"] %in%
                      c(ridx[ridx$spl == 1 & ridx$type == "val", "idx"],
                        ridx[ridx$spl == 1 & ridx$type == "train", "idx"]) |

                        ridx[ridx$spl == 2 & ridx$type == "test", "idx"] %in%
                        c(ridx[ridx$spl == 2 & ridx$type == "val", "idx"],
                          ridx[ridx$spl == 2 & ridx$type == "train", "idx"])
                      ))

      ensemble(mod = tmod, fml = fml, tab_dat = tab_dat, im = im, ridx = ridx,
               splits = 2, ensembles = 2, nn = nn, input_shape = c(28, 28, 1),
               bs = 5, lr = 5e-5, epochs = 2, loss = "nll", ws = TRUE,
               augment = FALSE, out_dir = "./", fname = "ttest")

      cdf_test <- list_cdfs(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "test")
      cdf_val <- list_cdfs(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "val")
      cdf_train <- list_cdfs(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "train")
      trafo_test <- list_trafos(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "test")
      trafo_val <- list_trafos(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "val")
      trafo_train <- list_trafos(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "train")
      terms_test <- list_terms(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "test")
      terms_val <- list_terms(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "val")
      terms_train <- list_terms(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2, type = "train")
      cdf1_test <- cdf_test[[1]][[1]]
      cdf1_val <- cdf_test[[1]][[1]]
      cdf1_train <- cdf_test[[1]][[1]]
      trafo1_test <- trafo_test[[1]][[1]]
      trafo1_val <- trafo_test[[1]][[1]]
      trafo1_train <- trafo_test[[1]][[1]]
      terms1_test <- terms_test[[1]][[1]]
      terms1_val <- terms_test[[1]][[1]]
      terms1_train <- terms_test[[1]][[1]]
      K <- ncol(cdf1_test)

      expect_true(is.list(cdf_test) & is.list(cdf_val) & is.list(cdf_train))
      expect_true(all(cdf1_test[, K] == 1) &
                  all(cdf1_val[, K] == 1) &
                  all(cdf1_train[, K] == 1))
      expect_false(any(is.na(cdf1_test)) &
                   any(is.na(cdf1_val)) &
                   any(is.na(cdf1_train)))

      expect_true(is.list(trafo_test) & is.list(trafo_val) & is.list(trafo_train))
      expect_true(ncol(trafo1_test) == K - 1 &
                  ncol(trafo1_val) == K - 1 &
                  ncol(trafo1_train) == K - 1)
      expect_false(any(is.na(trafo1_test)) &
                   any(is.na(trafo1_val)) &
                   any(is.na(trafo1_train)))

      expect_true(is.list(terms_test) & is.list(terms_val) & is.list(terms_train))
      expect_true(ncol(terms1_test) == K  &
                  ncol(terms1_val) == K &
                  ncol(terms1_train) == K)
      expect_false(any(is.na(terms1_test)) &
                   any(is.na(terms1_val)) &
                   any(is.na(terms1_train)))

      if (tmod %in% c("silscs", "cils")) {
        lors <- list_lors(in_dir = "./", fname = "ttest", splits = 2, ensembles = 2)
        lor1 <- lors[[1]][[1]]
        expect_true(is.list(lors))
        expect_true(ncol(lor1) == 2)
        expect_false(is.null(colnames(lor1)))
        expect_false(any(is.na(lor1)))
      }
    })
  })
})

unlink("ttest*.csv")
unlink("ckpts", recursive = TRUE)

