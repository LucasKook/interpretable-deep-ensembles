#' Fit unconditional (SI) or linear shift model (SI-LSx) by minimizing the RPS
#' @param fml model formula (e.g. \code{y ~ 1} for fitting the unconditional model).
#' @param tab_dat data frame of tabular data (including response variable). Response variable needs to be ordered
#' (also for binary responses).
#' @param ridx row indices used for each split (output from \code{get_ridx}).
#' @param out_dir directory to save results.
#' @param fname unique file name.
#' @export
fit_ref_lossrps <- function(mod = c("sils", "si"), fml, tab_dat, ridx,
                            splits = 6,
                            bs = ncol(tab_dat), lr = 0.1, decay = 1e-4,
                            optimizer = optimizer_adam(learning_rate = lr, decay = decay),
                            epochs = 800,
                            out_dir, fname) {

  m_dir <- file.path(out_dir, "ckpts/")
  mod <- match.arg(mod)
  l <- k_ontram_rps(K)
  mf <- model.frame(fml, data = tab_dat)
  if (ncol(mf) != 1) { # check if tabular predictors
    x <- ontram:::.rm_int(model.matrix(fml, data = mf))
  } else {
    x <- NULL
  }
  y <- mf[, 1L, drop = FALSE]
  y_name <- colnames(y)
  stopifnot(is.ordered(tab_dat[, y_name])) # otherwise one-hot encoding won't work
  fmly <- ~ 0 + eval(parse(text = y_name))
  mfy <- model.frame(fmly, data = tab_dat)
  y <- model.matrix(fmly, data = mfy)

  for (spl in seq_len(splits)) {
    message("Split: ", spl)

    ### Train, test, validation set for current split ###

    rtrain <- ridx[ridx$spl == spl & ridx$type == "train", "idx"]
    rval <- ridx[ridx$spl == spl & ridx$type == "val", "idx"]
    rtest <- ridx[ridx$spl == spl & ridx$type == "test", "idx"]

    if (!is.null(x)) {
      x_train <- x[rtrain, , drop = FALSE]
      x_val <- x[rval, , drop = FALSE]
      x_test <- x[rtest, , drop = FALSE]
    } else {
      x_train <- NULL
      x_val <- NULL
      x_test <- NULL
    }
    y_train <- y[rtrain, , drop = FALSE]
    y_val <- y[rval, , drop = FALSE]
    y_test <- y[rtest, , drop = FALSE]

    K <- sum(colSums(y_train) != 0)

    ### Prepare input ###

    one_train <- matrix(1, nrow = nrow(y_train))
    one_val <- matrix(1, nrow = nrow(y_val))
    one_test <- matrix(1, nrow = nrow(y_test))

    if (mod == "si") {
      inp_train <- list(one_train, one_train)
      inp_val <- list(one_val, one_val)
      inp_test <- list(one_test, one_test)
    } else {
      inp_train <- list(one_train, x_train)
      inp_val <- list(one_val, x_val)
      inp_test <- list(one_test, x_test)
    }

    ### Prepare model ###

    if (mod == "si") {
      mod_si <- k_mod_baseline(K)
      m <- k_ontram(mod_si, complex_intercept = FALSE)
    } else {
      mod_si <- k_mod_baseline(K)
      mod_ls <- mod_shift(ncol(x_train))
      m <- k_ontram(mod_si, mod_ls, complex_intercept = FALSE)
    }
    m <- compile(m, optimizer = optimizer, loss = l)

    ### Fit model ###

    mpath <- file.path(m_dir, paste0(fname, "_m_spl", spl, ".hdf5"))
    hpath <- file.path(out_dir, paste0(fname, "_hist_spl", spl, ".csv"))

    h <- h <- fit(m,
                  x = inp_train, y = y_train,
                  validation_data = list(inp_val, y_val),
                  shuffle = TRUE, batch_size = bs, epochs = epochs,
                  view_metrics = TRUE)
    save_k_hist(h, hpath)
    save_model_weights_hdf5(m, mpath)

    ### Evaluate model ###

    ctrainpath <- file.path(out_dir, paste0(fname, "_cdftrain_spl", spl, ".csv"))
    cvalpath <- file.path(out_dir, paste0(fname, "_cdfval_spl", spl, ".csv"))
    ctestpath <- file.path(out_dir, paste0(fname, "_cdftest_spl", spl, ".csv"))
    ttrainpath <- file.path(out_dir, paste0(fname, "_trafotrain_spl", spl, ".csv"))
    tvalpath <- file.path(out_dir, paste0(fname, "_trafoval_spl", spl, ".csv"))
    ttestpath <- file.path(out_dir, paste0(fname, "_trafotest_spl", spl, ".csv"))
    rtrainpath <- file.path(out_dir, paste0(fname, "_rawtrain_spl", spl, ".csv"))
    rvalpath <- file.path(out_dir, paste0(fname, "_rawval_spl", spl, ".csv"))
    rtestpath <- file.path(out_dir, paste0(fname, "_rawtest_spl", spl, ".csv"))
    lorpath <- file.path(out_dir, paste0(fname, "_lor_spl", spl, ".csv"))

    #### CDFs

    cdf_train <- as.data.frame(predict(m, x = inp_train, type = "distribution", batch_size = bs))
    rownames(cdf_train) <- rtrain
    cdf_val <- as.data.frame(predict(m, x = inp_val, type = "distribution", batch_size = bs))
    rownames(cdf_val) <- rval
    cdf_test <- as.data.frame(predict(m, x = inp_test, type = "distribution", batch_size = bs))
    rownames(cdf_test) <- rtest
    write.csv(cdf_train, file = ctrainpath)
    write.csv(cdf_val, file = cvalpath)
    write.csv(cdf_test, file = ctestpath)

    #### Trafo (theta - xB - eta(B))

    trafo_train <- as.data.frame(predict(m, x = inp_train, type = "trafo", batch_size = bs))
    rownames(trafo_train) <- rtrain
    trafo_val <- as.data.frame(predict(m, x = inp_val, type = "trafo", batch_size = bs))
    rownames(trafo_val) <- rval
    trafo_test <- as.data.frame(predict(m, x = inp_test, type = "trafo", batch_size = bs))
    rownames(trafo_test) <- rtest
    write.csv(trafo_train, file = ttrainpath)
    write.csv(trafo_val, file = tvalpath)
    write.csv(trafo_test, file = ttestpath)

    #### Raw (theta, xB - eta(B))

    raw_train <- as.data.frame(predict(m, x = inp_train, type = "terms", batch_size = bs))
    rownames(raw_train) <- rtrain
    raw_val <- as.data.frame(predict(m, x = inp_val, type = "terms", batch_size = bs))
    rownames(raw_val) <- rval
    raw_test <- as.data.frame(predict(m, x = inp_test, type = "terms", batch_size = bs))
    rownames(raw_test) <- rtest
    write.csv(raw_train, file = rtrainpath)
    write.csv(raw_val, file = rvalpath)
    write.csv(raw_test, file = rtestpath)

    #### Log odds ratios

    if (mod == "sils") {
      lor <- as.data.frame(t(unlist(get_weights(m$list_of_shift_models))))
      colnames(lor) <- colnames(x_train)
      write.csv(lor, file = lorpath, row.names = FALSE)
    }
  }
}
