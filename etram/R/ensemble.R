
#' Perform ensembling experiment
#' @description Depending on the model and response type a \code{k_ontram}, \code{Polr} or \code{glm} model is fit.
#' The function handles ordinal and binary responses.
#' @param fml model formula (e.g. \code{y ~ 1} for fitting a CI model).
#' @param tab_dat data frame of tabular data (including response variable). Response variable needs to be ordered
#' (also for binary responses).
#' @param ridx row indices used for each split (output from \code{get_ridx}).
#' @param nn function to build neural network used for modeling complex intercept or complex shift term.
#' Arguments \code{input_shape}, \code{output_shape} and \code{mbl} are required (see \code{\link{cnn_stroke}} for an example).
#' @param single_split numeric. Number of split that should be fitted.
#' @param single_ens numeric. Number of ensemble member that should be fitted.
#' @param ws logical. Whether to initialize intercepts (simple and complex) and linear shift terms using the parameters estimated by a
#' \code{Polr} or \code{glm} model (depending on the response type).
#' @param augment logical. Whether to perform 2D image augmentation during training.
#' @param aug_params parameters used for image augmentation.
#' @param train_batchwise logical. If \code{TRUE} model is trained batchwise (used to save memory).
#' @param out_dir directory to save results.
#' @param fname unique file name.
#' @export
ensemble <- function(mod = c("silscs", "sics", "cils", "ci", "si", "sils"), fml, tab_dat, im = NULL, ridx,
                     splits = 6, ensembles = 5, single_split = NULL, single_ens = NULL,
                     nn = NULL, input_shape = NULL,
                     bs = 1, lr = 5e-5, optimizer = optimizer_adam(learning_rate = lr), epochs = 10,
                     loss = c("nll", "rps"),
                     ws = FALSE, augment = FALSE, aug_params = list(rotation_range = 20,
                                                                    width_shift_range = 0.2,
                                                                    height_shift_range = 0.2,
                                                                    zoom_range = 0.15,
                                                                    shear_range = 0.15,
                                                                    fill_mode = "nearest"),
                     train_batchwise = FALSE,
                     out_dir, fname) {

  # stopifnot(nrow(ridx) == (splits * nrow(tab_dat)))
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  if (mod %in% c("silscs", "sics", "cils", "ci") & !dir.exists(file.path(out_dir, "ckpts/"))) {
    dir.create(file.path(out_dir, "ckpts/"))
  }
  m_dir <- file.path(out_dir, "ckpts/")

  mod <- match.arg(mod)
  loss <- match.arg(loss)

  ### Split tab_dat into response and predictors ###

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

  if (is.null(single_split)) {
    start_spl <- 1
    end_spl <- splits
  } else {
    start_spl <- end_spl <- single_split
  }

  for (spl in start_spl:end_spl) {

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
    if (!is.null(im)) {
      im_train <- ontram:::.batch_subset(im, rtrain, dim(im))
      im_val <- ontram:::.batch_subset(im, rval, dim(im))
      im_test <- ontram:::.batch_subset(im, rtest, dim(im))
    }
    y_train <- y[rtrain, , drop = FALSE]
    y_val <- y[rval, , drop = FALSE]
    y_test <- y[rtest, , drop = FALSE]

    K <- sum(colSums(y_train) != 0)
    binary <- K == 2

    if (binary) { # otherwise if subset of ordinal response ncol(y) == K
      classes <- which(colSums(y_train) != 0)
      y_train <- y_train[, classes, drop = FALSE]
      y_val <- y_val[, classes, drop = FALSE]
      y_test <- y_test[, classes, drop = FALSE]
    }

    ### Fit ensembles ###

    if (!(mod %in% c("si", "sils"))) {

      if (is.null(single_ens)) {
        start_ens <- 1
        end_ens <- ensembles
      } else {
        start_ens <- end_ens <- single_ens
      }

      for (ens in start_ens:end_ens) {
        message("Split: ", spl, ", Ensemble: ", ens)
        k_clear_session()

        ### Prepare input ###

        inp <- get_input(mod = mod,
                         x_train = x_train, x_val = x_val, x_test = x_test,
                         im_train = im_train, im_val = im_val, im_test = im_test)
        inp_train <- inp$inp_train
        inp_val <- inp$inp_val
        inp_test <- inp$inp_test

        ### Prepare model ###

        m <- get_model(mod = mod,
                       y_dim = K, x_dim = ncol(x_train),
                       nn = nn, input_shape = input_shape)
        if (ws) {
          if (mod %in% c("silscs", "sics", "cils", "ci")) {
            m <- warm_mod(m, mod = mod, x = x_train, y = y_train, binary = binary)
          } else if (mod %in% c("si", "sils")) {
            message("Warmstart not possible for this model type")
          }
        }
        if (loss == "rps") {
          l <- k_ontram_rps(K)
        } else if (loss == "nll") {
          l <- k_ontram_loss(K)
        }
        m <- compile(m, optimizer = optimizer, loss = l)

        ### Fit model ###

        mpath <- file.path(m_dir, paste0(fname, "_m_spl", spl, "_ens", ens, ".hdf5"))
        hpath <- file.path(out_dir, paste0(fname, "_hist_spl", spl, "_ens", ens, ".csv"))

        if (!augment) {
          if (!train_batchwise) {
            h <- fit(m,
                     x = inp_train, y = y_train,
                     validation_data = list(inp_val, y_val),
                     shuffle = TRUE, batch_size = bs, epochs = epochs,
                     callbacks = list(callback_model_checkpoint(mpath,
                                                                monitor = "val_loss",
                                                                save_best_only = TRUE,
                                                                save_weights_only = TRUE)),
                     view_metrics = FALSE)
            save_k_hist(h, hpath)
            load_model_weights_hdf5(m, mpath)
          } else if (train_batchwise) {
            gen_train <- generator(mod = mod, im = im_train, x = x_train, y = y_train,
                                   batch_size = bs,
                                   shuffle = TRUE)
            gen_val <- generator(mod = mod, im = im_val, x = x_val, y = y_val,
                                 batch_size = bs,
                                 shuffle = TRUE)
            h <- fit(m,
                     x = gen_train,
                     validation_data = gen_val,
                     epochs = epochs,
                     steps_per_epoch = ceiling(nrow(y_train)/bs),
                     validation_steps = ceiling(nrow(y_val)/bs),
                     callbacks = list(callback_model_checkpoint(mpath,
                                                                monitor = "val_loss",
                                                                save_best_only = TRUE,
                                                                save_weights_only = TRUE)),
                     view_metrics = FALSE)
            save_k_hist(h, hpath)
            load_model_weights_hdf5(m, mpath)
          }

        } else if (augment) {
          im_gen <- do.call(image_data_generator, aug_params)
          if (mod %in% c("cils", "ci")) {
            mim_as_mbl <- TRUE
          } else {
            mim_as_mbl <- FALSE
          }
          f <- fit_k_ontram_augmented_data(m,
                                           im_train = im_train, im_val = im_val,
                                           x_train = x_train, x_val = x_val,
                                           y_train = y_train, y_val = y_val,
                                           generator = im_gen, epochs = epochs, bs = bs,
                                           mim_as_mbl = mim_as_mbl,
                                           patience = 1, filepath = mpath)
          h <- f$hist
          save_k_hist(h, hpath)
          load_model_weights_hdf5(m, mpath)
        }

        ### Evaluate model ###

        ctrainpath <- file.path(out_dir, paste0(fname, "_cdftrain_spl", spl, "_ens", ens, ".csv"))
        cvalpath <- file.path(out_dir, paste0(fname, "_cdfval_spl", spl, "_ens", ens, ".csv"))
        ctestpath <- file.path(out_dir, paste0(fname, "_cdftest_spl", spl, "_ens", ens, ".csv"))
        ttrainpath <- file.path(out_dir, paste0(fname, "_trafotrain_spl", spl, "_ens", ens, ".csv"))
        tvalpath <- file.path(out_dir, paste0(fname, "_trafoval_spl", spl, "_ens", ens, ".csv"))
        ttestpath <- file.path(out_dir, paste0(fname, "_trafotest_spl", spl, "_ens", ens, ".csv"))
        rtrainpath <- file.path(out_dir, paste0(fname, "_rawtrain_spl", spl, "_ens", ens, ".csv"))
        rvalpath <- file.path(out_dir, paste0(fname, "_rawval_spl", spl, "_ens", ens, ".csv"))
        rtestpath <- file.path(out_dir, paste0(fname, "_rawtest_spl", spl, "_ens", ens, ".csv"))
        lorpath <- file.path(out_dir, paste0(fname, "_lor_spl", spl, "_ens", ens, ".csv"))

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

        if (mod %in% c("silscs", "cils")) {
          if (mod == "silscs") {
            lor <- as.data.frame(t(unlist(get_weights(m$list_of_shift_models[[0]]))))
            colnames(lor) <- colnames(x_train)
          } else if (mod == "cils") {
            lor <- as.data.frame(t(unlist(get_weights(m$list_of_shift_models))))
            colnames(lor) <- colnames(x_train)
          }
          write.csv(lor, file = lorpath, row.names = FALSE)
        }
        gc() # garbage collection
        rm(m, h)
      }
    } else if (mod %in% c("si", "sils")) {

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

      #### Prepare data sets to fit Polr and GLM

      # train
      ytab <- apply(y_train, 1, which.max)
      df_train <- as.data.frame(cbind(ytab, x_train))
      # val
      ytab <- apply(y_val, 1, which.max)
      df_val <- as.data.frame(cbind(ytab, x_val))
      # test
      ytab <- apply(y_test, 1, which.max)
      df_test <- as.data.frame(cbind(ytab, x_test))

      if (!binary) {

        #### Y as ordered factor

        df_train$ytab <- factor(df_train$ytab, ordered = TRUE)
        df_val$ytab <- factor(df_val$ytab, ordered = TRUE)
        df_test$ytab <- factor(df_test$ytab, ordered = TRUE)

        if (mod == "si") {

          m <- Polr(ytab ~ 1, data = df_train)

          #### CDFs

          # cdf will be the same for train, val, test data (prevalence)
          cdf_train <- as.data.frame(cbind(0, matrix(rep(predict(m, newdata = df_train[, -1], type = "distribution"),
                                                         each = nrow(y_train)), nrow = nrow(y_train), byrow = FALSE)))
          rownames(cdf_train) <- rtrain
          cdf_val <- as.data.frame(cbind(0, matrix(rep(predict(m, newdata = df_val[, -1], type = "distribution"),
                                                       each = nrow(y_val)), nrow = nrow(y_val), byrow = FALSE)))
          rownames(cdf_val) <- rval
          cdf_test <- as.data.frame(cbind(0, matrix(rep(predict(m, newdata = df_test[, -1], type = "distribution"),
                                                        each = nrow(y_test)), nrow = nrow(y_test), byrow = FALSE)))
          rownames(cdf_test) <- rtest

          #### Trafo (theta)

          thetas <- coef(m, with_baseline = TRUE)[1:K-1]
          trafo_train <- as.data.frame(matrix(rep(thetas, each = nrow(df_train)), nrow = nrow(df_train)))
          rownames(trafo_train) <- rtrain

          trafo_val <- as.data.frame(matrix(rep(thetas, each = nrow(df_val)), nrow = nrow(df_val)))
          rownames(trafo_val) <- rval

          trafo_test <- as.data.frame(matrix(rep(thetas, each = nrow(df_test)), nrow = nrow(df_test)))
          rownames(trafo_test) <- rtest

          #### Raw (theta, 0)

          raw_train <- cbind(trafo_train, 0)
          rownames(raw_train) <- rtrain

          raw_val <- cbind(trafo_val, 0)
          rownames(raw_val) <- rval

          raw_test <- cbind(trafo_test, 0)
          rownames(raw_test) <- rtest

        } else if (mod == "sils") {

          m <- Polr(ytab ~ ., data = df_train)

          #### CDFs

          cdf_train <- as.data.frame(cbind(0, t(predict(m, newdata = df_train[, -1, drop = FALSE], type = "distribution"))))
          rownames(cdf_train) <- rtrain
          cdf_val <- as.data.frame(cbind(0, t(predict(m, newdata = df_val[, -1, drop = FALSE], type = "distribution"))))
          rownames(cdf_val) <- rval
          cdf_test <- as.data.frame(cbind(0, t(predict(m, newdata = df_test[, -1, drop = FALSE], type = "distribution"))))
          rownames(cdf_test) <- rtest

          #### Trafo (theta - xB)

          thetas <- coef(m, with_baseline = TRUE)[1:K-1]
          b <- coef(m)
          thetas_train <- matrix(rep(thetas, each = nrow(df_train)), nrow = nrow(df_train))
          xb_train <- x_train %*% b
          trafo_train <- as.data.frame(apply(thetas_train, 2, function(x) x - xb_train))
          rownames(trafo_train) <- rtrain

          thetas_val <- matrix(rep(thetas, each = nrow(df_val)), nrow = nrow(df_val))
          xb_val <- x_val %*% b
          trafo_val <- as.data.frame(apply(thetas_val, 2, function(x) x - xb_val))
          rownames(trafo_val) <- rval

          thetas_test <- matrix(rep(thetas, each = nrow(df_test)), nrow = nrow(df_test))
          xb_test <- x_test %*% b
          trafo_test <- as.data.frame(apply(thetas_test, 2, function(x) x - xb_test))
          rownames(trafo_test) <- rtest

          #### Raw (theta, xB)

          raw_train <- as.data.frame(cbind(thetas_train, xb_train))
          rownames(raw_train) <- rtrain

          raw_val <- as.data.frame(cbind(thetas_val, xb_val))
          rownames(raw_val) <- rval

          raw_test <- as.data.frame(cbind(thetas_test, xb_test))
          rownames(raw_test) <- rtest

          #### Log odds ratios

          lor <- as.data.frame(t(coef(m)))
          colnames(lor) <- colnames(x_train)
          write.csv(lor, file = lorpath, row.names = FALSE)

        }
      } else if (binary) {

        #### Prepare data sets to fit GLM

        df_train$ytab <- factor(df_train$ytab, levels = c(1, 2))
        df_val$ytab <- factor(df_val$ytab, levels = c(1, 2))
        df_test$ytab <- factor(df_test$ytab, levels = c(1, 2))

        if (mod == "si") {

          m <- glm(ytab ~ 1, family = "binomial", data = df_train)

          #### CDFs

          p1 <- predict(m, newdata = df_train, type = "response")
          cdf_train <- as.data.frame(cbind(0, 1 - p1, 1))
          rownames(cdf_train) <- rtrain

          p1 <- predict(m, newdata = df_val, type = "response")
          cdf_val <- as.data.frame(cbind(0, 1 - p1, 1))
          rownames(cdf_val) <- rval

          p1 <- predict(m, newdata = df_test, type = "response")
          cdf_test <- as.data.frame(cbind(0, 1 - p1, 1))
          rownames(cdf_test) <- rtest

          #### Trafo (theta)

          theta <- coef(m)[1]
          trafo_train <- as.data.frame(rep(theta, nrow(df_train)))
          rownames(trafo_train) <- rtrain
          trafo_val <- as.data.frame(rep(theta, nrow(df_val)))
          rownames(trafo_val) <- rval
          trafo_test <- as.data.frame(rep(theta, nrow(df_test)))
          rownames(trafo_test) <- rtest

          #### Raw (theta, 0)

          raw_train <- cbind(trafo_train, 0)
          rownames(raw_train) <- rtrain
          raw_val <- cbind(trafo_val, 0)
          rownames(raw_val) <- rval
          raw_test <- cbind(trafo_test, 0)
          rownames(raw_test) <- rtest

        } else if (mod == "sils") {

          m <- glm(ytab ~ ., family = "binomial", data = df_train)

          #### CDFs

          p1 <- predict(m, newdata = df_train, type = "response")
          cdf_train <- as.data.frame(cbind(0, 1 - p1, 1))
          rownames(cdf_train) <- rtrain

          p1 <- predict(m, newdata = df_val, type = "response")
          cdf_val <- as.data.frame(cbind(0, 1 - p1, 1))
          rownames(cdf_val) <- rval

          p1 <- predict(m, newdata = df_test, type = "response")
          cdf_test <- as.data.frame(cbind(0, 1 - p1, 1))
          rownames(cdf_test) <- rtest

          #### Trafo (theta - xB)

          theta <- coef(m)[1]
          theta_train <- as.data.frame(rep(theta, nrow(df_train)))
          theta_val <- as.data.frame(rep(theta, nrow(df_val)))
          theta_test <- as.data.frame(rep(theta, nrow(df_test)))

          trafo_train <- as.data.frame(theta_train - predict(m, newdata = df_train, type = "link"))
          rownames(trafo_train) <- rtrain

          trafo_val <- as.data.frame(theta_val - predict(m, newdata = df_val, type = "link"))
          rownames(trafo_val) <- rval

          trafo_test <- as.data.frame(theta_test - predict(m, newdata = df_test, type = "link"))
          rownames(trafo_test) <- rtest

          #### Raw (theta, xB)

          xb_train <- predict(m, newdata = df_train, type = "link")
          raw_train <- as.data.frame(cbind(theta_train, xb_train))
          rownames(raw_train) <- rtrain

          xb_val <- predict(m, newdata = df_val, type = "link")
          raw_val <- as.data.frame(cbind(theta_val, xb_val))
          rownames(raw_val) <- rval

          xb_test <- predict(m, newdata = df_test, type = "link")
          raw_test <- as.data.frame(cbind(theta_test, xb_test))
          rownames(raw_test) <- rtest

          #### Log odds ratio

          lor <- as.data.frame(t(coef(m)[2:length(coef(m))]))
          colnames(lor) <- colnames(x_train)
          write.csv(lor, file = lorpath, row.names = FALSE)
        }
      }
      write.csv(cdf_train, file = ctrainpath)
      write.csv(cdf_val, file = cvalpath)
      write.csv(cdf_test, file = ctestpath)
      write.csv(trafo_train, file = ttrainpath)
      write.csv(trafo_val, file = tvalpath)
      write.csv(trafo_test, file = ttestpath)
      write.csv(raw_train, file = rtrainpath)
      write.csv(raw_val, file = rvalpath)
      write.csv(raw_test, file = rtestpath)
    }
    rm(y_train, y_val, y_test)
    if (!is.null(x)) {
      rm(x_train, x_val, x_test)
    }
    if (!is.null(im)) {
      rm(im_train, im_val, im_test)
    }
  }
}
