### Benchmark transformation ensembles
### Lucas Kook, March 2023

set.seed(24101968)

# Dependencies ------------------------------------------------------------

library("caret")
library("parallel")
library("deeptrafo")
library("tidyverse")
library("xtable")

# Params ------------------------------------------------------------------

ncores <- 5
nfolds <- 5
nens <- 5
nep <- 1e4

odir <- file.path("benchmark", "results", make.names(Sys.time()))

### Create output directory
if(!dir.exists(odir))
  dir.create(odir, recursive = TRUE)

### Data sets to iterate over
datas <- c("airfoil", "concrete", "diabetes", "energy",
           "fish", "forest_fire", "ltfsid",
           "naval_compressor", "naval_turbine",
           "real", "wine", "yacht")

# FUNs --------------------------------------------------------------------

### lapply for loops
mylapply <- lapply
# mylapply <- function(...) mclapply(..., mc.cores = ncores)

### Load data
data_reader <- function(
    name = c("airfoil", "concrete", "diabetes", "energy", "fish", "forest_fire",
             "ltfsid", "naval_compressor", "naval_turbine", "real", "wine",
             "yacht")
) {
    name <- match.arg(name)
    read.table(paste0("benchmark/benchmark_data/", name, ".data"))
}

### Benchmark
benchmark_per_dataset <- function(name, folds = nfolds){

  data <- data_reader(name)

  X <- data[,1:(ncol(data)-1)]
  # Exclude columns with too few unique values
  X <- X[,which(apply(X, 2, function(x) length(unique(x))) > 1)]
  # Scale
  X <- scale(X)
  y <- data[,ncol(data)]

  # Wine has an ordinal outcome
  if (name == "wine")
    y <- ordered(y)

  # CV folds
  set.seed(1)
  folds <- createFolds(y, k = folds)

  # Fit all folds in parallel
  res <- mylapply(folds, function(testind) {

    # Prepare data
    trainind <- setdiff(1:nrow(X), testind)
    trainX <- X[trainind,]
    trainY <- if (name == "wine") y[trainind] else as.numeric(y[trainind])
    testX <- X[testind,]
    testY <- if (name == "wine") y[testind] else as.numeric(y[testind])

    dtrain <- list(y = trainY, x = trainX)
    dtest <- list(y = testY, x = testX)

    # Define neural network architecture
    nn <- \(x) x %>%
      layer_dense(input_shape = ncol(trainX), units = 16L, activation = "relu") %>%
      layer_dense(input_shape = ncol(trainX), units = 16L, activation = "relu") %>%
      layer_dense(input_shape = ncol(trainX), units = 8L, activation = "relu") %>%
      layer_dense(1L)

    # Fit transformation ensemble with callbacks
    ens <- trafoensemble(
      y ~ nn(x), data = dtrain, list_of_deep_models = list(nn = nn),
      validation_split = 0.1, optimizer = optimizer_adam(1e-3, decay = 1e-4),
      callbacks = list(callback_early_stopping(patience = 50, restore_best_weights = TRUE)),
      epochs = nep, n_ensemble = nens, seed = sample.int(1e6, 5), verbose = TRUE
    )

    lltrain <- logLik(ens, newdata = dtrain, convert_fun = mean)
    preds <- do.call("cbind", predict(ens, newdata = dtrain, type = "pdf"))
    lltrain$classical <- -mean(log(apply(preds, 1, mean)))

    lltest <- logLik(ens, newdata = dtest, convert_fun = mean)
    test_preds <- do.call("cbind", predict(ens, newdata = dtest, type = "pdf"))
    lltest$classical <- -mean(log(apply(test_preds, 1, mean)))

    # Extract NLL from train and test data
    res_fold_i <- suppressMessages(suppressWarnings(
      data.frame(cbind(train = unlist(lltrain),
                       test = unlist(lltest))) %>%
        rownames_to_column("ens")
    ))

  })

  res

}

# Run ---------------------------------------------------------------------

for (nam in datas) {
  res <- benchmark_per_dataset(nam)
  saveRDS(res, file = file.path(odir, paste0(nam, ".RDS")))
}

# Results -----------------------------------------------------------------

ndigits <- 2
lf <- list.files(odir)
tab_raw <- do.call("rbind", lapply(seq_along(lf), function(i) {

  table_for_data_i <- bind_rows(readRDS(file.path(odir, lf[i])), .id = "fold")
  sumtab <- table_for_data_i %>%
    group_by(ens) %>%
    summarize_at(c("train", "test"), list(mean = mean, sd = sd))
  dftrain <- as.data.frame(t(paste0(signif(sumtab$train_mean, ndigits), " (",
                               signif(sumtab$train_sd, ndigits), ")")))
  dftest <- as.data.frame(t(paste0(signif(sumtab$test_mean, ndigits), " (",
                               signif(sumtab$test_sd, ndigits), ")")))
  rownames(dftrain) <- rownames(dftest) <- gsub("(.*)\\.RDS", "\\1", lf[i])
  colnames(dftrain) <- colnames(dftest) <- sumtab$ens
  dftest
}))

dsn <- tools::toTitleCase(rownames(tab_raw))
dsn[dsn=="Forest_fire"] <- "ForestF"
dsn[dsn=="Naval_compressor"] <- "NavalC"
dsn[dsn=="Naval_turbine"] <- "NavalT"

rownames(tab_raw) <- dsn
tab_raw %>% xtable()
