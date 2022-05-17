# Merge all val, test CDF of all models
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(etram)

# Directories -------------------------------------------------------------

in_dir <- out_dir <- "experiments/results/DE/MNIST/"

# Params ------------------------------------------------------------------

splits <- 6
ensembles <- 5

fname_cinll <- "mnist_ci_lossnll_wsno_augno"
fname_cirps <- "mnist_ci_lossrps_wsno_augno"

# Functions ---------------------------------------------------------------

bindr <- function(pat1 = "met", pat2 = NULL) {
  obj <- ls(pattern = pat1, envir = .GlobalEnv)
  if (!is.null(pat2)) {
    obj <- grep(pattern = pat2, x = obj, value = TRUE)
  }
  ret <- bind_rows(mget(obj, envir = .GlobalEnv))
  return(ret)
}

# Read data ---------------------------------------------------------------

dat <- load_data("mnist")
tab_dat <- dat$tab_dat

ridx <- get_ridx(in_dir, fname = "mnist")

# Combine single CDFs -----------------------------------------------------

### CI

## NLL

cdfval_cinll <- list_cdfs(in_dir, fname_cinll, splits, ensembles, "val")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdfval_cinll[[s]][[m]]$mem <- m
  }
  cdfval_cinll[[s]] <- do.call("rbind", cdfval_cinll[[s]])
  cdfval_cinll[[s]]$spl <- s
  cdfval_cinll[[s]]$type <- "val"
  cdfval_cinll[[s]]$loss <- "nll"
  cdfval_cinll[[s]]$mod <- "ci"
}
cdfval_cinll <- do.call("rbind", cdfval_cinll)

cdftest_cinll <- list_cdfs(in_dir, fname_cinll, splits, ensembles, "test")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdftest_cinll[[s]][[m]]$mem <- m
  }
  cdftest_cinll[[s]] <- do.call("rbind", cdftest_cinll[[s]])
  cdftest_cinll[[s]]$spl <- s
  cdftest_cinll[[s]]$type <- "test"
  cdftest_cinll[[s]]$loss <- "nll"
  cdftest_cinll[[s]]$mod <- "ci"
}
cdftest_cinll <- do.call("rbind", cdftest_cinll)

## RPS

cdfval_cirps <- list_cdfs(in_dir, fname_cirps, splits, ensembles, "val")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdfval_cirps[[s]][[m]]$mem <- m
  }
  cdfval_cirps[[s]] <- do.call("rbind", cdfval_cirps[[s]])
  cdfval_cirps[[s]]$spl <- s
  cdfval_cirps[[s]]$type <- "val"
  cdfval_cirps[[s]]$loss <- "rps"
  cdfval_cirps[[s]]$mod <- "ci"
}
cdfval_cirps <- do.call("rbind", cdfval_cirps)

cdftest_cirps <- list_cdfs(in_dir, fname_cirps, splits, ensembles, "test")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdftest_cirps[[s]][[m]]$mem <- m
  }
  cdftest_cirps[[s]] <- do.call("rbind", cdftest_cirps[[s]])
  cdftest_cirps[[s]]$spl <- s
  cdftest_cirps[[s]]$type <- "test"
  cdftest_cirps[[s]]$loss <- "rps"
  cdftest_cirps[[s]]$mod <- "ci"
}
cdftest_cirps <- do.call("rbind", cdftest_cirps)

cdf_ci <- bind_rows(cdfval_cinll, cdfval_cirps, cdftest_cinll, cdftest_cirps)
write.csv(cdf_ci, paste0(out_dir, "mnist_merged_cdf_ci.csv"), row.names = FALSE)

# Combine true Y ----------------------------------------------------------

y <- model.matrix(~ 0 + y, data = tab_dat)

yval_1 <- y[ridx[ridx$spl == 1 & ridx$type == "val", "idx"], ]
yval_2 <- y[ridx[ridx$spl == 2 & ridx$type == "val", "idx"], ]
yval_3 <- y[ridx[ridx$spl == 3 & ridx$type == "val", "idx"], ]
yval_4 <- y[ridx[ridx$spl == 4 & ridx$type == "val", "idx"], ]
yval_5 <- y[ridx[ridx$spl == 5 & ridx$type == "val", "idx"], ]
yval_6 <- y[ridx[ridx$spl == 6 & ridx$type == "val", "idx"], ]
y_true_val_all <- list(yval_1, yval_2, yval_3, yval_4, yval_5, yval_6)

for (s in seq_len(splits)) {
  y_true_val_all[[s]] <- as.data.frame(y_true_val_all[[s]])
  y_true_val_all[[s]]$spl <- s
  y_true_val_all[[s]]$type <- "val"
}
y_true_val_all <- do.call("rbind", y_true_val_all)

ytest_1 <- y[ridx[ridx$spl == 1 & ridx$type == "test", "idx"], ]
ytest_2 <- y[ridx[ridx$spl == 2 & ridx$type == "test", "idx"], ]
ytest_3 <- y[ridx[ridx$spl == 3 & ridx$type == "test", "idx"], ]
ytest_4 <- y[ridx[ridx$spl == 4 & ridx$type == "test", "idx"], ]
ytest_5 <- y[ridx[ridx$spl == 5 & ridx$type == "test", "idx"], ]
ytest_6 <- y[ridx[ridx$spl == 6 & ridx$type == "test", "idx"], ]
y_true_test_all <- list(ytest_1, ytest_2, ytest_3, ytest_4, ytest_5, ytest_6)

for (s in seq_len(splits)) {
  y_true_test_all[[s]] <- as.data.frame(y_true_test_all[[s]])
  y_true_test_all[[s]]$spl <- s
  y_true_test_all[[s]]$type <- "test"
}
y_true_test_all <- do.call("rbind", y_true_test_all)

# Combine -----------------------------------------------------------------

y_all <- bindr("y_true")
write.csv(y_all, paste0(out_dir, "mnist_merged_y.csv"), row.names = FALSE)
