# Merge all val, test CDF of all models
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(readr)
library(rhdf5)
library(etram)

# Directories -------------------------------------------------------------

im_path <- "~/../data/mela_all/data/train_mela_images/"
path <- "~/../data/mela_all/data/train.csv"
in_dir <- out_dir <- "experiments/results/DE/melanoma/"

# Params ------------------------------------------------------------------

splits <- 6
ensembles <- 5

fname_cilsnll <- "mela_silscs_lossnll_wsyes_augno"
fname_cilsrps <- "mela_silscs_lossrps_wsyes_augno"
fname_cinll <- "mela_ci_lossnll_wsyes_augno"
fname_cirps <- "mela_ci_lossrps_wsyes_augno"

fname_sinll <- "mela_si"
fname_silsnll <- "mela_sils"
fname_sirps <- "mela_si_rps"
fname_silsrps <- "mela_sils_rps"

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

dat <- load_data("melanoma", path = path, im_path = im_path)
tab_dat <- dat$tab_dat

ridx <- get_ridx(in_dir, fname = "melanoma")

# Combine single CDFs -----------------------------------------------------

### CI-LS

## NLL

cdfval_cilsnll <- list_cdfs(in_dir, fname_cilsnll, splits, ensembles, "val")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdfval_cilsnll[[s]][[m]]$mem <- m
  }
  cdfval_cilsnll[[s]] <- do.call("rbind", cdfval_cilsnll[[s]])
  cdfval_cilsnll[[s]]$spl <- s
  cdfval_cilsnll[[s]]$type <- "val"
  cdfval_cilsnll[[s]]$loss <- "nll"
  cdfval_cilsnll[[s]]$mod <- "cils"
}
cdfval_cilsnll <- do.call("rbind", cdfval_cilsnll)

cdftest_cilsnll <- list_cdfs(in_dir, fname_cilsnll, splits, ensembles, "test")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdftest_cilsnll[[s]][[m]]$mem <- m
  }
  cdftest_cilsnll[[s]] <- do.call("rbind", cdftest_cilsnll[[s]])
  cdftest_cilsnll[[s]]$spl <- s
  cdftest_cilsnll[[s]]$type <- "test"
  cdftest_cilsnll[[s]]$loss <- "nll"
  cdftest_cilsnll[[s]]$mod <- "cils"
}
cdftest_cilsnll <- do.call("rbind", cdftest_cilsnll)

## RPS

cdfval_cilsrps <- list_cdfs(in_dir, fname_cilsrps, splits, ensembles, "val")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdfval_cilsrps[[s]][[m]]$mem <- m
  }
  cdfval_cilsrps[[s]] <- do.call("rbind", cdfval_cilsrps[[s]])
  cdfval_cilsrps[[s]]$spl <- s
  cdfval_cilsrps[[s]]$type <- "val"
  cdfval_cilsrps[[s]]$loss <- "rps"
  cdfval_cilsrps[[s]]$mod <- "cils"
}
cdfval_cilsrps <- do.call("rbind", cdfval_cilsrps)

cdftest_cilsrps <- list_cdfs(in_dir, fname_cilsrps, splits, ensembles, "test")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdftest_cilsrps[[s]][[m]]$mem <- m
  }
  cdftest_cilsrps[[s]] <- do.call("rbind", cdftest_cilsrps[[s]])
  cdftest_cilsrps[[s]]$spl <- s
  cdftest_cilsrps[[s]]$type <- "test"
  cdftest_cilsrps[[s]]$loss <- "rps"
  cdftest_cilsrps[[s]]$mod <- "cils"
}
cdftest_cilsrps <- do.call("rbind", cdftest_cilsrps)

cdf_cils <- bind_rows(cdfval_cilsnll, cdfval_cilsrps, cdftest_cilsnll, cdftest_cilsrps)
write.csv(cdf_cils, paste0(out_dir, "mela_merged_cdf_cils.csv"), row.names = FALSE)


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
write.csv(cdf_ci, paste0(out_dir, "mela_merged_cdf_ci.csv"), row.names = FALSE)


### SI

## NLL

cdfval_sinll <- list_cdfs(in_dir, fname_sinll, splits, ensembles = NULL, "val")
for (s in seq_len(splits)) {
  for (m in 1) {
    cdfval_sinll[[s]][[m]]$mem <- m
  }
  cdfval_sinll[[s]] <- do.call("rbind", cdfval_sinll[[s]])
  cdfval_sinll[[s]]$spl <- s
  cdfval_sinll[[s]]$type <- "val"
  cdfval_sinll[[s]]$loss <- "nll"
  cdfval_sinll[[s]]$mod <- "si"
}
cdfval_sinll <- do.call("rbind", cdfval_sinll)

cdftest_sinll <- list_cdfs(in_dir, fname_sinll, splits, ensembles = NULL, "test")
for (s in seq_len(splits)) {
  for (m in 1) {
    cdftest_sinll[[s]][[m]]$mem <- m
  }
  cdftest_sinll[[s]] <- do.call("rbind", cdftest_sinll[[s]])
  cdftest_sinll[[s]]$spl <- s
  cdftest_sinll[[s]]$type <- "test"
  cdftest_sinll[[s]]$loss <- "nll"
  cdftest_sinll[[s]]$mod <- "si"
}
cdftest_sinll <- do.call("rbind", cdftest_sinll)

## RPS

cdfval_sirps <- list_cdfs(in_dir, fname_sirps, splits, ensembles = NULL, "val")
for (s in seq_len(splits)) {
  for (m in 1) {
    cdfval_sirps[[s]][[m]]$mem <- m
  }
  cdfval_sirps[[s]] <- do.call("rbind", cdfval_sirps[[s]])
  cdfval_sirps[[s]]$spl <- s
  cdfval_sirps[[s]]$type <- "val"
  cdfval_sirps[[s]]$loss <- "rps"
  cdfval_sirps[[s]]$mod <- "si"
}
cdfval_sirps <- do.call("rbind", cdfval_sirps)

cdftest_sirps <- list_cdfs(in_dir, fname_sirps, splits, ensembles = NULL, "test")
for (s in seq_len(splits)) {
  for (m in 1) {
    cdftest_sirps[[s]][[m]]$mem <- m
  }
  cdftest_sirps[[s]] <- do.call("rbind", cdftest_sirps[[s]])
  cdftest_sirps[[s]]$spl <- s
  cdftest_sirps[[s]]$type <- "test"
  cdftest_sirps[[s]]$loss <- "rps"
  cdftest_sirps[[s]]$mod <- "si"
}
cdftest_sirps <- do.call("rbind", cdftest_sirps)

cdf_si <- bind_rows(cdfval_sinll, cdfval_sirps, cdftest_sinll, cdftest_sirps)
write.csv(cdf_si, paste0(out_dir, "mela_merged_cdf_si.csv"), row.names = FALSE)


### SI-LS

## NLL

cdfval_silsnll <- list_cdfs(in_dir, fname_silsnll, splits, ensembles = NULL, "val")
for (s in seq_len(splits)) {
  for (m in 1) {
    cdfval_silsnll[[s]][[m]]$mem <- m
  }
  cdfval_silsnll[[s]] <- do.call("rbind", cdfval_silsnll[[s]])
  cdfval_silsnll[[s]]$spl <- s
  cdfval_silsnll[[s]]$type <- "val"
  cdfval_silsnll[[s]]$loss <- "nll"
  cdfval_silsnll[[s]]$mod <- "sils"
}
cdfval_silsnll <- do.call("rbind", cdfval_silsnll)

cdftest_silsnll <- list_cdfs(in_dir, fname_silsnll, splits, ensembles = NULL, "test")
for (s in seq_len(splits)) {
  for (m in 1) {
    cdftest_silsnll[[s]][[m]]$mem <- m
  }
  cdftest_silsnll[[s]] <- do.call("rbind", cdftest_silsnll[[s]])
  cdftest_silsnll[[s]]$spl <- s
  cdftest_silsnll[[s]]$type <- "test"
  cdftest_silsnll[[s]]$loss <- "nll"
  cdftest_silsnll[[s]]$mod <- "sils"
}
cdftest_silsnll <- do.call("rbind", cdftest_silsnll)

## RPS

cdfval_silsrps <- list_cdfs(in_dir, fname_silsrps, splits, ensembles = NULL, "val")
for (s in seq_len(splits)) {
  for (m in 1) {
    cdfval_silsrps[[s]][[m]]$mem <- m
  }
  cdfval_silsrps[[s]] <- do.call("rbind", cdfval_silsrps[[s]])
  cdfval_silsrps[[s]]$spl <- s
  cdfval_silsrps[[s]]$type <- "val"
  cdfval_silsrps[[s]]$loss <- "rps"
  cdfval_silsrps[[s]]$mod <- "sils"
}
cdfval_silsrps <- do.call("rbind", cdfval_silsrps)

cdftest_silsrps <- list_cdfs(in_dir, fname_silsrps, splits, ensembles = NULL, "test")
for (s in seq_len(splits)) {
  for (m in 1) {
    cdftest_silsrps[[s]][[m]]$mem <- m
  }
  cdftest_silsrps[[s]] <- do.call("rbind", cdftest_silsrps[[s]])
  cdftest_silsrps[[s]]$spl <- s
  cdftest_silsrps[[s]]$type <- "test"
  cdftest_silsrps[[s]]$loss <- "rps"
  cdftest_silsrps[[s]]$mod <- "sils"
}
cdftest_silsrps <- do.call("rbind", cdftest_silsrps)

cdf_sils <- bind_rows(cdfval_silsnll, cdfval_silsrps, cdftest_silsnll, cdftest_silsrps)
write.csv(cdf_sils, paste0(out_dir, "mela_merged_cdf_sils.csv"), row.names = FALSE)

# Combine true Y ----------------------------------------------------------

y <- model.matrix(~ 0 + target, data = tab_dat)

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
write.csv(y_all, paste0(out_dir, "mela_merged_y.csv"), row.names = FALSE)
