# Merge all val, test CDF of all models
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(readr)
library(rhdf5)
library(etram)

# Directories -------------------------------------------------------------

path <- "~/../data/UTKFace/UTKFace.h5"
in_dir <- out_dir <- "experiments/results/DE/UTKFace/"

# Params ------------------------------------------------------------------

splits <- 6
ensembles <- 5

fname_silscsnll <- "utkface_silscs_lossnll_wsyes_augno"
fname_silscsrps <- "utkface_silscs_lossrps_wsyes_augno"
fname_cilsnll <- "utkface_cils_lossnll_wsyes_augno"
fname_cilsrps <- "utkface_cils_lossrps_wsyes_augno"
fname_sicsnll <- "utkface_sics_lossnll_wsyes_augno"
fname_sicsrps <- "utkface_sics_lossrps_wsyes_augno"
fname_cinll <- "utkface_ci_lossnll_wsyes_augno"
fname_cirps <- "utkface_ci_lossrps_wsyes_augno"

fname_sinll <- "utkface_si"
fname_sirps <- "utkface_si_rps"
fname_silsnll <- "utkface_sils"
fname_silsrps <- "utkface_sils_rps"

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

dat <- load_data("utkface", path = path)
tab_dat <- dat$tab_dat

ridx <- get_ridx(in_dir, fname = "utkface")

# Combine single CDFs -----------------------------------------------------

### SI-LS-CS

## NLL

cdfval_silscsnll <- list_cdfs(in_dir, fname_silscsnll, splits, ensembles, "val")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdfval_silscsnll[[s]][[m]]$mem <- m
  }
  cdfval_silscsnll[[s]] <- do.call("rbind", cdfval_silscsnll[[s]])
  cdfval_silscsnll[[s]]$spl <- s
  cdfval_silscsnll[[s]]$type <- "val"
  cdfval_silscsnll[[s]]$loss <- "nll"
  cdfval_silscsnll[[s]]$mod <- "silscs"
}
cdfval_silscsnll <- do.call("rbind", cdfval_silscsnll)

cdftest_silscsnll <- list_cdfs(in_dir, fname_silscsnll, splits, ensembles, "test")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdftest_silscsnll[[s]][[m]]$mem <- m
  }
  cdftest_silscsnll[[s]] <- do.call("rbind", cdftest_silscsnll[[s]])
  cdftest_silscsnll[[s]]$spl <- s
  cdftest_silscsnll[[s]]$type <- "test"
  cdftest_silscsnll[[s]]$loss <- "nll"
  cdftest_silscsnll[[s]]$mod <- "silscs"
}
cdftest_silscsnll <- do.call("rbind", cdftest_silscsnll)

## RPS

cdfval_silscsrps <- list_cdfs(in_dir, fname_silscsrps, splits, ensembles, "val")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdfval_silscsrps[[s]][[m]]$mem <- m
  }
  cdfval_silscsrps[[s]] <- do.call("rbind", cdfval_silscsrps[[s]])
  cdfval_silscsrps[[s]]$spl <- s
  cdfval_silscsrps[[s]]$type <- "val"
  cdfval_silscsrps[[s]]$loss <- "rps"
  cdfval_silscsrps[[s]]$mod <- "silscs"
}
cdfval_silscsrps <- do.call("rbind", cdfval_silscsrps)

cdftest_silscsrps <- list_cdfs(in_dir, fname_silscsrps, splits, ensembles, "test")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdftest_silscsrps[[s]][[m]]$mem <- m
  }
  cdftest_silscsrps[[s]] <- do.call("rbind", cdftest_silscsrps[[s]])
  cdftest_silscsrps[[s]]$spl <- s
  cdftest_silscsrps[[s]]$type <- "test"
  cdftest_silscsrps[[s]]$loss <- "rps"
  cdftest_silscsrps[[s]]$mod <- "silscs"
}
cdftest_silscsrps <- do.call("rbind", cdftest_silscsrps)

cdf_silscs <- bind_rows(cdfval_silscsnll, cdfval_silscsrps, cdftest_silscsnll, cdftest_silscsrps)
write.csv(cdf_silscs, paste0(out_dir, "utkface_merged_cdf_silscs.csv"), row.names = FALSE)


### SI-CS

## NLL

cdfval_sicsnll <- list_cdfs(in_dir, fname_sicsnll, splits, ensembles, "val")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdfval_sicsnll[[s]][[m]]$mem <- m
  }
  cdfval_sicsnll[[s]] <- do.call("rbind", cdfval_sicsnll[[s]])
  cdfval_sicsnll[[s]]$spl <- s
  cdfval_sicsnll[[s]]$type <- "val"
  cdfval_sicsnll[[s]]$loss <- "nll"
  cdfval_sicsnll[[s]]$mod <- "sics"
}
cdfval_sicsnll <- do.call("rbind", cdfval_sicsnll)

cdftest_sicsnll <- list_cdfs(in_dir, fname_sicsnll, splits, ensembles, "test")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdftest_sicsnll[[s]][[m]]$mem <- m
  }
  cdftest_sicsnll[[s]] <- do.call("rbind", cdftest_sicsnll[[s]])
  cdftest_sicsnll[[s]]$spl <- s
  cdftest_sicsnll[[s]]$type <- "test"
  cdftest_sicsnll[[s]]$loss <- "nll"
  cdftest_sicsnll[[s]]$mod <- "sics"
}
cdftest_sicsnll <- do.call("rbind", cdftest_sicsnll)

## RPS

cdfval_sicsrps <- list_cdfs(in_dir, fname_sicsrps, splits, ensembles, "val")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdfval_sicsrps[[s]][[m]]$mem <- m
  }
  cdfval_sicsrps[[s]] <- do.call("rbind", cdfval_sicsrps[[s]])
  cdfval_sicsrps[[s]]$spl <- s
  cdfval_sicsrps[[s]]$type <- "val"
  cdfval_sicsrps[[s]]$loss <- "rps"
  cdfval_sicsrps[[s]]$mod <- "sics"
}
cdfval_sicsrps <- do.call("rbind", cdfval_sicsrps)

cdftest_sicsrps <- list_cdfs(in_dir, fname_sicsrps, splits, ensembles, "test")
for (s in seq_len(splits)) {
  for (m in seq_len(ensembles)) {
    cdftest_sicsrps[[s]][[m]]$mem <- m
  }
  cdftest_sicsrps[[s]] <- do.call("rbind", cdftest_sicsrps[[s]])
  cdftest_sicsrps[[s]]$spl <- s
  cdftest_sicsrps[[s]]$type <- "test"
  cdftest_sicsrps[[s]]$loss <- "rps"
  cdftest_sicsrps[[s]]$mod <- "sics"
}
cdftest_sicsrps <- do.call("rbind", cdftest_sicsrps)

cdf_sics <- bind_rows(cdfval_sicsnll, cdfval_sicsrps, cdftest_sicsnll, cdftest_sicsrps)
write.csv(cdf_sics, paste0(out_dir, "utkface_merged_cdf_sics.csv"), row.names = FALSE)


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
write.csv(cdf_cils, paste0(out_dir, "utkface_merged_cdf_cils.csv"), row.names = FALSE)


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
write.csv(cdf_ci, paste0(out_dir, "utkface_merged_cdf_ci.csv"), row.names = FALSE)


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
write.csv(cdf_si, paste0(out_dir, "utkface_merged_cdf_si.csv"), row.names = FALSE)


### SI-LS

#! Has different colnames than other data sets

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
colnames(cdfval_silsnll) <- colnames(cdftest_sirps)

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
colnames(cdftest_silsnll) <- colnames(cdftest_sirps)

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
colnames(cdfval_silsrps) <- colnames(cdftest_sirps)

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
colnames(cdftest_silsrps) <- colnames(cdftest_sirps)

cdf_sils <- bind_rows(cdfval_silsnll, cdfval_silsrps, cdftest_silsnll, cdftest_silsrps)
write.csv(cdf_sils, paste0(out_dir, "utkface_merged_cdf_sils.csv"), row.names = FALSE)

# Combine true Y ----------------------------------------------------------

y <- model.matrix(~ 0 + age_group, data = tab_dat)

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
write.csv(y_all, paste0(out_dir, "utkface_merged_y.csv"), row.names = FALSE)
