# Bootstrap confidence intervals
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(etram)
library(boot)

# Directories -------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/UTKFace/"

# Params ------------------------------------------------------------------

K <- 7

fname_silscsnll <- "utkface_silscs_lossnll_wsyes_augno"
fname_silscsrps <- "utkface_silscs_lossrps_wsyes_augno"
fname_cilsnll <- "utkface_cils_lossnll_wsyes_augno"
fname_cilsrps <- "utkface_cils_lossrps_wsyes_augno"
fname_sicsnll <- "utkface_sics_lossnll_wsyes_augno"
fname_sicsrps <- "utkface_sics_lossrps_wsyes_augno"
fname_cinll <- "utkface_ci_lossnll_wsyes_augno"
fname_cirps <- "utkface_ci_lossrps_wsyes_augno"

fname_si <- "utkface_si"
fname_sirps <- "utkface_si_rps"
fname_sils <- "utkface_sils"
fname_silsrps <- "utkface_sils_rps"

# Load results ------------------------------------------------------------

## Load results of reference model SI-LS
met_sils <-  read.csv(file = paste0(in_dir, "met_", fname_sils, ".csv"))
met_silsrps <- read.csv(file = paste0(in_dir, "met_", fname_silsrps, ".csv"))

## all CDF
cdf_files <- list.files(path = in_dir,
                        pattern = paste0("utkface_merged_cdf.*\\.csv$"))
cdf_files <- lapply(cdf_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})
all_cdf <- do.call("rbind", cdf_files)

## all Y
all_y <- read.csv(paste0(in_dir, "utkface_merged_y.csv"))

## CDFs all splits

### SI-LS-CS
cdftest_silscsnll <- load_lys_cdf_all(all_cdf, m = "silscs", K = K, l = "nll", t = "test")
cdftest_silscsrps <- load_lys_cdf_all(all_cdf, m = "silscs", K = K, l = "rps", t = "test")

### SI-CS
cdftest_sicsnll <- load_lys_cdf_all(all_cdf, m = "sics", K = K, l = "nll", t = "test")
cdftest_sicsrps <- load_lys_cdf_all(all_cdf, m = "sics", K = K, l = "rps", t = "test")

### CI-LS
cdftest_cilsnll <- load_lys_cdf_all(all_cdf, m = "cils", K = K, l = "nll", t = "test")
cdftest_cilsrps <- load_lys_cdf_all(all_cdf, m = "cils", K = K, l = "rps", t = "test")

### CI
cdftest_cinll <- load_lys_cdf_all(all_cdf, m = "ci", K = K, l = "nll", t = "test")
cdftest_cirps <- load_lys_cdf_all(all_cdf, m = "ci", K = K, l = "rps", t = "test")

### SI
cdftest_si <- load_lys_cdf_all(all_cdf, m = "si", K = K, l = "nll", t = "test")
cdftest_sirps <- load_lys_cdf_all(all_cdf, m = "si", K = K, l = "rps", t = "test")

### SI-LS
cdftest_sils <- load_lys_cdf_all(all_cdf, m = "sils", K = K, l = "nll", t = "test")
cdftest_silsrps <- load_lys_cdf_all(all_cdf, m = "sils", K = K, l = "rps", t = "test")

## Y true

y_true_all <- load_y_true_all(all_y, K = K, t = "test")

### Weights

# SI-LS-CS
w_silscsnll <- read.csv(paste0(in_dir, "w_", fname_silscsnll, ".csv"))
w_silscsnll_l <- extract_w(w_silscsnll, meth = "linear")
w_silscsnll_ll <- extract_w(w_silscsnll, meth = "log-linear")
w_silscsnll_trf <- extract_w(w_silscsnll, meth = "trafo")

w_silscsrps <- read.csv(paste0(in_dir, "w_", fname_silscsrps, ".csv"))
w_silscsrps_l <- extract_w(w_silscsrps, meth = "linear")
w_silscsrps_ll <- extract_w(w_silscsrps, meth = "log-linear")
w_silscsrps_trf <- extract_w(w_silscsrps, meth = "trafo")

# SI-CS
w_sicsnll <- read.csv(paste0(in_dir, "w_", fname_sicsnll, ".csv"))
w_sicsnll_l <- extract_w(w_sicsnll, meth = "linear")
w_sicsnll_ll <- extract_w(w_sicsnll, meth = "log-linear")
w_sicsnll_trf <- extract_w(w_sicsnll, meth = "trafo")

w_sicsrps <- read.csv(paste0(in_dir, "w_", fname_sicsrps, ".csv"))
w_sicsrps_l <- extract_w(w_sicsrps, meth = "linear")
w_sicsrps_ll <- extract_w(w_sicsrps, meth = "log-linear")
w_sicsrps_trf <- extract_w(w_sicsrps, meth = "trafo")

# CI-LS
w_cilsnll <- read.csv(paste0(in_dir, "w_", fname_cilsnll, ".csv"))
w_cilsnll_l <- extract_w(w_cilsnll, meth = "linear")
w_cilsnll_ll <- extract_w(w_cilsnll, meth = "log-linear")
w_cilsnll_trf <- extract_w(w_cilsnll, meth = "trafo")

w_cilsrps <- read.csv(paste0(in_dir, "w_", fname_cilsrps, ".csv"))
w_cilsrps_l <- extract_w(w_cilsrps, meth = "linear")
w_cilsrps_ll <- extract_w(w_cilsrps, meth = "log-linear")
w_cilsrps_trf <- extract_w(w_cilsrps, meth = "trafo")

# CI
w_cinll <- read.csv(paste0(in_dir, "w_", fname_cinll, ".csv"))
w_cinll_l <- extract_w(w_cinll, meth = "linear")
w_cinll_ll <- extract_w(w_cinll, meth = "log-linear")
w_cinll_trf <- extract_w(w_cinll, meth = "trafo")

w_cirps <- read.csv(paste0(in_dir, "w_", fname_cirps, ".csv"))
w_cirps_l <- extract_w(w_cirps, meth = "linear")
w_cirps_ll <- extract_w(w_cirps, meth = "log-linear")
w_cirps_trf <- extract_w(w_cirps, meth = "trafo")


# Construct ensembles (equal weights) -------------------------------------

# SI-LS-CS
ens_silscsnll_l <- lapply(lapply(cdftest_silscsnll, get_ensemble, type = "linear"), list)
ens_silscsnll_ll <- lapply(lapply(cdftest_silscsnll, get_ensemble, type = "log-linear"), list)
ens_silscsnll_trf <- lapply(lapply(cdftest_silscsnll, get_ensemble, type = "trafo"), list)

ens_silscsrps_l <- lapply(lapply(cdftest_silscsrps, get_ensemble, type = "linear"), list)
ens_silscsrps_ll <- lapply(lapply(cdftest_silscsrps, get_ensemble, type = "log-linear"), list)
ens_silscsrps_trf <- lapply(lapply(cdftest_silscsrps, get_ensemble, type = "trafo"), list)

# SI-CS
ens_sicsnll_l <- lapply(lapply(cdftest_sicsnll, get_ensemble, type = "linear"), list)
ens_sicsnll_ll <- lapply(lapply(cdftest_sicsnll, get_ensemble, type = "log-linear"), list)
ens_sicsnll_trf <- lapply(lapply(cdftest_sicsnll, get_ensemble, type = "trafo"), list)

ens_sicsrps_l <- lapply(lapply(cdftest_sicsrps, get_ensemble, type = "linear"), list)
ens_sicsrps_ll <- lapply(lapply(cdftest_sicsrps, get_ensemble, type = "log-linear"), list)
ens_sicsrps_trf <- lapply(lapply(cdftest_sicsrps, get_ensemble, type = "trafo"), list)

# CI-LS
ens_cilsnll_l <- lapply(lapply(cdftest_cilsnll, get_ensemble, type = "linear"), list)
ens_cilsnll_ll <- lapply(lapply(cdftest_cilsnll, get_ensemble, type = "log-linear"), list)
ens_cilsnll_trf <- lapply(lapply(cdftest_cilsnll, get_ensemble, type = "trafo"), list)

ens_cilsrps_l <- lapply(lapply(cdftest_cilsrps, get_ensemble, type = "linear"), list)
ens_cilsrps_ll <- lapply(lapply(cdftest_cilsrps, get_ensemble, type = "log-linear"), list)
ens_cilsrps_trf <- lapply(lapply(cdftest_cilsrps, get_ensemble, type = "trafo"), list)

# CI
ens_cinll_l <- lapply(lapply(cdftest_cinll, get_ensemble, type = "linear"), list)
ens_cinll_ll <- lapply(lapply(cdftest_cinll, get_ensemble, type = "log-linear"), list)
ens_cinll_trf <- lapply(lapply(cdftest_cinll, get_ensemble, type = "trafo"), list)

ens_cirps_l <- lapply(lapply(cdftest_cirps, get_ensemble, type = "linear"), list)
ens_cirps_ll <- lapply(lapply(cdftest_cirps, get_ensemble, type = "log-linear"), list)
ens_cirps_trf <- lapply(lapply(cdftest_cirps, get_ensemble, type = "trafo"), list)

# Construct ensembles (weighted) ------------------------------------------

# SI-LS-CS
ens_w_silscsnll_l <- lapply(get_weighted_ens(lys_cdf_all = cdftest_silscsnll, 
                                             lys_weigths = w_silscsnll_l, type = "linear"), list)
ens_w_silscsnll_ll <- lapply(get_weighted_ens(lys_cdf_all = cdftest_silscsnll, 
                                              lys_weigths = w_silscsnll_ll, type = "log-linear"), list)
ens_w_silscsnll_trf <- lapply(get_weighted_ens(lys_cdf_all = cdftest_silscsnll,
                                               lys_weigths = w_silscsnll_trf, type = "trafo"), list)

ens_w_silscsrps_l <- lapply(get_weighted_ens(lys_cdf_all = cdftest_silscsrps, 
                                             lys_weigths = w_silscsrps_l, type = "linear"), list)
ens_w_silscsrps_ll <- lapply(get_weighted_ens(lys_cdf_all = cdftest_silscsrps, 
                                              lys_weigths = w_silscsrps_ll, type = "log-linear"), list)
ens_w_silscsrps_trf <- lapply(get_weighted_ens(lys_cdf_all = cdftest_silscsrps,
                                               lys_weigths = w_silscsrps_trf, type = "trafo"), list)
# SI-CS
ens_w_sicsnll_l <- lapply(get_weighted_ens(lys_cdf_all = cdftest_sicsnll, 
                                           lys_weigths = w_sicsnll_l, type = "linear"), list)
ens_w_sicsnll_ll <- lapply(get_weighted_ens(lys_cdf_all = cdftest_sicsnll, 
                                            lys_weigths = w_sicsnll_ll, type = "log-linear"), list)
ens_w_sicsnll_trf <- lapply(get_weighted_ens(lys_cdf_all = cdftest_sicsnll,
                                             lys_weigths = w_sicsnll_trf, type = "trafo"), list)

ens_w_sicsrps_l <- lapply(get_weighted_ens(lys_cdf_all = cdftest_sicsrps, 
                                           lys_weigths = w_sicsrps_l, type = "linear"), list)
ens_w_sicsrps_ll <- lapply(get_weighted_ens(lys_cdf_all = cdftest_sicsrps, 
                                            lys_weigths = w_sicsrps_ll, type = "log-linear"), list)
ens_w_sicsrps_trf <- lapply(get_weighted_ens(lys_cdf_all = cdftest_sicsrps,
                                             lys_weigths = w_sicsrps_trf, type = "trafo"), list)

# CI-LS
ens_w_cilsnll_l <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cilsnll, 
                                           lys_weigths = w_cilsnll_l, type = "linear"), list)
ens_w_cilsnll_ll <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cilsnll, 
                                            lys_weigths = w_cilsnll_ll, type = "log-linear"), list)
ens_w_cilsnll_trf <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cilsnll,
                                             lys_weigths = w_cilsnll_trf, type = "trafo"), list)

ens_w_cilsrps_l <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cilsrps, 
                                           lys_weigths = w_cilsrps_l, type = "linear"), list)
ens_w_cilsrps_ll <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cilsrps, 
                                            lys_weigths = w_cilsrps_ll, type = "log-linear"), list)
ens_w_cilsrps_trf <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cilsrps,
                                             lys_weigths = w_cilsrps_trf, type = "trafo"), list)


# CI
ens_w_cinll_l <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cinll, 
                                         lys_weigths = w_cinll_l, type = "linear"), list)
ens_w_cinll_ll <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cinll, 
                                          lys_weigths = w_cinll_ll, type = "log-linear"), list)
ens_w_cinll_trf <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cinll,
                                           lys_weigths = w_cinll_trf, type = "trafo"), list)

ens_w_cirps_l <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cirps, 
                                         lys_weigths = w_cirps_l, type = "linear"), list)
ens_w_cirps_ll <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cirps, 
                                          lys_weigths = w_cirps_ll, type = "log-linear"), list)
ens_w_cirps_trf <- lapply(get_weighted_ens(lys_cdf_all = cdftest_cirps,
                                           lys_weigths = w_cirps_trf, type = "trafo"), list)

# equal weights -----------------------------------------------------------

## NLL

# SI-LS-CS

# 'normal average'
boot_silscsnll_avg <- get_bs(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_silscsnll_avg$mod <- "silscs"
boot_silscsnll_avg$method <- "avg"
write.csv(boot_silscsnll_avg, file = paste0(out_dir, "boot_silscsnll_avg.csv"), row.names = FALSE)

# linear ensemble
boot_silscsnll_l <- get_bs(lys_cdf_all = ens_silscsnll_l, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_silscsnll_l$mod <- "silscs"
boot_silscsnll_l$method <- "linear"
write.csv(boot_silscsnll_l, file = paste0(out_dir, "boot_silscsnll_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_silscsnll_ll <- get_bs(lys_cdf_all = ens_silscsnll_ll, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_silscsnll_ll$mod <- "silscs"
boot_silscsnll_ll$method <- "log-linear"
write.csv(boot_silscsnll_ll, file = paste0(out_dir, "boot_silscsnll_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_silscsnll_trf <- get_bs(lys_cdf_all = ens_silscsnll_trf, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_silscsnll_trf$mod <- "silscs"
boot_silscsnll_trf$method <- "trafo"
write.csv(boot_silscsnll_trf, file = paste0(out_dir, "boot_silscsnll_trf.csv"), row.names = FALSE)

# SI-CS

# 'normal average'
boot_sicsnll_avg <- get_bs(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_sicsnll_avg$mod <- "sics"
boot_sicsnll_avg$method <- "avg"
write.csv(boot_sicsnll_avg, file = paste0(out_dir, "boot_sicsnll_avg.csv"), row.names = FALSE)

# linear ensemble
boot_sicsnll_l <- get_bs(lys_cdf_all = ens_sicsnll_l, y_true_all = y_true_all, R = 1000,
                         binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_sicsnll_l$mod <- "sics"
boot_sicsnll_l$method <- "linear"
write.csv(boot_sicsnll_l, file = paste0(out_dir, "boot_sicsnll_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_sicsnll_ll <- get_bs(lys_cdf_all = ens_sicsnll_ll, y_true_all = y_true_all, R = 1000,
                          binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_sicsnll_ll$mod <- "sics"
boot_sicsnll_ll$method <- "log-linear"
write.csv(boot_sicsnll_ll, file = paste0(out_dir, "boot_sicsnll_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_sicsnll_trf <- get_bs(lys_cdf_all = ens_sicsnll_trf, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_sicsnll_trf$mod <- "sics"
boot_sicsnll_trf$method <- "trafo"
write.csv(boot_sicsnll_trf, file = paste0(out_dir, "boot_sicsnll_trf.csv"), row.names = FALSE)

# CI-LS

# 'normal average'
boot_cilsnll_avg <- get_bs(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_cilsnll_avg$mod <- "cils"
boot_cilsnll_avg$method <- "avg"
write.csv(boot_cilsnll_avg, file = paste0(out_dir, "boot_cilsnll_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cilsnll_l <- get_bs(lys_cdf_all = ens_cilsnll_l, y_true_all = y_true_all, R = 1000,
                         binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_cilsnll_l$mod <- "cils"
boot_cilsnll_l$method <- "linear"
write.csv(boot_cilsnll_l, file = paste0(out_dir, "boot_cilsnll_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cilsnll_ll <- get_bs(lys_cdf_all = ens_cilsnll_ll, y_true_all = y_true_all, R = 1000,
                          binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_cilsnll_ll$mod <- "cils"
boot_cilsnll_ll$method <- "log-linear"
write.csv(boot_cilsnll_ll, file = paste0(out_dir, "boot_cilsnll_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cilsnll_trf <- get_bs(lys_cdf_all = ens_cilsnll_trf, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_cilsnll_trf$mod <- "cils"
boot_cilsnll_trf$method <- "trafo"
write.csv(boot_cilsnll_trf, file = paste0(out_dir, "boot_cilsnll_trf.csv"), row.names = FALSE)

# CI

# 'normal average'
boot_cinll_avg <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000,
                         binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_cinll_avg$mod <- "ci"
boot_cinll_avg$method <- "avg"
write.csv(boot_cinll_avg, file = paste0(out_dir, "boot_cinll_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cinll_l <- get_bs(lys_cdf_all = ens_cinll_l, y_true_all = y_true_all, R = 1000, 
                       binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_cinll_l$mod <- "ci"
boot_cinll_l$method <- "linear"
write.csv(boot_cinll_l, file = paste0(out_dir, "boot_cinll_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cinll_ll <- get_bs(lys_cdf_all = ens_cinll_ll, y_true_all = y_true_all, R = 1000, 
                        binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_cinll_ll$mod <- "ci"
boot_cinll_ll$method <- "log-linear"
write.csv(boot_cinll_ll, file = paste0(out_dir, "boot_cinll_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cinll_trf <- get_bs(lys_cdf_all = ens_cinll_trf, y_true_all = y_true_all, R = 1000,
                         binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_cinll_trf$mod <- "ci"
boot_cinll_trf$method <- "trafo"
write.csv(boot_cinll_trf, file = paste0(out_dir, "boot_cinll_trf.csv"), row.names = FALSE)

# SI-LS

boot_sils <- get_bs(lys_cdf_all = cdftest_sils, y_true_all = y_true_all, R = 1000, 
                    binary = FALSE, ncpus = 7) # no ref 
boot_sils$mod <- "sils"
boot_sils$method <- NA
write.csv(boot_sils, file = paste0(out_dir, "boot_silsnll.csv"), row.names = FALSE)

# SI

boot_si <- get_bs(lys_cdf_all = cdftest_si, y_true_all = y_true_all, R = 1000, 
                  binary = FALSE, ncpus = 7, met_ref = met_sils)
boot_si$mod <- "si"
boot_si$method <- NA
write.csv(boot_si, file = paste0(out_dir, "boot_sinll.csv"), row.names = FALSE)


## RPS

# SI-LS-CS

# 'normal average'
boot_silscsrps_avg <- get_bs(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_silscsrps_avg$mod <- "silscs"
boot_silscsrps_avg$method <- "avg"
write.csv(boot_silscsrps_avg, file = paste0(out_dir, "boot_silscsrps_avg.csv"), row.names = FALSE)

# linear ensemble
boot_silscsrps_l <- get_bs(lys_cdf_all = ens_silscsrps_l, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_silscsrps_l$mod <- "silscs"
boot_silscsrps_l$method <- "linear"
write.csv(boot_silscsrps_l, file = paste0(out_dir, "boot_silscsrps_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_silscsrps_ll <- get_bs(lys_cdf_all = ens_silscsrps_ll, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_silscsrps_ll$mod <- "silscs"
boot_silscsrps_ll$method <- "log-linear"
write.csv(boot_silscsrps_ll, file = paste0(out_dir, "boot_silscsrps_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_silscsrps_trf <- get_bs(lys_cdf_all = ens_silscsrps_trf, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_silscsrps_trf$mod <- "silscs"
boot_silscsrps_trf$method <- "trafo"
write.csv(boot_silscsrps_trf, file = paste0(out_dir, "boot_silscsrps_trf.csv"), row.names = FALSE)

# SI-CS

# 'normal average'
boot_sicsrps_avg <- get_bs(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_sicsrps_avg$mod <- "sics"
boot_sicsrps_avg$method <- "avg"
write.csv(boot_sicsrps_avg, file = paste0(out_dir, "boot_sicsrps_avg.csv"), row.names = FALSE)

# linear ensemble
boot_sicsrps_l <- get_bs(lys_cdf_all = ens_sicsrps_l, y_true_all = y_true_all, R = 1000,
                         binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_sicsrps_l$mod <- "sics"
boot_sicsrps_l$method <- "linear"
write.csv(boot_sicsrps_l, file = paste0(out_dir, "boot_sicsrps_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_sicsrps_ll <- get_bs(lys_cdf_all = ens_sicsrps_ll, y_true_all = y_true_all, R = 1000,
                          binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_sicsrps_ll$mod <- "sics"
boot_sicsrps_ll$method <- "log-linear"
write.csv(boot_sicsrps_ll, file = paste0(out_dir, "boot_sicsrps_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_sicsrps_trf <- get_bs(lys_cdf_all = ens_sicsrps_trf, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_sicsrps_trf$mod <- "sics"
boot_sicsrps_trf$method <- "trafo"
write.csv(boot_sicsrps_trf, file = paste0(out_dir, "boot_sicsrps_trf.csv"), row.names = FALSE)

# CI-LS

# 'normal average'
boot_cilsrps_avg <- get_bs(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_cilsrps_avg$mod <- "cils"
boot_cilsrps_avg$method <- "avg"
write.csv(boot_cilsrps_avg, file = paste0(out_dir, "boot_cilsrps_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cilsrps_l <- get_bs(lys_cdf_all = ens_cilsrps_l, y_true_all = y_true_all, R = 1000,
                         binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_cilsrps_l$mod <- "cils"
boot_cilsrps_l$method <- "linear"
write.csv(boot_cilsrps_l, file = paste0(out_dir, "boot_cilsrps_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cilsrps_ll <- get_bs(lys_cdf_all = ens_cilsrps_ll, y_true_all = y_true_all, R = 1000,
                          binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_cilsrps_ll$mod <- "cils"
boot_cilsrps_ll$method <- "log-linear"
write.csv(boot_cilsrps_ll, file = paste0(out_dir, "boot_cilsrps_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cilsrps_trf <- get_bs(lys_cdf_all = ens_cilsrps_trf, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_cilsrps_trf$mod <- "cils"
boot_cilsrps_trf$method <- "trafo"
write.csv(boot_cilsrps_trf, file = paste0(out_dir, "boot_cilsrps_trf.csv"), row.names = FALSE)

# CI

# 'normal average'
boot_cirps_avg <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000,
                         binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_cirps_avg$mod <- "ci"
boot_cirps_avg$method <- "avg"
write.csv(boot_cirps_avg, file = paste0(out_dir, "boot_cirps_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cirps_l <- get_bs(lys_cdf_all = ens_cirps_l, y_true_all = y_true_all, R = 1000, 
                       binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_cirps_l$mod <- "ci"
boot_cirps_l$method <- "linear"
write.csv(boot_cirps_l, file = paste0(out_dir, "boot_cirps_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cirps_ll <- get_bs(lys_cdf_all = ens_cirps_ll, y_true_all = y_true_all, R = 1000, 
                        binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_cirps_ll$mod <- "ci"
boot_cirps_ll$method <- "log-linear"
write.csv(boot_cirps_ll, file = paste0(out_dir, "boot_cirps_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cirps_trf <- get_bs(lys_cdf_all = ens_cirps_trf, y_true_all = y_true_all, R = 1000,
                         binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_cirps_trf$mod <- "ci"
boot_cirps_trf$method <- "trafo"
write.csv(boot_cirps_trf, file = paste0(out_dir, "boot_cirps_trf.csv"), row.names = FALSE)

# SI-LS

boot_sils <- get_bs(lys_cdf_all = cdftest_sils, y_true_all = y_true_all, R = 1000, 
                    binary = FALSE, ncpus = 7) # no ref 
boot_sils$mod <- "sils"
boot_sils$method <- NA
write.csv(boot_sils, file = paste0(out_dir, "boot_silsrps.csv"), row.names = FALSE)

# SI

boot_si <- get_bs(lys_cdf_all = cdftest_si, y_true_all = y_true_all, R = 1000, 
                  binary = FALSE, ncpus = 7, met_ref = met_silsrps)
boot_si$mod <- "si"
boot_si$method <- NA
write.csv(boot_si, file = paste0(out_dir, "boot_sirps.csv"), row.names = FALSE)


# weighted ----------------------------------------------------------------

## NLL

# SI-LS-CS

# linear average
wboot_silscsnll_avgl <- get_bs(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, R = 1000,
                               binary = FALSE, ncpus = 7,
                               weights = w_silscsnll_l, met_ref = met_sils)
wboot_silscsnll_avgl$mod <- "silscs"
wboot_silscsnll_avgl$method <- "avg"
write.csv(wboot_silscsnll_avgl, file = paste0(out_dir, "wboot_silscsnll_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_silscsnll_avgll <- get_bs(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, R = 1000,
                                binary = FALSE, ncpus = 7,
                                weights = w_silscsnll_ll, met_ref = met_sils)
wboot_silscsnll_avgll$mod <- "silscs"
wboot_silscsnll_avgll$method <- "avgll"
write.csv(wboot_silscsnll_avgll, file = paste0(out_dir, "wboot_silscsnll_avgll.csv"), row.names = FALSE)

# trafo average
wboot_silscsnll_avgtrf <- get_bs(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, R = 1000,
                                 binary = FALSE, ncpus = 7,
                                 weights = w_silscsnll_trf, met_ref = met_sils)
wboot_silscsnll_avgtrf$mod <- "silscs"
wboot_silscsnll_avgtrf$method <- "avgtrf"
write.csv(wboot_silscsnll_avgtrf, file = paste0(out_dir, "wboot_silscsnll_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_silscsnll_l <- get_bs(lys_cdf_all = ens_w_silscsnll_l, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_silscsnll_l$mod <- "silscs"
wboot_silscsnll_l$method <- "linear"
write.csv(wboot_silscsnll_l, file = paste0(out_dir, "wboot_silscsnll_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_silscsnll_ll <- get_bs(lys_cdf_all = ens_w_silscsnll_ll, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_silscsnll_ll$mod <- "silscs"
wboot_silscsnll_ll$method <- "log-linear"
write.csv(wboot_silscsnll_ll, file = paste0(out_dir, "wboot_silscsnll_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_silscsnll_trf <- get_bs(lys_cdf_all = ens_w_silscsnll_trf, y_true_all = y_true_all, R = 1000,
                              binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_silscsnll_trf$mod <- "silscs"
wboot_silscsnll_trf$method <- "trafo"
write.csv(wboot_silscsnll_trf, file = paste0(out_dir, "wboot_silscsnll_trf.csv"), row.names = FALSE)

# SI-CS

# linear average
wboot_sicsnll_avgl <- get_bs(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7,
                             weights = w_sicsnll_l, met_ref = met_sils)
wboot_sicsnll_avgl$mod <- "sics"
wboot_sicsnll_avgl$method <- "avg"
write.csv(wboot_sicsnll_avgl, file = paste0(out_dir, "wboot_sicsnll_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_sicsnll_avgll <- get_bs(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, R = 1000,
                              binary = FALSE, ncpus = 7,
                              weights = w_sicsnll_ll, met_ref = met_sils)
wboot_sicsnll_avgll$mod <- "sics"
wboot_sicsnll_avgll$method <- "avgll"
write.csv(wboot_sicsnll_avgll, file = paste0(out_dir, "wboot_sicsnll_avgll.csv"), row.names = FALSE)

# trafo average
wboot_sicsnll_avgtrf <- get_bs(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, R = 1000,
                               binary = FALSE, ncpus = 7,
                               weights = w_sicsnll_trf, met_ref = met_sils)
wboot_sicsnll_avgtrf$mod <- "sics"
wboot_sicsnll_avgtrf$method <- "avgtrf"
write.csv(wboot_sicsnll_avgtrf, file = paste0(out_dir, "wboot_sicsnll_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_sicsnll_l <- get_bs(lys_cdf_all = ens_w_sicsnll_l, y_true_all = y_true_all, R = 1000,
                          binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_sicsnll_l$mod <- "sics"
wboot_sicsnll_l$method <- "linear"
write.csv(wboot_sicsnll_l, file = paste0(out_dir, "wboot_sicsnll_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_sicsnll_ll <- get_bs(lys_cdf_all = ens_w_sicsnll_ll, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_sicsnll_ll$mod <- "sics"
wboot_sicsnll_ll$method <- "log-linear"
write.csv(wboot_sicsnll_ll, file = paste0(out_dir, "wboot_sicsnll_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_sicsnll_trf <- get_bs(lys_cdf_all = ens_w_sicsnll_trf, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_sicsnll_trf$mod <- "sics"
wboot_sicsnll_trf$method <- "trafo"
write.csv(wboot_sicsnll_trf, file = paste0(out_dir, "wboot_sicsnll_trf.csv"), row.names = FALSE)

# CI-LS

# linear average
wboot_cilsnll_avgl <- get_bs(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7,
                             weights = w_cilsnll_l, met_ref = met_sils)
wboot_cilsnll_avgl$mod <- "cils"
wboot_cilsnll_avgl$method <- "avg"
write.csv(wboot_cilsnll_avgl, file = paste0(out_dir, "wboot_cilsnll_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cilsnll_avgll <- get_bs(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, R = 1000,
                              binary = FALSE, ncpus = 7,
                              weights = w_cilsnll_ll, met_ref = met_sils)
wboot_cilsnll_avgll$mod <- "cils"
wboot_cilsnll_avgll$method <- "avgll"
write.csv(wboot_cilsnll_avgll, file = paste0(out_dir, "wboot_cilsnll_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cilsnll_avgtrf <- get_bs(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, R = 1000,
                               binary = FALSE, ncpus = 7,
                               weights = w_cilsnll_trf, met_ref = met_sils)
wboot_cilsnll_avgtrf$mod <- "cils"
wboot_cilsnll_avgtrf$method <- "avgtrf"
write.csv(wboot_cilsnll_avgtrf, file = paste0(out_dir, "wboot_cilsnll_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cilsnll_l <- get_bs(lys_cdf_all = ens_w_cilsnll_l, y_true_all = y_true_all, R = 1000,
                          binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_cilsnll_l$mod <- "cils"
wboot_cilsnll_l$method <- "linear"
write.csv(wboot_cilsnll_l, file = paste0(out_dir, "wboot_cilsnll_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cilsnll_ll <- get_bs(lys_cdf_all = ens_w_cilsnll_ll, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_cilsnll_ll$mod <- "cils"
wboot_cilsnll_ll$method <- "log-linear"
write.csv(wboot_cilsnll_ll, file = paste0(out_dir, "wboot_cilsnll_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cilsnll_trf <- get_bs(lys_cdf_all = ens_w_cilsnll_trf, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_cilsnll_trf$mod <- "cils"
wboot_cilsnll_trf$method <- "trafo"
write.csv(wboot_cilsnll_trf, file = paste0(out_dir, "wboot_cilsnll_trf.csv"), row.names = FALSE)

# CI

# linear average
wboot_cinll_avgl <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7,
                           weights = w_cinll_l, met_ref = met_sils)
wboot_cinll_avgl$mod <- "ci"
wboot_cinll_avgl$method <- "avg"
write.csv(wboot_cinll_avgl, file = paste0(out_dir, "wboot_cinll_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cinll_avgll <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7,
                            weights = w_cinll_ll, met_ref = met_sils)
wboot_cinll_avgll$mod <- "ci"
wboot_cinll_avgll$method <- "avgll"
write.csv(wboot_cinll_avgll, file = paste0(out_dir, "wboot_cinll_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cinll_avgtrf <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000, 
                             binary = FALSE, ncpus = 7,
                             weights = w_cinll_trf, met_ref = met_sils)
wboot_cinll_avgtrf$mod <- "ci"
wboot_cinll_avgtrf$method <- "avgtrf"
write.csv(wboot_cinll_avgtrf, file = paste0(out_dir, "wboot_cinll_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cinll_l <- get_bs(lys_cdf_all = ens_w_cinll_l, y_true_all = y_true_all, R = 1000, 
                        binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_cinll_l$mod <- "ci"
wboot_cinll_l$method <- "linear"
write.csv(wboot_cinll_l, file = paste0(out_dir, "wboot_cinll_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cinll_ll <- get_bs(lys_cdf_all = ens_w_cinll_ll, y_true_all = y_true_all, R = 1000, 
                         binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_cinll_ll$mod <- "ci"
wboot_cinll_ll$method <- "log-linear"
write.csv(wboot_cinll_ll, file = paste0(out_dir, "wboot_cinll_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cinll_trf <- get_bs(lys_cdf_all = ens_w_cinll_trf, y_true_all = y_true_all, R = 1000, 
                          binary = FALSE, ncpus = 7, met_ref = met_sils)
wboot_cinll_trf$mod <- "ci"
wboot_cinll_trf$method <- "trafo"
write.csv(wboot_cinll_trf, file = paste0(out_dir, "wboot_cinll_trf.csv"), row.names = FALSE)


## RPS

# SI-LS-CS

# linear average
wboot_silscsrps_avgl <- get_bs(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, R = 1000,
                               binary = FALSE, ncpus = 7,
                               weights = w_silscsrps_l, met_ref = met_silsrps)
wboot_silscsrps_avgl$mod <- "silscs"
wboot_silscsrps_avgl$method <- "avg"
write.csv(wboot_silscsrps_avgl, file = paste0(out_dir, "wboot_silscsrps_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_silscsrps_avgll <- get_bs(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, R = 1000,
                                binary = FALSE, ncpus = 7,
                                weights = w_silscsrps_ll, met_ref = met_silsrps)
wboot_silscsrps_avgll$mod <- "silscs"
wboot_silscsrps_avgll$method <- "avgll"
write.csv(wboot_silscsrps_avgll, file = paste0(out_dir, "wboot_silscsrps_avgll.csv"), row.names = FALSE)

# trafo average
wboot_silscsrps_avgtrf <- get_bs(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, R = 1000,
                                 binary = FALSE, ncpus = 7,
                                 weights = w_silscsrps_trf, met_ref = met_silsrps)
wboot_silscsrps_avgtrf$mod <- "silscs"
wboot_silscsrps_avgtrf$method <- "avgtrf"
write.csv(wboot_silscsrps_avgtrf, file = paste0(out_dir, "wboot_silscsrps_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_silscsrps_l <- get_bs(lys_cdf_all = ens_w_silscsrps_l, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_silscsrps_l$mod <- "silscs"
wboot_silscsrps_l$method <- "linear"
write.csv(wboot_silscsrps_l, file = paste0(out_dir, "wboot_silscsrps_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_silscsrps_ll <- get_bs(lys_cdf_all = ens_w_silscsrps_ll, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_silscsrps_ll$mod <- "silscs"
wboot_silscsrps_ll$method <- "log-linear"
write.csv(wboot_silscsrps_ll, file = paste0(out_dir, "wboot_silscsrps_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_silscsrps_trf <- get_bs(lys_cdf_all = ens_w_silscsrps_trf, y_true_all = y_true_all, R = 1000,
                              binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_silscsrps_trf$mod <- "silscs"
wboot_silscsrps_trf$method <- "trafo"
write.csv(wboot_silscsrps_trf, file = paste0(out_dir, "wboot_silscsrps_trf.csv"), row.names = FALSE)

# SI-CS

# linear average
wboot_sicsrps_avgl <- get_bs(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7,
                             weights = w_sicsrps_l, met_ref = met_silsrps)
wboot_sicsrps_avgl$mod <- "sics"
wboot_sicsrps_avgl$method <- "avg"
write.csv(wboot_sicsrps_avgl, file = paste0(out_dir, "wboot_sicsrps_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_sicsrps_avgll <- get_bs(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, R = 1000,
                              binary = FALSE, ncpus = 7,
                              weights = w_sicsrps_ll, met_ref = met_silsrps)
wboot_sicsrps_avgll$mod <- "sics"
wboot_sicsrps_avgll$method <- "avgll"
write.csv(wboot_sicsrps_avgll, file = paste0(out_dir, "wboot_sicsrps_avgll.csv"), row.names = FALSE)

# trafo average
wboot_sicsrps_avgtrf <- get_bs(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, R = 1000,
                               binary = FALSE, ncpus = 7,
                               weights = w_sicsrps_trf, met_ref = met_silsrps)
wboot_sicsrps_avgtrf$mod <- "sics"
wboot_sicsrps_avgtrf$method <- "avgtrf"
write.csv(wboot_sicsrps_avgtrf, file = paste0(out_dir, "wboot_sicsrps_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_sicsrps_l <- get_bs(lys_cdf_all = ens_w_sicsrps_l, y_true_all = y_true_all, R = 1000,
                          binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_sicsrps_l$mod <- "sics"
wboot_sicsrps_l$method <- "linear"
write.csv(wboot_sicsrps_l, file = paste0(out_dir, "wboot_sicsrps_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_sicsrps_ll <- get_bs(lys_cdf_all = ens_w_sicsrps_ll, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_sicsrps_ll$mod <- "sics"
wboot_sicsrps_ll$method <- "log-linear"
write.csv(wboot_sicsrps_ll, file = paste0(out_dir, "wboot_sicsrps_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_sicsrps_trf <- get_bs(lys_cdf_all = ens_w_sicsrps_trf, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_sicsrps_trf$mod <- "sics"
wboot_sicsrps_trf$method <- "trafo"
write.csv(wboot_sicsrps_trf, file = paste0(out_dir, "wboot_sicsrps_trf.csv"), row.names = FALSE)

# CI-LS

# linear average
wboot_cilsrps_avgl <- get_bs(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, R = 1000,
                             binary = FALSE, ncpus = 7,
                             weights = w_cilsrps_l, met_ref = met_silsrps)
wboot_cilsrps_avgl$mod <- "cils"
wboot_cilsrps_avgl$method <- "avg"
write.csv(wboot_cilsrps_avgl, file = paste0(out_dir, "wboot_cilsrps_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cilsrps_avgll <- get_bs(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, R = 1000,
                              binary = FALSE, ncpus = 7,
                              weights = w_cilsrps_ll, met_ref = met_silsrps)
wboot_cilsrps_avgll$mod <- "cils"
wboot_cilsrps_avgll$method <- "avgll"
write.csv(wboot_cilsrps_avgll, file = paste0(out_dir, "wboot_cilsrps_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cilsrps_avgtrf <- get_bs(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, R = 1000,
                               binary = FALSE, ncpus = 7,
                               weights = w_cilsrps_trf, met_ref = met_silsrps)
wboot_cilsrps_avgtrf$mod <- "cils"
wboot_cilsrps_avgtrf$method <- "avgtrf"
write.csv(wboot_cilsrps_avgtrf, file = paste0(out_dir, "wboot_cilsrps_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cilsrps_l <- get_bs(lys_cdf_all = ens_w_cilsrps_l, y_true_all = y_true_all, R = 1000,
                          binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_cilsrps_l$mod <- "cils"
wboot_cilsrps_l$method <- "linear"
write.csv(wboot_cilsrps_l, file = paste0(out_dir, "wboot_cilsrps_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cilsrps_ll <- get_bs(lys_cdf_all = ens_w_cilsrps_ll, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_cilsrps_ll$mod <- "cils"
wboot_cilsrps_ll$method <- "log-linear"
write.csv(wboot_cilsrps_ll, file = paste0(out_dir, "wboot_cilsrps_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cilsrps_trf <- get_bs(lys_cdf_all = ens_w_cilsrps_trf, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_cilsrps_trf$mod <- "cils"
wboot_cilsrps_trf$method <- "trafo"
write.csv(wboot_cilsrps_trf, file = paste0(out_dir, "wboot_cilsrps_trf.csv"), row.names = FALSE)

# CI

# linear average
wboot_cirps_avgl <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000,
                           binary = FALSE, ncpus = 7,
                           weights = w_cirps_l, met_ref = met_silsrps)
wboot_cirps_avgl$mod <- "ci"
wboot_cirps_avgl$method <- "avg"
write.csv(wboot_cirps_avgl, file = paste0(out_dir, "wboot_cirps_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cirps_avgll <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000,
                            binary = FALSE, ncpus = 7,
                            weights = w_cirps_ll, met_ref = met_silsrps)
wboot_cirps_avgll$mod <- "ci"
wboot_cirps_avgll$method <- "avgll"
write.csv(wboot_cirps_avgll, file = paste0(out_dir, "wboot_cirps_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cirps_avgtrf <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, 
                             binary = FALSE, ncpus = 7,
                             weights = w_cirps_trf, met_ref = met_silsrps)
wboot_cirps_avgtrf$mod <- "ci"
wboot_cirps_avgtrf$method <- "avgtrf"
write.csv(wboot_cirps_avgtrf, file = paste0(out_dir, "wboot_cirps_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cirps_l <- get_bs(lys_cdf_all = ens_w_cirps_l, y_true_all = y_true_all, R = 1000, 
                        binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_cirps_l$mod <- "ci"
wboot_cirps_l$method <- "linear"
write.csv(wboot_cirps_l, file = paste0(out_dir, "wboot_cirps_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cirps_ll <- get_bs(lys_cdf_all = ens_w_cirps_ll, y_true_all = y_true_all, R = 1000, 
                         binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_cirps_ll$mod <- "ci"
wboot_cirps_ll$method <- "log-linear"
write.csv(wboot_cirps_ll, file = paste0(out_dir, "wboot_cirps_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cirps_trf <- get_bs(lys_cdf_all = ens_w_cirps_trf, y_true_all = y_true_all, R = 1000, 
                          binary = FALSE, ncpus = 7, met_ref = met_silsrps)
wboot_cirps_trf$mod <- "ci"
wboot_cirps_trf$method <- "trafo"
write.csv(wboot_cirps_trf, file = paste0(out_dir, "wboot_cirps_trf.csv"), row.names = FALSE)


