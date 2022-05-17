# Bootstrap confidence intervals
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(etram)
library(boot)

# Directories -------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/melanoma/"

# Params ------------------------------------------------------------------

K <- 2

fname_cilsnll <- "mela_silscs_lossnll_wsyes_augno"
fname_cilsrps <- "mela_silscs_lossrps_wsyes_augno"
fname_cinll <- "mela_ci_lossnll_wsyes_augno"
fname_cirps <- "mela_ci_lossrps_wsyes_augno"

fname_si <- "mela_si"
fname_sirps <- "mela_si_rps"
fname_sils <- "mela_sils"
fname_silsrps <- "mela_sils_rps"

# Load results ------------------------------------------------------------

## Load results of reference model SI-LS
met_sils <-  read.csv(file = paste0(in_dir, "met_", fname_sils, ".csv"))
met_silsrps <- read.csv(file = paste0(in_dir, "met_", fname_silsrps, ".csv"))

## all CDF
cdf_files <- list.files(path = in_dir,
                        pattern = paste0("mela_merged_cdf.*\\.csv$"))
cdf_files <- lapply(cdf_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})
all_cdf <- do.call("rbind", cdf_files)

## all Y
all_y <- read.csv(paste0(in_dir, "mela_merged_y.csv"))

## CDFs all splits

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

# CI-LS

# 'normal average'
boot_cilsnll_avg <- get_bs(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, R = 1000,
                           binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_cilsnll_avg$mod <- "cils"
boot_cilsnll_avg$method <- "avg"
write.csv(boot_cilsnll_avg, file = paste0(out_dir, "boot_cilsnll_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cilsnll_l <- get_bs(lys_cdf_all = ens_cilsnll_l, y_true_all = y_true_all, R = 1000, 
                         binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_cilsnll_l$mod <- "cils"
boot_cilsnll_l$method <- "linear"
write.csv(boot_cilsnll_l, file = paste0(out_dir, "boot_cilsnll_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cilsnll_ll <- get_bs(lys_cdf_all = ens_cilsnll_ll, y_true_all = y_true_all, R = 1000, 
                          binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_cilsnll_ll$mod <- "cils"
boot_cilsnll_ll$method <- "log-linear"
write.csv(boot_cilsnll_ll, file = paste0(out_dir, "boot_cilsnll_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cilsnll_trf <- get_bs(lys_cdf_all = ens_cilsnll_trf, y_true_all = y_true_all, R = 1000, 
                           binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_cilsnll_trf$mod <- "cils"
boot_cilsnll_trf$method <- "trafo"
write.csv(boot_cilsnll_trf, file = paste0(out_dir, "boot_cilsnll_trf.csv"), row.names = FALSE)

# CI

# 'normal average'
boot_cinll_avg <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000,
                         binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_cinll_avg$mod <- "ci"
boot_cinll_avg$method <- "avg"
write.csv(boot_cinll_avg, file = paste0(out_dir, "boot_cinll_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cinll_l <- get_bs(lys_cdf_all = ens_cinll_l, y_true_all = y_true_all, R = 1000, 
                       binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_cinll_l$mod <- "ci"
boot_cinll_l$method <- "linear"
write.csv(boot_cinll_l, file = paste0(out_dir, "boot_cinll_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cinll_ll <- get_bs(lys_cdf_all = ens_cinll_ll, y_true_all = y_true_all, R = 1000, 
                        binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_cinll_ll$mod <- "ci"
boot_cinll_ll$method <- "log-linear"
write.csv(boot_cinll_ll, file = paste0(out_dir, "boot_cinll_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cinll_trf <- get_bs(lys_cdf_all = ens_cinll_trf, y_true_all = y_true_all, R = 1000,
                         binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_cinll_trf$mod <- "ci"
boot_cinll_trf$method <- "trafo"
write.csv(boot_cinll_trf, file = paste0(out_dir, "boot_cinll_trf.csv"), row.names = FALSE)

# SI-LS

boot_sils <- get_bs(lys_cdf_all = cdftest_sils, y_true_all = y_true_all, R = 1000, 
                    binary = TRUE, ncpus = 7) # no ref 
boot_sils$mod <- "sils"
boot_sils$method <- NA
write.csv(boot_sils, file = paste0(out_dir, "boot_silsnll.csv"), row.names = FALSE)

# SI

boot_si <- get_bs(lys_cdf_all = cdftest_si, y_true_all = y_true_all, R = 1000, 
                  binary = TRUE, ncpus = 7, met_ref = met_sils)
boot_si$mod <- "si"
boot_si$method <- NA
write.csv(boot_si, file = paste0(out_dir, "boot_sinll.csv"), row.names = FALSE)


## RPS

# CI-LS

# 'normal average'
boot_cilsrps_avg <- get_bs(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, R = 1000,
                           binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_cilsrps_avg$mod <- "cils"
boot_cilsrps_avg$method <- "avg"
write.csv(boot_cilsrps_avg, file = paste0(out_dir, "boot_cilsrps_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cilsrps_l <- get_bs(lys_cdf_all = ens_cilsrps_l, y_true_all = y_true_all, R = 1000, 
                         binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_cilsrps_l$mod <- "cils"
boot_cilsrps_l$method <- "linear"
write.csv(boot_cilsrps_l, file = paste0(out_dir, "boot_cilsrps_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cilsrps_ll <- get_bs(lys_cdf_all = ens_cilsrps_ll, y_true_all = y_true_all, R = 1000, 
                          binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_cilsrps_ll$mod <- "cils"
boot_cilsrps_ll$method <- "log-linear"
write.csv(boot_cilsrps_ll, file = paste0(out_dir, "boot_cilsrps_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cilsrps_trf <- get_bs(lys_cdf_all = ens_cilsrps_trf, y_true_all = y_true_all, R = 1000, 
                           binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_cilsrps_trf$mod <- "cils"
boot_cilsrps_trf$method <- "trafo"
write.csv(boot_cilsrps_trf, file = paste0(out_dir, "boot_cilsrps_trf.csv"), row.names = FALSE)

# CI

# 'normal average'
boot_cirps_avg <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, 
                         binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_cirps_avg$mod <- "ci"
boot_cirps_avg$method <- "avg"
write.csv(boot_cirps_avg, file = paste0(out_dir, "boot_cirps_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cirps_l <- get_bs(lys_cdf_all = ens_cirps_l, y_true_all = y_true_all, R = 1000, 
                       binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_cirps_l$mod <- "ci"
boot_cirps_l$method <- "linear"
write.csv(boot_cirps_l, file = paste0(out_dir, "boot_cirps_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cirps_ll <- get_bs(lys_cdf_all = ens_cirps_ll, y_true_all = y_true_all, R = 1000, 
                        binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_cirps_ll$mod <- "ci"
boot_cirps_ll$method <- "log-linear"
write.csv(boot_cirps_ll, file = paste0(out_dir, "boot_cirps_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cirps_trf <- get_bs(lys_cdf_all = ens_cirps_trf, y_true_all = y_true_all, R = 1000, 
                         binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_cirps_trf$mod <- "ci"
boot_cirps_trf$method <- "trafo"
write.csv(boot_cirps_trf, file = paste0(out_dir, "boot_cirps_trf.csv"), row.names = FALSE)

# SI-LS

boot_silsrps <- get_bs(lys_cdf_all = cdftest_silsrps, y_true_all = y_true_all, R = 1000, 
                       binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_silsrps$mod <- "sils"
boot_silsrps$method <- NA
write.csv(boot_silsrps, file = paste0(out_dir, "boot_silsrps.csv"), row.names = FALSE)

# SI

boot_sirps <- get_bs(lys_cdf_all = cdftest_sirps, y_true_all = y_true_all, R = 1000, 
                     binary = TRUE, ncpus = 7, met_ref = met_silsrps)
boot_sirps$mod <- "si"
boot_sirps$method <- NA
write.csv(boot_sirps, file = paste0(out_dir, "boot_sirps.csv"), row.names = FALSE)


# weighted ----------------------------------------------------------------

## NLL

# CI-LS

# linear average
wboot_cilsnll_avgl <- get_bs(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, R = 1000, 
                             binary = TRUE, ncpus = 7, ncpus = 7,
                             weights = w_cilsnll_l, met_ref = met_sils)
wboot_cilsnll_avgl$mod <- "cils"
wboot_cilsnll_avgl$method <- "avg"
write.csv(wboot_cilsnll_avgl, file = paste0(out_dir, "wboot_cilsnll_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cilsnll_avgll <- get_bs(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, R = 1000, 
                              binary = TRUE, ncpus = 7,
                              weights = w_cilsnll_ll, met_ref = met_sils)
wboot_cilsnll_avgll$mod <- "cils"
wboot_cilsnll_avgll$method <- "avgll"
write.csv(wboot_cilsnll_avgll, file = paste0(out_dir, "wboot_cilsnll_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cilsnll_avgtrf <- get_bs(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, R = 1000,
                               binary = TRUE, ncpus = 7,
                               weights = w_cilsnll_trf, met_ref = met_sils)
wboot_cilsnll_avgtrf$mod <- "cils"
wboot_cilsnll_avgtrf$method <- "avgtrf"
write.csv(wboot_cilsnll_avgtrf, file = paste0(out_dir, "wboot_cilsnll_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cilsnll_l <- get_bs(lys_cdf_all = ens_w_cilsnll_l, y_true_all = y_true_all, R = 1000, 
                          binary = TRUE, ncpus = 7, met_ref = met_sils)
wboot_cilsnll_l$mod <- "cils"
wboot_cilsnll_l$method <- "linear"
write.csv(wboot_cilsnll_l, file = paste0(out_dir, "wboot_cilsnll_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cilsnll_ll <- get_bs(lys_cdf_all = ens_w_cilsnll_ll, y_true_all = y_true_all, R = 1000, 
                           binary = TRUE, ncpus = 7, met_ref = met_sils)
wboot_cilsnll_ll$mod <- "cils"
wboot_cilsnll_ll$method <- "log-linear"
write.csv(wboot_cilsnll_ll, file = paste0(out_dir, "wboot_cilsnll_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cilsnll_trf <- get_bs(lys_cdf_all = ens_w_cilsnll_trf, y_true_all = y_true_all, R = 1000, 
                            binary = TRUE, ncpus = 7, met_ref = met_sils)
wboot_cilsnll_trf$mod <- "cils"
wboot_cilsnll_trf$method <- "trafo"
write.csv(wboot_cilsnll_trf, file = paste0(out_dir, "wboot_cilsnll_trf.csv"), row.names = FALSE)

# CI

# linear average
wboot_cinll_avgl <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000, 
                           binary = TRUE, ncpus = 7,
                           weights = w_cinll_l, met_ref = met_sils)
wboot_cinll_avgl$mod <- "ci"
wboot_cinll_avgl$method <- "avg"
write.csv(wboot_cinll_avgl, file = paste0(out_dir, "wboot_cinll_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cinll_avgll <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000, 
                            binary = TRUE, ncpus = 7,
                            weights = w_cinll_ll, met_ref = met_sils)
wboot_cinll_avgll$mod <- "ci"
wboot_cinll_avgll$method <- "avgll"
write.csv(wboot_cinll_avgll, file = paste0(out_dir, "wboot_cinll_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cinll_avgtrf <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000, 
                             binary = TRUE, ncpus = 7,
                             weights = w_cinll_trf, met_ref = met_sils)
wboot_cinll_avgtrf$mod <- "ci"
wboot_cinll_avgtrf$method <- "avgtrf"
write.csv(wboot_cinll_avgtrf, file = paste0(out_dir, "wboot_cinll_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cinll_l <- get_bs(lys_cdf_all = ens_w_cinll_l, y_true_all = y_true_all, R = 1000, 
                        binary = TRUE, ncpus = 7, met_ref = met_sils)
wboot_cinll_l$mod <- "ci"
wboot_cinll_l$method <- "linear"
write.csv(wboot_cinll_l, file = paste0(out_dir, "wboot_cinll_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cinll_ll <- get_bs(lys_cdf_all = ens_w_cinll_ll, y_true_all = y_true_all, R = 1000, 
                         binary = TRUE, ncpus = 7, met_ref = met_sils)
wboot_cinll_ll$mod <- "ci"
wboot_cinll_ll$method <- "log-linear"
write.csv(wboot_cinll_ll, file = paste0(out_dir, "wboot_cinll_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cinll_trf <- get_bs(lys_cdf_all = ens_w_cinll_trf, y_true_all = y_true_all, R = 1000, 
                          binary = TRUE, ncpus = 7, met_ref = met_sils)
wboot_cinll_trf$mod <- "ci"
wboot_cinll_trf$method <- "trafo"
write.csv(wboot_cinll_trf, file = paste0(out_dir, "wboot_cinll_trf.csv"), row.names = FALSE)


## RPS

# CI-LS

# linear average
wboot_cilsrps_avgl <- get_bs(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, R = 1000, 
                             binary = TRUE, ncpus = 7,
                             weights = w_cilsrps_l, met_ref = met_silsrps)
wboot_cilsrps_avgl$mod <- "cils"
wboot_cilsrps_avgl$method <- "avg"
write.csv(wboot_cilsrps_avgl, file = paste0(out_dir, "wboot_cilsrps_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cilsrps_avgll <- get_bs(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, R = 1000, 
                              binary = TRUE, ncpus = 7,
                              weights = w_cilsrps_ll, met_ref = met_silsrps)
wboot_cilsrps_avgll$mod <- "cils"
wboot_cilsrps_avgll$method <- "avgll"
write.csv(wboot_cilsrps_avgll, file = paste0(out_dir, "wboot_cilsrps_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cilsrps_avgtrf <- get_bs(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, R = 1000, 
                               binary = TRUE, ncpus = 7,
                               weights = w_cilsrps_trf, met_ref = met_silsrps)
wboot_cilsrps_avgtrf$mod <- "cils"
wboot_cilsrps_avgtrf$method <- "avgtrf"
write.csv(wboot_cilsrps_avgtrf, file = paste0(out_dir, "wboot_cilsrps_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cilsrps_l <- get_bs(lys_cdf_all = ens_w_cilsrps_l, y_true_all = y_true_all, R = 1000, 
                          binary = TRUE, ncpus = 7, met_ref = met_silsrps)
wboot_cilsrps_l$mod <- "cils"
wboot_cilsrps_l$method <- "linear"
write.csv(wboot_cilsrps_l, file = paste0(out_dir, "wboot_cilsrps_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cilsrps_ll <- get_bs(lys_cdf_all = ens_w_cilsrps_ll, y_true_all = y_true_all, R = 1000, 
                           binary = TRUE, ncpus = 7, met_ref = met_silsrps)
wboot_cilsrps_ll$mod <- "cils"
wboot_cilsrps_ll$method <- "log-linear"
write.csv(wboot_cilsrps_ll, file = paste0(out_dir, "wboot_cilsrps_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cilsrps_trf <- get_bs(lys_cdf_all = ens_w_cilsrps_trf, y_true_all = y_true_all, R = 1000, 
                            binary = TRUE, ncpus = 7, met_ref = met_silsrps)
wboot_cilsrps_trf$mod <- "cils"
wboot_cilsrps_trf$method <- "trafo"
write.csv(wboot_cilsrps_trf, file = paste0(out_dir, "wboot_cilsrps_trf.csv"), row.names = FALSE)

# CI

# linear average
wboot_cirps_avgl <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, 
                           binary = TRUE, ncpus = 7,
                           weights = w_cirps_l, met_ref = met_silsrps)
wboot_cirps_avgl$mod <- "ci"
wboot_cirps_avgl$method <- "avg"
write.csv(wboot_cirps_avgl, file = paste0(out_dir, "wboot_cirps_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cirps_avgll <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, 
                            binary = TRUE, ncpus = 7,
                            weights = w_cirps_ll, met_ref = met_silsrps)
wboot_cirps_avgll$mod <- "ci"
wboot_cirps_avgll$method <- "avgll"
write.csv(wboot_cirps_avgll, file = paste0(out_dir, "wboot_cirps_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cirps_avgtrf <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, 
                             binary = TRUE, ncpus = 7,
                             weights = w_cirps_trf, met_ref = met_silsrps)
wboot_cirps_avgtrf$mod <- "ci"
wboot_cirps_avgtrf$method <- "avgtrf"
write.csv(wboot_cirps_avgtrf, file = paste0(out_dir, "wboot_cirps_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cirps_l <- get_bs(lys_cdf_all = ens_w_cirps_l, y_true_all = y_true_all, R = 1000, 
                        binary = TRUE, ncpus = 7, met_ref = met_silsrps)
wboot_cirps_l$mod <- "ci"
wboot_cirps_l$method <- "linear"
write.csv(wboot_cirps_l, file = paste0(out_dir, "wboot_cirps_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cirps_ll <- get_bs(lys_cdf_all = ens_w_cirps_ll, y_true_all = y_true_all, R = 1000, 
                         binary = TRUE, ncpus = 7, met_ref = met_silsrps)
wboot_cirps_ll$mod <- "ci"
wboot_cirps_ll$method <- "log-linear"
write.csv(wboot_cirps_ll, file = paste0(out_dir, "wboot_cirps_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cirps_trf <- get_bs(lys_cdf_all = ens_w_cirps_trf, y_true_all = y_true_all, R = 1000,
                          binary = TRUE, ncpus = 7, met_ref = met_silsrps)
wboot_cirps_trf$mod <- "ci"
wboot_cirps_trf$method <- "trafo"
write.csv(wboot_cirps_trf, file = paste0(out_dir, "wboot_cirps_trf.csv"), row.names = FALSE)
