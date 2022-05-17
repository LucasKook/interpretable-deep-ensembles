# Bootstrap confidence intervals
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(etram)
library(boot)

# Directories -------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/MNIST/"

# Params ------------------------------------------------------------------

K <- 10

fname_cinll <- "mnist_ci_lossnll_wsno_augno"
fname_cirps <- "mnist_ci_lossrps_wsno_augno"

# Load results ------------------------------------------------------------

## all CDF, Y
all_cdf <- read.csv(paste0(in_dir, "mnist_merged_cdf_ci.csv"))
all_y <- read.csv(paste0(in_dir, "mnist_merged_y.csv"))

## CDFs all splits

### CI
cdftest_cinll <- load_lys_cdf_all(all_cdf, m = "ci", K = K, l = "nll", t = "test")
cdftest_cirps <- load_lys_cdf_all(all_cdf, m = "ci", K = K, l = "rps", t = "test")

## Y true

y_true_all <- load_y_true_all(all_y, K = K, t = "test")

### Weights

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

# CI
ens_cinll_l <- lapply(lapply(cdftest_cinll, get_ensemble, type = "linear"), list)
ens_cinll_ll <- lapply(lapply(cdftest_cinll, get_ensemble, type = "log-linear"), list)
ens_cinll_trf <- lapply(lapply(cdftest_cinll, get_ensemble, type = "trafo"), list)

ens_cirps_l <- lapply(lapply(cdftest_cirps, get_ensemble, type = "linear"), list)
ens_cirps_ll <- lapply(lapply(cdftest_cirps, get_ensemble, type = "log-linear"), list)
ens_cirps_trf <- lapply(lapply(cdftest_cirps, get_ensemble, type = "trafo"), list)

# Construct ensembles (weighted) ------------------------------------------

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

# 'normal average'
boot_cinll_avg <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000, binary = FALSE)
boot_cinll_avg$mod <- "ci"
boot_cinll_avg$method <- "avg"
write.csv(boot_cinll_avg, file = paste0(out_dir, "boot_cinll_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cinll_l <- get_bs(lys_cdf_all = ens_cinll_l, y_true_all = y_true_all, R = 1000, binary = FALSE)
boot_cinll_l$mod <- "ci"
boot_cinll_l$method <- "linear"
write.csv(boot_cinll_l, file = paste0(out_dir, "boot_cinll_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cinll_ll <- get_bs(lys_cdf_all = ens_cinll_ll, y_true_all = y_true_all, R = 1000, binary = FALSE)
boot_cinll_ll$mod <- "ci"
boot_cinll_ll$method <- "log-linear"
write.csv(boot_cinll_ll, file = paste0(out_dir, "boot_cinll_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cinll_trf <- get_bs(lys_cdf_all = ens_cinll_trf, y_true_all = y_true_all, R = 1000, binary = FALSE)
boot_cinll_trf$mod <- "ci"
boot_cinll_trf$method <- "trafo"
write.csv(boot_cinll_trf, file = paste0(out_dir, "boot_cinll_trf.csv"), row.names = FALSE)


## RPS

# 'normal average'
boot_cirps_avg <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, binary = FALSE)
boot_cirps_avg$mod <- "ci"
boot_cirps_avg$method <- "avg"
write.csv(boot_cirps_avg, file = paste0(out_dir, "boot_cirps_avg.csv"), row.names = FALSE)

# linear ensemble
boot_cirps_l <- get_bs(lys_cdf_all = ens_cirps_l, y_true_all = y_true_all, R = 1000, binary = FALSE)
boot_cirps_l$mod <- "ci"
boot_cirps_l$method <- "linear"
write.csv(boot_cirps_l, file = paste0(out_dir, "boot_cirps_l.csv"), row.names = FALSE)

# log-linear ensemble
boot_cirps_ll <- get_bs(lys_cdf_all = ens_cirps_ll, y_true_all = y_true_all, R = 1000, binary = FALSE)
boot_cirps_ll$mod <- "ci"
boot_cirps_ll$method <- "log-linear"
write.csv(boot_cirps_ll, file = paste0(out_dir, "boot_cirps_ll.csv"), row.names = FALSE)

# trafo ensemble
boot_cirps_trf <- get_bs(lys_cdf_all = ens_cirps_trf, y_true_all = y_true_all, R = 1000, binary = FALSE)
boot_cirps_trf$mod <- "ci"
boot_cirps_trf$method <- "trafo"
write.csv(boot_cirps_trf, file = paste0(out_dir, "boot_cirps_trf.csv"), row.names = FALSE)


# weighted ----------------------------------------------------------------

## NLL

# linear average
wboot_cinll_avgl <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000, binary = FALSE,
                           weights = w_cinll_l)
wboot_cinll_avgl$mod <- "ci"
wboot_cinll_avgl$method <- "avg"
write.csv(wboot_cinll_avgl, file = paste0(out_dir, "wboot_cinll_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cinll_avgll <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000, binary = FALSE,
                            weights = w_cinll_ll)
wboot_cinll_avgll$mod <- "ci"
wboot_cinll_avgll$method <- "avgll"
write.csv(wboot_cinll_avgll, file = paste0(out_dir, "wboot_cinll_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cinll_avgtrf <- get_bs(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, R = 1000, binary = FALSE,
                             weights = w_cinll_trf)
wboot_cinll_avgtrf$mod <- "ci"
wboot_cinll_avgtrf$method <- "avgtrf"
write.csv(wboot_cinll_avgtrf, file = paste0(out_dir, "wboot_cinll_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cinll_l <- get_bs(lys_cdf_all = ens_w_cinll_l, y_true_all = y_true_all, R = 1000, binary = FALSE)
wboot_cinll_l$mod <- "ci"
wboot_cinll_l$method <- "linear"
write.csv(wboot_cinll_l, file = paste0(out_dir, "wboot_cinll_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cinll_ll <- get_bs(lys_cdf_all = ens_w_cinll_ll, y_true_all = y_true_all, R = 1000, binary = FALSE)
wboot_cinll_ll$mod <- "ci"
wboot_cinll_ll$method <- "log-linear"
write.csv(wboot_cinll_ll, file = paste0(out_dir, "wboot_cinll_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cinll_trf <- get_bs(lys_cdf_all = ens_w_cinll_trf, y_true_all = y_true_all, R = 1000, binary = FALSE)
wboot_cinll_trf$mod <- "ci"
wboot_cinll_trf$method <- "trafo"
write.csv(wboot_cinll_trf, file = paste0(out_dir, "wboot_cinll_trf.csv"), row.names = FALSE)


## RPS

# linear average
wboot_cirps_avgl <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, binary = FALSE,
                           weights = w_cirps_l)
wboot_cirps_avgl$mod <- "ci"
wboot_cirps_avgl$method <- "avg"
write.csv(wboot_cirps_avgl, file = paste0(out_dir, "wboot_cirps_avgl.csv"), row.names = FALSE)

# log-linear average
wboot_cirps_avgll <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, binary = FALSE,
                            weights = w_cirps_ll)
wboot_cirps_avgll$mod <- "ci"
wboot_cirps_avgll$method <- "avgll"
write.csv(wboot_cirps_avgll, file = paste0(out_dir, "wboot_cirps_avgll.csv"), row.names = FALSE)

# trafo average
wboot_cirps_avgtrf <- get_bs(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, R = 1000, binary = FALSE,
                             weights = w_cirps_trf)
wboot_cirps_avgtrf$mod <- "ci"
wboot_cirps_avgtrf$method <- "avgtrf"
write.csv(wboot_cirps_avgtrf, file = paste0(out_dir, "wboot_cirps_avgtrf.csv"), row.names = FALSE)

# linear ensemble
wboot_cirps_l <- get_bs(lys_cdf_all = ens_w_cirps_l, y_true_all = y_true_all, R = 1000, binary = FALSE)
wboot_cirps_l$mod <- "ci"
wboot_cirps_l$method <- "linear"
write.csv(wboot_cirps_l, file = paste0(out_dir, "wboot_cirps_l.csv"), row.names = FALSE)

# log-linear ensemble
wboot_cirps_ll <- get_bs(lys_cdf_all = ens_w_cirps_ll, y_true_all = y_true_all, R = 1000, binary = FALSE)
wboot_cirps_ll$mod <- "ci"
wboot_cirps_ll$method <- "log-linear"
write.csv(wboot_cirps_ll, file = paste0(out_dir, "wboot_cirps_ll.csv"), row.names = FALSE)

# trafo ensemble
wboot_cirps_trf <- get_bs(lys_cdf_all = ens_w_cirps_trf, y_true_all = y_true_all, R = 1000, binary = FALSE)
wboot_cirps_trf$mod <- "ci"
wboot_cirps_trf$method <- "trafo"
write.csv(wboot_cirps_trf, file = paste0(out_dir, "wboot_cirps_trf.csv"), row.names = FALSE)
