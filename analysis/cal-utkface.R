# Calibration utkface data
# Andrea Goetschi
# April 2022

# Deps --------------------------------------------------------------------

library(etram)
library(caret)
library(readr)
library(rhdf5)
library(Hmisc)

# Params ------------------------------------------------------------------

ensembles <- 5
spl <- 6

source("experiments/functions/functions_DE.R")

path <- "~/../data/UTKFace/UTKFace.h5"
in_dir <- out_dir <- "experiments/results/DE/UTKFace/"

fname_silscsnll <- "utkface_silscs_lossnll_wsyes_augno"
fname_silscsrps <- "utkface_silscs_lossrps_wsyes_augno"
fname_cilsnll <- "utkface_cils_lossnll_wsyes_augno"
fname_cilsrps <- "utkface_cils_lossrps_wsyes_augno"
fname_sicsnll <- "utkface_sics_lossnll_wsyes_augno"
fname_sicsrps <- "utkface_sics_lossrps_wsyes_augno"
fname_cinll <- "utkface_ci_lossnll_wsyes_augno"
fname_cirps <- "utkface_ci_lossrps_wsyes_augno"

fname_sils <- "utkface_sils"
fname_si <- "utkface_si"


# Functions ---------------------------------------------------------------

cut_fun <- function(x) {
  q1 <- quantile(x, 0.05)
  q2 <- quantile(x, 0.95)
  cstep <- (q2 - q1) / 9
  between <- seq(from = q1, to = q2, by = cstep)
  cuts <- c(0, between)
  return(cuts)
}


# Read data ---------------------------------------------------------------

dat <- load_data("utkface", path = path)
tab_dat <- dat$tab_dat
y <- model.matrix(~ 0 + age_group, data = tab_dat)

ridx <- get_ridx(in_dir, fname = "utkface")

# Load results ------------------------------------------------------------

## CDFs all spl

### SI-LS-CS
cdftest_silscsnll <- list_cdfs(in_dir, fname_silscsnll, spl, ensembles, "test")
cdftest_silscsrps <- list_cdfs(in_dir, fname_silscsrps, spl, ensembles, "test")

### CI-LS
cdftest_cilsnll <- list_cdfs(in_dir, fname_cilsnll, spl, ensembles, "test")
cdftest_cilsrps <- list_cdfs(in_dir, fname_cilsrps, spl, ensembles, "test")

### SI-CS
cdftest_sicsnll <- list_cdfs(in_dir, fname_sicsnll, spl, ensembles, "test")
cdftest_sicsrps <- list_cdfs(in_dir, fname_sicsrps, spl, ensembles, "test")

### CI
cdftest_cinll <- list_cdfs(in_dir, fname_cinll, spl, ensembles, "test")
cdftest_cirps <- list_cdfs(in_dir, fname_cirps, spl, ensembles, "test")

### SI
cdftest_si <- list_cdfs(in_dir, fname_si, spl, ensembles = NULL, "test")

### SILS
cdftest_sils <- list_cdfs(in_dir, fname_sils, spl, ensembles = NULL, "test")


## Y true

ytest_1 <- y[ridx[ridx$spl == 1 & ridx$type == "test", "idx"], ]
ytest_2 <- y[ridx[ridx$spl == 2 & ridx$type == "test", "idx"], ]
ytest_3 <- y[ridx[ridx$spl == 3 & ridx$type == "test", "idx"], ]
ytest_4 <- y[ridx[ridx$spl == 4 & ridx$type == "test", "idx"], ]
ytest_5 <- y[ridx[ridx$spl == 5 & ridx$type == "test", "idx"], ]
ytest_6 <- y[ridx[ridx$spl == 6 & ridx$type == "test", "idx"], ]
y_true_all <- list(ytest_1, ytest_2, ytest_3, ytest_4, ytest_5, ytest_6)

## Weights 

### SILSCS
w_silscsnll <- read.csv(paste0(in_dir, "w_", fname_silscsnll, ".csv"))
w_silscsnll_l <- extract_w(w_silscsnll, meth = "linear")
w_silscsnll_ll <- extract_w(w_silscsnll, meth = "log-linear")
w_silscsnll_trf <- extract_w(w_silscsnll, meth = "trafo")

w_silscsrps <- read.csv(paste0(in_dir, "w_", fname_silscsrps, ".csv"))
w_silscsrps_l <- extract_w(w_silscsrps, meth = "linear")
w_silscsrps_ll <- extract_w(w_silscsrps, meth = "log-linear")
w_silscsrps_trf <- extract_w(w_silscsrps, meth = "trafo")

### CILS
w_cilsnll <- read.csv(paste0(in_dir, "w_", fname_cilsnll, ".csv"))
w_cilsnll_l <- extract_w(w_cilsnll, meth = "linear")
w_cilsnll_ll <- extract_w(w_cilsnll, meth = "log-linear")
w_cilsnll_trf <- extract_w(w_cilsnll, meth = "trafo")

w_cilsrps <- read.csv(paste0(in_dir, "w_", fname_cilsrps, ".csv"))
w_cilsrps_l <- extract_w(w_cilsrps, meth = "linear")
w_cilsrps_ll <- extract_w(w_cilsrps, meth = "log-linear")
w_cilsrps_trf <- extract_w(w_cilsrps, meth = "trafo")

### SICS
w_sicsnll <- read.csv(paste0(in_dir, "w_", fname_sicsnll, ".csv"))
w_sicsnll_l <- extract_w(w_sicsnll, meth = "linear")
w_sicsnll_ll <- extract_w(w_sicsnll, meth = "log-linear")
w_sicsnll_trf <- extract_w(w_sicsnll, meth = "trafo")

w_sicsrps <- read.csv(paste0(in_dir, "w_", fname_sicsrps, ".csv"))
w_sicsrps_l <- extract_w(w_sicsrps, meth = "linear")
w_sicsrps_ll <- extract_w(w_sicsrps, meth = "log-linear")
w_sicsrps_trf <- extract_w(w_sicsrps, meth = "trafo")

### CI
w_cinll <- read.csv(paste0(in_dir, "w_", fname_cinll, ".csv"))
w_cinll_l <- extract_w(w_cinll, meth = "linear")
w_cinll_ll <- extract_w(w_cinll, meth = "log-linear")
w_cinll_trf <- extract_w(w_cinll, meth = "trafo")

w_cirps <- read.csv(paste0(in_dir, "w_", fname_cirps, ".csv"))
w_cirps_l <- extract_w(w_cirps, meth = "linear")
w_cirps_ll <- extract_w(w_cirps, meth = "log-linear")
w_cirps_trf <- extract_w(w_cirps, meth = "trafo")


# Ensembles ---------------------------------------------------------------

# Equal weights -----------------------------------------------------------


## SILSCS

ens_silscsnll_l <- lapply(cdftest_silscsnll, get_ensemble, type = "linear")
ens_silscsnll_ll <- lapply(cdftest_silscsnll, get_ensemble, type = "log-linear")
ens_silscsnll_trf <- lapply(cdftest_silscsnll, get_ensemble, type = "trafo")

ens_silscsrps_l <- lapply(cdftest_silscsrps, get_ensemble, type = "linear")
ens_silscsrps_ll <- lapply(cdftest_silscsrps, get_ensemble, type = "log-linear")
ens_silscsrps_trf <- lapply(cdftest_silscsrps, get_ensemble, type = "trafo")

## CILS

ens_cilsnll_l <- lapply(cdftest_cilsnll, get_ensemble, type = "linear")
ens_cilsnll_ll <- lapply(cdftest_cilsnll, get_ensemble, type = "log-linear")
ens_cilsnll_trf <- lapply(cdftest_cilsnll, get_ensemble, type = "trafo")

ens_cilsrps_l <- lapply(cdftest_cilsrps, get_ensemble, type = "linear")
ens_cilsrps_ll <- lapply(cdftest_cilsrps, get_ensemble, type = "log-linear")
ens_cilsrps_trf <- lapply(cdftest_cilsrps, get_ensemble, type = "trafo")

## SICS

ens_sicsnll_l <- lapply(cdftest_sicsnll, get_ensemble, type = "linear")
ens_sicsnll_ll <- lapply(cdftest_sicsnll, get_ensemble, type = "log-linear")
ens_sicsnll_trf <- lapply(cdftest_sicsnll, get_ensemble, type = "trafo")

ens_sicsrps_l <- lapply(cdftest_sicsrps, get_ensemble, type = "linear")
ens_sicsrps_ll <- lapply(cdftest_sicsrps, get_ensemble, type = "log-linear")
ens_sicsrps_trf <- lapply(cdftest_sicsrps, get_ensemble, type = "trafo")

## CI

ens_cinll_l <- lapply(cdftest_cinll, get_ensemble, type = "linear")
ens_cinll_ll <- lapply(cdftest_cinll, get_ensemble, type = "log-linear")
ens_cinll_trf <- lapply(cdftest_cinll, get_ensemble, type = "trafo")

ens_cirps_l <- lapply(cdftest_cirps, get_ensemble, type = "linear")
ens_cirps_ll <- lapply(cdftest_cirps, get_ensemble, type = "log-linear")
ens_cirps_trf <- lapply(cdftest_cirps, get_ensemble, type = "trafo")


# Weighted ----------------------------------------------------------------

## SILSCS

ens_w_silscsnll_l <- get_weighted_ens(lys_cdf_all = cdftest_silscsnll, lys_weigths = w_silscsnll_l, type = "linear")
ens_w_silscsnll_ll <- get_weighted_ens(lys_cdf_all = cdftest_silscsnll, lys_weigths = w_silscsnll_ll, type = "log-linear")
ens_w_silscsnll_trf <- get_weighted_ens(lys_cdf_all = cdftest_silscsnll, lys_weigths = w_silscsnll_trf, type = "trafo")

ens_w_silscsrps_l <- get_weighted_ens(lys_cdf_all = cdftest_silscsrps, lys_weigths = w_silscsrps_l, type = "linear")
ens_w_silscsrps_ll <- get_weighted_ens(lys_cdf_all = cdftest_silscsrps, lys_weigths = w_silscsrps_ll, type = "log-linear")
ens_w_silscsrps_trf <- get_weighted_ens(lys_cdf_all = cdftest_silscsrps, lys_weigths = w_silscsrps_trf, type = "trafo")

## CILS

ens_w_cilsnll_l <- get_weighted_ens(lys_cdf_all = cdftest_cilsnll, lys_weigths = w_cilsnll_l, type = "linear")
ens_w_cilsnll_ll <- get_weighted_ens(lys_cdf_all = cdftest_cilsnll, lys_weigths = w_cilsnll_ll, type = "log-linear")
ens_w_cilsnll_trf <- get_weighted_ens(lys_cdf_all = cdftest_cilsnll, lys_weigths = w_cilsnll_trf, type = "trafo")

ens_w_cilsrps_l <- get_weighted_ens(lys_cdf_all = cdftest_cilsrps, lys_weigths = w_cilsrps_l, type = "linear")
ens_w_cilsrps_ll <- get_weighted_ens(lys_cdf_all = cdftest_cilsrps, lys_weigths = w_cilsrps_ll, type = "log-linear")
ens_w_cilsrps_trf <- get_weighted_ens(lys_cdf_all = cdftest_cilsrps, lys_weigths = w_cilsrps_trf, type = "trafo")

## SICS

ens_w_sicsnll_l <- get_weighted_ens(lys_cdf_all = cdftest_sicsnll, lys_weigths = w_sicsnll_l, type = "linear")
ens_w_sicsnll_ll <- get_weighted_ens(lys_cdf_all = cdftest_sicsnll, lys_weigths = w_sicsnll_ll, type = "log-linear")
ens_w_sicsnll_trf <- get_weighted_ens(lys_cdf_all = cdftest_sicsnll, lys_weigths = w_sicsnll_trf, type = "trafo")

ens_w_sicsrps_l <- get_weighted_ens(lys_cdf_all = cdftest_sicsrps, lys_weigths = w_sicsrps_l, type = "linear")
ens_w_sicsrps_ll <- get_weighted_ens(lys_cdf_all = cdftest_sicsrps, lys_weigths = w_sicsrps_ll, type = "log-linear")
ens_w_sicsrps_trf <- get_weighted_ens(lys_cdf_all = cdftest_sicsrps, lys_weigths = w_sicsrps_trf, type = "trafo")

## CI

ens_w_cinll_l <- get_weighted_ens(lys_cdf_all = cdftest_cinll, lys_weigths = w_cinll_l, type = "linear")
ens_w_cinll_ll <- get_weighted_ens(lys_cdf_all = cdftest_cinll, lys_weigths = w_cinll_ll, type = "log-linear")
ens_w_cinll_trf <- get_weighted_ens(lys_cdf_all = cdftest_cinll, lys_weigths = w_cinll_trf, type = "trafo")

ens_w_cirps_l <- get_weighted_ens(lys_cdf_all = cdftest_cirps, lys_weigths = w_cirps_l, type = "linear")
ens_w_cirps_ll <- get_weighted_ens(lys_cdf_all = cdftest_cirps, lys_weigths = w_cirps_ll, type = "log-linear")
ens_w_cirps_trf <- get_weighted_ens(lys_cdf_all = cdftest_cirps, lys_weigths = w_cirps_trf, type = "trafo")


# Calibration for all splits (ensembles) ----------------------------------

# Equal weights -----------------------------------------------------------

## SILSCS

cal_spl_silscsnll_l <- comb_ens_cal(lys_cdf_ens = ens_silscsnll_l, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_silscsnll_l$method <- "linear"
cal_spl_silscsnll_ll <- comb_ens_cal(lys_cdf_ens = ens_silscsnll_ll, y_true_all = y_true_all,
                                     emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_silscsnll_ll$method <- "log-linear"
cal_spl_silscsnll_trf <- comb_ens_cal(lys_cdf_ens = ens_silscsnll_trf, y_true_all = y_true_all,
                                      emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_silscsnll_trf$method <- "trafo"
c_spl_silscsnll <- bindr(pat1 = "cal_spl_silscsnll")
c_spl_silscsnll$mod <- "silscs"

cal_spl_silscsrps_l <- comb_ens_cal(lys_cdf_ens = ens_silscsrps_l, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_silscsrps_l$method <- "linear"
cal_spl_silscsrps_ll <- comb_ens_cal(lys_cdf_ens = ens_silscsrps_ll, y_true_all = y_true_all,
                                     emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_silscsrps_ll$method <- "log-linear"
cal_spl_silscsrps_trf <- comb_ens_cal(lys_cdf_ens = ens_silscsrps_trf, y_true_all = y_true_all,
                                      emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_silscsrps_trf$method <- "trafo"
c_spl_silscsrps <- bindr(pat1 = "cal_spl_silscsrps")
c_spl_silscsrps$mod <- "silscs"

## CILS

cal_spl_cilsnll_l <- comb_ens_cal(lys_cdf_ens = ens_cilsnll_l, y_true_all = y_true_all,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cilsnll_l$method <- "linear"
cal_spl_cilsnll_ll <- comb_ens_cal(lys_cdf_ens = ens_cilsnll_ll, y_true_all = y_true_all,
                                   emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cilsnll_ll$method <- "log-linear"
cal_spl_cilsnll_trf <- comb_ens_cal(lys_cdf_ens = ens_cilsnll_trf, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cilsnll_trf$method <- "trafo"
c_spl_cilsnll <- bindr(pat1 = "cal_spl_cilsnll")
c_spl_cilsnll$mod <- "cils"

cal_spl_cilsrps_l <- comb_ens_cal(lys_cdf_ens = ens_cilsrps_l, y_true_all = y_true_all,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cilsrps_l$method <- "linear"
cal_spl_cilsrps_ll <- comb_ens_cal(lys_cdf_ens = ens_cilsrps_ll, y_true_all = y_true_all,
                                   emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cilsrps_ll$method <- "log-linear"
cal_spl_cilsrps_trf <- comb_ens_cal(lys_cdf_ens = ens_cilsrps_trf, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cilsrps_trf$method <- "trafo"
c_spl_cilsrps <- bindr(pat1 = "cal_spl_cilsrps")
c_spl_cilsrps$mod <- "cils"

## SICS

cal_spl_sicsnll_l <- comb_ens_cal(lys_cdf_ens = ens_sicsnll_l, y_true_all = y_true_all,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_sicsnll_l$method <- "linear"
cal_spl_sicsnll_ll <- comb_ens_cal(lys_cdf_ens = ens_sicsnll_ll, y_true_all = y_true_all,
                                   emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_sicsnll_ll$method <- "log-linear"
cal_spl_sicsnll_trf <- comb_ens_cal(lys_cdf_ens = ens_sicsnll_trf, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_sicsnll_trf$method <- "trafo"
c_spl_sicsnll <- bindr(pat1 = "cal_spl_sicsnll")
c_spl_sicsnll$mod <- "sics"

cal_spl_sicsrps_l <- comb_ens_cal(lys_cdf_ens = ens_sicsrps_l, y_true_all = y_true_all,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_sicsrps_l$method <- "linear"
cal_spl_sicsrps_ll <- comb_ens_cal(lys_cdf_ens = ens_sicsrps_ll, y_true_all = y_true_all,
                                   emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_sicsrps_ll$method <- "log-linear"
cal_spl_sicsrps_trf <- comb_ens_cal(lys_cdf_ens = ens_sicsrps_trf, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_sicsrps_trf$method <- "trafo"
c_spl_sicsrps <- bindr(pat1 = "cal_spl_sicsrps")
c_spl_sicsrps$mod <- "sics"

## CI

cal_spl_cinll_l <- comb_ens_cal(lys_cdf_ens = ens_cinll_l, y_true_all = y_true_all,
                                emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cinll_l$method <- "linear"
cal_spl_cinll_ll <- comb_ens_cal(lys_cdf_ens = ens_cinll_ll, y_true_all = y_true_all,
                                 emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cinll_ll$method <- "log-linear"
cal_spl_cinll_trf <- comb_ens_cal(lys_cdf_ens = ens_cinll_trf, y_true_all = y_true_all,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cinll_trf$method <- "trafo"
c_spl_cinll <- bindr(pat1 = "cal_spl_cinll")
c_spl_cinll$mod <- "ci"

cal_spl_cirps_l <- comb_ens_cal(lys_cdf_ens = ens_cirps_l, y_true_all = y_true_all,
                                emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cirps_l$method <- "linear"
cal_spl_cirps_ll <- comb_ens_cal(lys_cdf_ens = ens_cirps_ll, y_true_all = y_true_all,
                                 emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cirps_ll$method <- "log-linear"
cal_spl_cirps_trf <- comb_ens_cal(lys_cdf_ens = ens_cirps_trf, y_true_all = y_true_all,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cirps_trf$method <- "trafo"
c_spl_cirps <- bindr(pat1 = "cal_spl_cirps")
c_spl_cirps$mod <- "ci"

###### Combine models

cal_spl_nll <- bindr(pat1 = "c_spl_", pat2 = "nll")
cal_spl_nll$weights <- "equal"
cal_spl_rps <- bindr(pat1 = "c_spl_", pat2 = "rps")
cal_spl_rps$weights <- "equal"


# Weighted ----------------------------------------------------------------

## SILSCS

cal_w_spl_silscsnll_l <- comb_ens_cal(lys_cdf_ens = ens_w_silscsnll_l, y_true_all = y_true_all,
                                      emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_silscsnll_l$method <- "linear"
cal_w_spl_silscsnll_ll <- comb_ens_cal(lys_cdf_ens = ens_w_silscsnll_ll, y_true_all = y_true_all,
                                       emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_silscsnll_ll$method <- "log-linear"
cal_w_spl_silscsnll_trf <- comb_ens_cal(lys_cdf_ens = ens_w_silscsnll_trf, y_true_all = y_true_all,
                                        emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_silscsnll_trf$method <- "trafo"
c_w_spl_silscsnll <- bindr(pat1 = "cal_w_spl_silscsnll")
c_w_spl_silscsnll$mod <- "silscs"

cal_w_spl_silscsrps_l <- comb_ens_cal(lys_cdf_ens = ens_w_silscsrps_l, y_true_all = y_true_all,
                                      emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_silscsrps_l$method <- "linear"
cal_w_spl_silscsrps_ll <- comb_ens_cal(lys_cdf_ens = ens_w_silscsrps_ll, y_true_all = y_true_all,
                                       emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_silscsrps_ll$method <- "log-linear"
cal_w_spl_silscsrps_trf <- comb_ens_cal(lys_cdf_ens = ens_w_silscsrps_trf, y_true_all = y_true_all,
                                        emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_silscsrps_trf$method <- "trafo"
c_w_spl_silscsrps <- bindr(pat1 = "cal_w_spl_silscsrps")
c_w_spl_silscsrps$mod <- "silscs"

## CILS

cal_w_spl_cilsnll_l <- comb_ens_cal(lys_cdf_ens = ens_w_cilsnll_l, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cilsnll_l$method <- "linear"
cal_w_spl_cilsnll_ll <- comb_ens_cal(lys_cdf_ens = ens_w_cilsnll_ll, y_true_all = y_true_all,
                                     emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cilsnll_ll$method <- "log-linear"
cal_w_spl_cilsnll_trf <- comb_ens_cal(lys_cdf_ens = ens_w_cilsnll_trf, y_true_all = y_true_all,
                                      emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cilsnll_trf$method <- "trafo"
c_w_spl_cilsnll <- bindr(pat1 = "cal_w_spl_cilsnll")
c_w_spl_cilsnll$mod <- "cils"

cal_w_spl_cilsrps_l <- comb_ens_cal(lys_cdf_ens = ens_w_cilsrps_l, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cilsrps_l$method <- "linear"
cal_w_spl_cilsrps_ll <- comb_ens_cal(lys_cdf_ens = ens_w_cilsrps_ll, y_true_all = y_true_all,
                                     emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cilsrps_ll$method <- "log-linear"
cal_w_spl_cilsrps_trf <- comb_ens_cal(lys_cdf_ens = ens_w_cilsrps_trf, y_true_all = y_true_all,
                                      emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cilsrps_trf$method <- "trafo"
c_w_spl_cilsrps <- bindr(pat1 = "cal_w_spl_cilsrps")
c_w_spl_cilsrps$mod <- "cils"

## SICS

cal_w_spl_sicsnll_l <- comb_ens_cal(lys_cdf_ens = ens_w_sicsnll_l, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_sicsnll_l$method <- "linear"
cal_w_spl_sicsnll_ll <- comb_ens_cal(lys_cdf_ens = ens_w_sicsnll_ll, y_true_all = y_true_all,
                                     emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_sicsnll_ll$method <- "log-linear"
cal_w_spl_sicsnll_trf <- comb_ens_cal(lys_cdf_ens = ens_w_sicsnll_trf, y_true_all = y_true_all,
                                      emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_sicsnll_trf$method <- "trafo"
c_w_spl_sicsnll <- bindr(pat1 = "cal_w_spl_sicsnll")
c_w_spl_sicsnll$mod <- "sics"

cal_w_spl_sicsrps_l <- comb_ens_cal(lys_cdf_ens = ens_w_sicsrps_l, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_sicsrps_l$method <- "linear"
cal_w_spl_sicsrps_ll <- comb_ens_cal(lys_cdf_ens = ens_w_sicsrps_ll, y_true_all = y_true_all,
                                     emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_sicsrps_ll$method <- "log-linear"
cal_w_spl_sicsrps_trf <- comb_ens_cal(lys_cdf_ens = ens_w_sicsrps_trf, y_true_all = y_true_all,
                                      emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_sicsrps_trf$method <- "trafo"
c_w_spl_sicsrps <- bindr(pat1 = "cal_w_spl_sicsrps")
c_w_spl_sicsrps$mod <- "sics"

## CI

cal_w_spl_cinll_l <- comb_ens_cal(lys_cdf_ens = ens_w_cinll_l, y_true_all = y_true_all,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cinll_l$method <- "linear"
cal_w_spl_cinll_ll <- comb_ens_cal(lys_cdf_ens = ens_w_cinll_ll, y_true_all = y_true_all,
                                   emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cinll_ll$method <- "log-linear"
cal_w_spl_cinll_trf <- comb_ens_cal(lys_cdf_ens = ens_w_cinll_trf, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cinll_trf$method <- "trafo"
c_w_spl_cinll <- bindr(pat1 = "cal_w_spl_cinll")
c_w_spl_cinll$mod <- "ci"

cal_w_spl_cirps_l <- comb_ens_cal(lys_cdf_ens = ens_w_cirps_l, y_true_all = y_true_all,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cirps_l$method <- "linear"
cal_w_spl_cirps_ll <- comb_ens_cal(lys_cdf_ens = ens_w_cirps_ll, y_true_all = y_true_all,
                                   emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cirps_ll$method <- "log-linear"
cal_w_spl_cirps_trf <- comb_ens_cal(lys_cdf_ens = ens_w_cirps_trf, y_true_all = y_true_all,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cirps_trf$method <- "trafo"
c_w_spl_cirps <- bindr(pat1 = "cal_w_spl_cirps")
c_w_spl_cirps$mod <- "ci"

###### Combine models

cal_w_spl_nll <- bindr(pat1 = "c_w_spl_", pat2 = "nll")
cal_w_spl_nll$weights <- "tuned"
cal_w_spl_rps <- bindr(pat1 = "c_w_spl_", pat2 = "rps")
cal_w_spl_rps$weights <- "tuned"

###### Combine weigted non-weighted

cal_spl_nll <- rbind(cal_spl_nll, cal_w_spl_nll)
cal_spl_nll$mod <- factor(cal_spl_nll$mod)
cal_spl_nll$weights <- factor(cal_spl_nll$weights)
cal_spl_nll$method <- factor(cal_spl_nll$method)

cal_spl_rps <- rbind(cal_spl_rps, cal_w_spl_rps)
cal_spl_rps$mod <- factor(cal_spl_rps$mod)
cal_spl_rps$weights <- factor(cal_spl_rps$weights)
cal_spl_rps$method <- factor(cal_spl_rps$method)


# Mean calibration of members per split -----------------------------------

# Equal weights (AVG-LIN) -------------------------------------------------

## SILSCS

avg_silscsnll <- comb_avg_cal(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all,
                              emp = TRUE, custom_cuts_fun = cut_fun)
avg_silscsnll$mod <- "silscs"
avg_silscsrps <- comb_avg_cal(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all,
                              emp = TRUE, custom_cuts_fun = cut_fun)
avg_silscsrps$mod <- "silscs"

## CILS

avg_cilsnll <- comb_avg_cal(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
                            emp = TRUE, custom_cuts_fun = cut_fun)
avg_cilsnll$mod <- "cils"
avg_cilsrps <- comb_avg_cal(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
                            emp = TRUE, custom_cuts_fun = cut_fun)
avg_cilsrps$mod <- "cils"

## SICS

avg_sicsnll <- comb_avg_cal(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all,
                            emp = TRUE, custom_cuts_fun = cut_fun)
avg_sicsnll$mod <- "sics"
avg_sicsrps <- comb_avg_cal(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all,
                            emp = TRUE, custom_cuts_fun = cut_fun)
avg_sicsrps$mod <- "sics"

## CI

avg_cinll <- comb_avg_cal(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                          emp = TRUE, custom_cuts_fun = cut_fun)
avg_cinll$mod <- "ci"
avg_cirps <- comb_avg_cal(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                          emp = TRUE, custom_cuts_fun = cut_fun)
avg_cirps$mod <- "ci"

###### Combine models

cal_avg_spl_nll <- bindr(pat1 = "avg_", pat2 = "nll")
cal_avg_spl_nll$weights <- "equal"
cal_avg_spl_nll$method <- "avg"

cal_avg_spl_rps <- bindr(pat1 = "avg_", pat2 = "rps")
cal_avg_spl_rps$weights <- "equal"
cal_avg_spl_rps$method <- "avg"


# Weighted (AVG-LIN, AVG-LOG, AVG-TRF) ------------------------------------

## SILSCS

avg_w_silscsnll_l <- comb_avg_cal(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all,
                                  weights = w_silscsnll_l, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_silscsnll_l$mod <- "silscs"
avg_w_silscsnll_l$method <- "avg"

avg_w_silscsnll_ll <- comb_avg_cal(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all,
                                   weights = w_silscsnll_ll, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_silscsnll_ll$mod <- "silscs"
avg_w_silscsnll_ll$method <- "avgll"

avg_w_silscsnll_trf <- comb_avg_cal(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all,
                                    weights = w_silscsnll_trf, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_silscsnll_trf$mod <- "silscs"
avg_w_silscsnll_trf$method <- "avgtrf"


avg_w_silscsrps_l <- comb_avg_cal(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all,
                                  weights = w_silscsrps_l, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_silscsrps_l$mod <- "silscs"
avg_w_silscsrps_l$method <- "avg"

avg_w_silscsrps_ll <- comb_avg_cal(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all,
                                   weights = w_silscsrps_ll, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_silscsrps_ll$mod <- "silscs"
avg_w_silscsrps_ll$method <- "avgll"

avg_w_silscsrps_trf <- comb_avg_cal(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all,
                                    weights = w_silscsrps_trf, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_silscsrps_trf$mod <- "silscs"
avg_w_silscsrps_trf$method <- "avgtrf"

## CILS

avg_w_cilsnll_l <- comb_avg_cal(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
                                weights = w_cilsnll_l, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cilsnll_l$mod <- "cils"
avg_w_cilsnll_l$method <- "avg"

avg_w_cilsnll_ll <- comb_avg_cal(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
                                 weights = w_cilsnll_ll, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cilsnll_ll$mod <- "cils"
avg_w_cilsnll_ll$method <- "avgll"

avg_w_cilsnll_trf <- comb_avg_cal(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
                                  weights = w_cilsnll_trf, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cilsnll_trf$mod <- "cils"
avg_w_cilsnll_trf$method <- "avgtrf"


avg_w_cilsrps_l <- comb_avg_cal(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
                                weights = w_cilsrps_l, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cilsrps_l$mod <- "cils"
avg_w_cilsrps_l$method <- "avg"

avg_w_cilsrps_ll <- comb_avg_cal(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
                                 weights = w_cilsrps_ll, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cilsrps_ll$mod <- "cils"
avg_w_cilsrps_ll$method <- "avgll"

avg_w_cilsrps_trf <- comb_avg_cal(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
                                  weights = w_cilsrps_trf, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cilsrps_trf$mod <- "cils"
avg_w_cilsrps_trf$method <- "avgtrf"

## SICS

avg_w_sicsnll_l <- comb_avg_cal(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all,
                                weights = w_sicsnll_l, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_sicsnll_l$mod <- "sics"
avg_w_sicsnll_l$method <- "avg"

avg_w_sicsnll_ll <- comb_avg_cal(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all,
                                 weights = w_sicsnll_ll, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_sicsnll_ll$mod <- "sics"
avg_w_sicsnll_ll$method <- "avgll"

avg_w_sicsnll_trf <- comb_avg_cal(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all,
                                  weights = w_sicsnll_trf, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_sicsnll_trf$mod <- "sics"
avg_w_sicsnll_trf$method <- "avgtrf"


avg_w_sicsrps_l <- comb_avg_cal(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all,
                                weights = w_sicsrps_l, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_sicsrps_l$mod <- "sics"
avg_w_sicsrps_l$method <- "avg"

avg_w_sicsrps_ll <- comb_avg_cal(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all,
                                 weights = w_sicsrps_ll, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_sicsrps_ll$mod <- "sics"
avg_w_sicsrps_ll$method <- "avgll"

avg_w_sicsrps_trf <- comb_avg_cal(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all,
                                  weights = w_sicsrps_trf, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_sicsrps_trf$mod <- "sics"
avg_w_sicsrps_trf$method <- "avgtrf"

## CI

avg_w_cinll_l <- comb_avg_cal(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                              weights = w_cinll_l, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cinll_l$mod <- "ci"
avg_w_cinll_l$method <- "avg"

avg_w_cinll_ll <- comb_avg_cal(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                               weights = w_cinll_ll, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cinll_ll$mod <- "ci"
avg_w_cinll_ll$method <- "avgll"

avg_w_cinll_trf <- comb_avg_cal(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                                weights = w_cinll_trf, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cinll_trf$mod <- "ci"
avg_w_cinll_trf$method <- "avgtrf"


avg_w_cirps_l <- comb_avg_cal(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                              weights = w_cirps_l, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cirps_l$mod <- "ci"
avg_w_cirps_l$method <- "avg"

avg_w_cirps_ll <- comb_avg_cal(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                               weights = w_cirps_ll, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cirps_ll$mod <- "ci"
avg_w_cirps_ll$method <- "avgll"

avg_w_cirps_trf <- comb_avg_cal(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                                weights = w_cirps_trf, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cirps_trf$mod <- "ci"
avg_w_cirps_trf$method <- "avgtrf"


###### Combine models

cal_w_avg_spl_nll <- bindr(pat1 = "avg_w_", pat2 = "nll")
cal_w_avg_spl_nll$weights <- "tuned"

cal_w_avg_spl_rps <- bindr(pat1 = "avg_w_", pat2 = "rps")
cal_w_avg_spl_rps$weights <- "tuned"

###### Combine weigted non-weighted

cal_avg_spl_nll <- rbind(cal_avg_spl_nll, cal_w_avg_spl_nll)
cal_avg_spl_nll$mod <- factor(cal_avg_spl_nll$mod)
cal_avg_spl_nll$weights <- factor(cal_avg_spl_nll$weights)
cal_avg_spl_nll$method <- factor(cal_avg_spl_nll$method)

cal_avg_spl_rps <- rbind(cal_avg_spl_rps, cal_w_avg_spl_rps)
cal_avg_spl_rps$mod <- factor(cal_avg_spl_rps$mod)
cal_avg_spl_rps$weights <- factor(cal_avg_spl_rps$weights)
cal_avg_spl_rps$method <- factor(cal_avg_spl_rps$method)


# Individual --------------------------------------------------------------

## SILSCS

lys_i_silscsnll <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = TRUE, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_silscsnll))) {
  lys_i_silscsnll[[i]]$spl <- factor(i)
}
indiv_silscsnll <- do.call("rbind", lys_i_silscsnll)
indiv_silscsnll$mod <- "silscs"


lys_i_silscsrps <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = TRUE, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_silscsrps))) {
  lys_i_silscsrps[[i]]$spl <- factor(i)
}
indiv_silscsrps <- do.call("rbind", lys_i_silscsrps)
indiv_silscsrps$mod <- "silscs"


## CILS

lys_i_cilsnll <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = TRUE, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_cilsnll))) {
  lys_i_cilsnll[[i]]$spl <- factor(i)
}
indiv_cilsnll <- do.call("rbind", lys_i_cilsnll)
indiv_cilsnll$mod <- "cils"


lys_i_cilsrps <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = TRUE, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_cilsrps))) {
  lys_i_cilsrps[[i]]$spl <- factor(i)
}
indiv_cilsrps <- do.call("rbind", lys_i_cilsrps)
indiv_cilsrps$mod <- "cils"

## SICS

lys_i_sicsnll <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = TRUE, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_sicsnll))) {
  lys_i_sicsnll[[i]]$spl <- factor(i)
}
indiv_sicsnll <- do.call("rbind", lys_i_sicsnll)
indiv_sicsnll$mod <- "sics"


lys_i_sicsrps <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = TRUE, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_sicsrps))) {
  lys_i_sicsrps[[i]]$spl <- factor(i)
}
indiv_sicsrps <- do.call("rbind", lys_i_sicsrps)
indiv_sicsrps$mod <- "sics"

## CI

lys_i_cinll <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = TRUE, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_cinll))) {
  lys_i_cinll[[i]]$spl <- factor(i)
}
indiv_cinll <- do.call("rbind", lys_i_cinll)
indiv_cinll$mod <- "ci"


lys_i_cirps <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = TRUE, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_cirps))) {
  lys_i_cirps[[i]]$spl <- factor(i)
}
indiv_cirps <- do.call("rbind", lys_i_cirps)
indiv_cirps$mod <- "ci"


###### Combine models

ind_nll <- bindr(pat1 = "indiv", pat2 = "nll")
ind_rps <- bindr(pat1 = "indiv", pat2 = "rps")


# Calibration of ref models -----------------------------------------------

# list of cdf (not nested list)
cdftest_si <- lapply(cdftest_si, function(x) do.call(rbind, x))
cdftest_sils <- lapply(cdftest_sils, function(x) do.call(rbind, x))

cal_spl_si <- comb_ens_cal(lys_cdf_ens = cdftest_si, y_true_all = y_true_all,
                           emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_si$mod <- "si"
cal_spl_si$weights <- "equal"
cal_spl_si$method <- NA
cal_w_spl_si <- cal_spl_si %>% mutate(weights = "tuned")

cal_spl_sils <- comb_ens_cal(lys_cdf_ens = cdftest_sils, y_true_all = y_true_all,
                             emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_sils$mod <- "sils"
cal_spl_sils$weights <- "equal"
cal_spl_sils$method <- NA
cal_w_spl_sils <- cal_spl_sils %>% mutate(weights = "tuned")

cal_spl_refs <- rbind(cal_spl_si, cal_w_spl_si,
                      cal_spl_sils, cal_w_spl_sils)


# Combine all calibrations of single splits -------------------------------

spl_nll <- rbind(cal_spl_nll, cal_avg_spl_nll, cal_spl_refs)
spl_rps <- rbind(cal_spl_rps, cal_avg_spl_rps, cal_spl_refs)


# Average cal across all splits (ensembles) -------------------------------

splitted_nll <- cal_spl_nll %>% group_split(method, weights, mod)
avg_meths_nll <- avg_across_spl(splitted_nll)

splitted_rps <- cal_spl_rps %>% group_split(method, weights, mod)
avg_meths_rps <- avg_across_spl(splitted_rps)

# Average cal across splits (avg methods) ---------------------------------

splitted_avg_nll <- cal_avg_spl_nll %>% group_split(method, weights, mod)
avg_avg_meths_nll <- avg_across_spl(splitted_avg_nll)

splitted_avg_rps <- cal_avg_spl_rps %>% group_split(method, weights, mod)
avg_avg_meths_rps <- avg_across_spl(splitted_avg_rps)

# Average across splits (refs) --------------------------------------------

splitted_refs <- cal_spl_refs %>% group_split(weights, mod)
avg_refs <- avg_across_spl(splitted_refs)

# Combine all average calibrations across splits --------------------------

avg_nll <- rbind(avg_meths_nll, avg_avg_meths_nll, avg_refs)
avg_rps <- rbind(avg_meths_rps, avg_avg_meths_rps, avg_refs)


# Save results ------------------------------------------------------------

write.csv(spl_nll, file = paste0(out_dir, "cal_splnll.csv"))
write.csv(spl_rps, file = paste0(out_dir, "cal_splrps.csv"))

write.csv(avg_nll, file = paste0(out_dir, "cal_avgnll.csv"))
write.csv(avg_rps, file = paste0(out_dir, "cal_avgrps.csv"))

write.csv(ind_nll, file = paste0(out_dir, "cal_indnll.csv"))
write.csv(ind_rps, file = paste0(out_dir, "cal_indrps.csv"))

