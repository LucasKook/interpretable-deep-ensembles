# Calibration mnist data
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

in_dir <- out_dir <- "experiments/results/DE/MNIST/"

fname_cinll <- "mnist_ci_lossnll_wsno_augno"
fname_cirps <- "mnist_ci_lossrps_wsno_augno"

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

dat <- load_data("mnist")
tab_dat <- dat$tab_dat
y <- model.matrix(~ 0 + y, data = tab_dat)

ridx <- get_ridx(in_dir, fname = "mnist")


# Load results ------------------------------------------------------------

## CDFs all spl

### CI
cdftest_cinll <- list_cdfs(in_dir, fname_cinll, spl, ensembles, "test")
cdftest_cirps <- list_cdfs(in_dir, fname_cirps, spl, ensembles, "test")

## Y true

ytest_1 <- y[ridx[ridx$spl == 1 & ridx$type == "test", "idx"], ]
ytest_2 <- y[ridx[ridx$spl == 2 & ridx$type == "test", "idx"], ]
ytest_3 <- y[ridx[ridx$spl == 3 & ridx$type == "test", "idx"], ]
ytest_4 <- y[ridx[ridx$spl == 4 & ridx$type == "test", "idx"], ]
ytest_5 <- y[ridx[ridx$spl == 5 & ridx$type == "test", "idx"], ]
ytest_6 <- y[ridx[ridx$spl == 6 & ridx$type == "test", "idx"], ]
y_true_all <- list(ytest_1, ytest_2, ytest_3, ytest_4, ytest_5, ytest_6)

## Weights 

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

## CI

ens_cinll_l <- lapply(cdftest_cinll, get_ensemble, type = "linear")
ens_cinll_ll <- lapply(cdftest_cinll, get_ensemble, type = "log-linear")
ens_cinll_trf <- lapply(cdftest_cinll, get_ensemble, type = "trafo")

ens_cirps_l <- lapply(cdftest_cirps, get_ensemble, type = "linear")
ens_cirps_ll <- lapply(cdftest_cirps, get_ensemble, type = "log-linear")
ens_cirps_trf <- lapply(cdftest_cirps, get_ensemble, type = "trafo")


# Weighted ----------------------------------------------------------------

## CI

ens_w_cinll_l <- get_weighted_ens(lys_cdf_all = cdftest_cinll, lys_weigths = w_cinll_l, type = "linear")
ens_w_cinll_ll <- get_weighted_ens(lys_cdf_all = cdftest_cinll, lys_weigths = w_cinll_ll, type = "log-linear")
ens_w_cinll_trf <- get_weighted_ens(lys_cdf_all = cdftest_cinll, lys_weigths = w_cinll_trf, type = "trafo")

ens_w_cirps_l <- get_weighted_ens(lys_cdf_all = cdftest_cirps, lys_weigths = w_cirps_l, type = "linear")
ens_w_cirps_ll <- get_weighted_ens(lys_cdf_all = cdftest_cirps, lys_weigths = w_cirps_ll, type = "log-linear")
ens_w_cirps_trf <- get_weighted_ens(lys_cdf_all = cdftest_cirps, lys_weigths = w_cirps_trf, type = "trafo")


# Calibration for all splits (ensembles) ----------------------------------

# Equal weights -----------------------------------------------------------

## CI

cal_spl_cinll_l <- comb_ens_cal(lys_cdf_ens = ens_cinll_l, y_true_all = y_true_all, cumulative = TRUE,
                                emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cinll_l$method <- "linear"
cal_spl_cinll_ll <- comb_ens_cal(lys_cdf_ens = ens_cinll_ll, y_true_all = y_true_all, cumulative = TRUE,
                                 emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cinll_ll$method <- "log-linear"
cal_spl_cinll_trf <- comb_ens_cal(lys_cdf_ens = ens_cinll_trf, y_true_all = y_true_all, cumulative = TRUE,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cinll_trf$method <- "trafo"
c_spl_cinll <- bindr(pat1 = "cal_spl_cinll")
c_spl_cinll$mod <- "ci"

cal_spl_cirps_l <- comb_ens_cal(lys_cdf_ens = ens_cirps_l, y_true_all = y_true_all, cumulative = TRUE,
                                emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cirps_l$method <- "linear"
cal_spl_cirps_ll <- comb_ens_cal(lys_cdf_ens = ens_cirps_ll, y_true_all = y_true_all, cumulative = TRUE,
                                 emp = TRUE, custom_cuts_fun = cut_fun)
cal_spl_cirps_ll$method <- "log-linear"
cal_spl_cirps_trf <- comb_ens_cal(lys_cdf_ens = ens_cirps_trf, y_true_all = y_true_all, cumulative = TRUE,
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

## CI

cal_w_spl_cinll_l <- comb_ens_cal(lys_cdf_ens = ens_w_cinll_l, y_true_all = y_true_all, cumulative = TRUE,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cinll_l$method <- "linear"
cal_w_spl_cinll_ll <- comb_ens_cal(lys_cdf_ens = ens_w_cinll_ll, y_true_all = y_true_all, cumulative = TRUE,
                                   emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cinll_ll$method <- "log-linear"
cal_w_spl_cinll_trf <- comb_ens_cal(lys_cdf_ens = ens_w_cinll_trf, y_true_all = y_true_all, cumulative = TRUE,
                                    emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cinll_trf$method <- "trafo"
c_w_spl_cinll <- bindr(pat1 = "cal_w_spl_cinll")
c_w_spl_cinll$mod <- "ci"

cal_w_spl_cirps_l <- comb_ens_cal(lys_cdf_ens = ens_w_cirps_l, y_true_all = y_true_all, cumulative = TRUE,
                                  emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cirps_l$method <- "linear"
cal_w_spl_cirps_ll <- comb_ens_cal(lys_cdf_ens = ens_w_cirps_ll, y_true_all = y_true_all, cumulative = TRUE,
                                   emp = TRUE, custom_cuts_fun = cut_fun)
cal_w_spl_cirps_ll$method <- "log-linear"
cal_w_spl_cirps_trf <- comb_ens_cal(lys_cdf_ens = ens_w_cirps_trf, y_true_all = y_true_all, cumulative = TRUE,
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

## CI

avg_cinll <- comb_avg_cal(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, cumulative = TRUE,
                          emp = TRUE, custom_cuts_fun = cut_fun)
avg_cinll$mod <- "ci"
avg_cirps <- comb_avg_cal(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, cumulative = TRUE,
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

## CI

avg_w_cinll_l <- comb_avg_cal(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                              weights = w_cinll_l, cumulative = TRUE, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cinll_l$mod <- "ci"
avg_w_cinll_l$method <- "avg"

avg_w_cinll_ll <- comb_avg_cal(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                               weights = w_cinll_ll, cumulative = TRUE, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cinll_ll$mod <- "ci"
avg_w_cinll_ll$method <- "avgll"

avg_w_cinll_trf <- comb_avg_cal(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                                weights = w_cinll_trf, cumulative = TRUE, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cinll_trf$mod <- "ci"
avg_w_cinll_trf$method <- "avgtrf"


avg_w_cirps_l <- comb_avg_cal(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                              weights = w_cirps_l, cumulative = TRUE, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cirps_l$mod <- "ci"
avg_w_cirps_l$method <- "avg"

avg_w_cirps_ll <- comb_avg_cal(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                               weights = w_cirps_ll, cumulative = TRUE, emp = TRUE, custom_cuts_fun = cut_fun)
avg_w_cirps_ll$mod <- "ci"
avg_w_cirps_ll$method <- "avgll"

avg_w_cirps_trf <- comb_avg_cal(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                                weights = w_cirps_trf, cumulative = TRUE, emp = TRUE, custom_cuts_fun = cut_fun)
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

## CI

lys_i_cinll <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts, emp) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = emp, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, emp = TRUE, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_cinll))) {
  lys_i_cinll[[i]]$spl <- factor(i)
}
indiv_cinll <- do.call("rbind", lys_i_cinll)
indiv_cinll$mod <- "ci"


lys_i_cirps <- mapply(function(lys_cdf_all, y_true_all, cumulative,  cuts, emp) {
  cal_indiv(lys_cdf = lys_cdf_all, y_true = y_true_all, 
            cumulative = cumulative, cuts = cuts, emp = emp, custom_cuts_fun = cut_fun)
}, lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
cumulative = TRUE, cuts = 11, emp = TRUE, SIMPLIFY = FALSE)

for (i in seq_len(length(lys_i_cirps))) {
  lys_i_cirps[[i]]$spl <- factor(i)
}
indiv_cirps <- do.call("rbind", lys_i_cirps)
indiv_cirps$mod <- "ci"


###### Combine models

ind_nll <- bindr(pat1 = "indiv", pat2 = "nll")
ind_rps <- bindr(pat1 = "indiv", pat2 = "rps")


# Combine all calibrations of single splits -------------------------------

spl_nll <- rbind(cal_spl_nll, cal_avg_spl_nll)
spl_rps <- rbind(cal_spl_rps, cal_avg_spl_rps)


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


# Combine all average calibrations across splits --------------------------

avg_nll <- rbind(avg_meths_nll, avg_avg_meths_nll)
avg_rps <- rbind(avg_meths_rps, avg_avg_meths_rps)

# Save results ------------------------------------------------------------

write.csv(spl_nll, file = paste0(out_dir, "cumcal_splnll.csv"))
write.csv(spl_rps, file = paste0(out_dir, "cumcal_splrps.csv"))

write.csv(avg_nll, file = paste0(out_dir, "cumcal_avgnll.csv"))
write.csv(avg_rps, file = paste0(out_dir, "cumcal_avgrps.csv"))

write.csv(ind_nll, file = paste0(out_dir, "cumcal_indnll.csv"))
write.csv(ind_rps, file = paste0(out_dir, "cumcal_indrps.csv"))
