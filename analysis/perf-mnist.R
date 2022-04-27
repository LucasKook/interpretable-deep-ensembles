# Performance measures for equally weighted ensembles (CI) fitted to mnist data
# Andrea Goetschi
# April 2022

# Dependencies ------------------------------------------------------------

library(etram)

# Params ------------------------------------------------------------------

source("experiments/functions/functions_DE.R")

in_dir <- out_dir <- "experiments/results/DE/MNIST/"

fname_cinll <- "mnist_ci_lossnll_wsno_augno"
fname_cirps <- "mnist_ci_lossrps_wsno_augno"

splits <- 6
ensembles <- 5

# Read data ---------------------------------------------------------------

dat <- dat <- load_data("mnist")
tab_dat <- dat$tab_dat
y <- model.matrix(~ 0 + y, data = tab_dat)

ridx <- get_ridx(in_dir, fname = "mnist")

# Load results ------------------------------------------------------------

## CDFs all splits

### CI
cdftest_cinll <- list_cdfs(in_dir, fname_cinll, splits, ensembles, "test")
cdfval_cinll <- list_cdfs(in_dir, fname_cinll, splits, ensembles, "val")
cdftest_cirps <- list_cdfs(in_dir, fname_cirps, splits, ensembles, "test")
cdfval_cirps <- list_cdfs(in_dir, fname_cirps, splits, ensembles, "val")

### Y true

ytest_1 <- y[ridx[ridx$spl == 1 & ridx$type == "test", "idx"], ]
ytest_2 <- y[ridx[ridx$spl == 2 & ridx$type == "test", "idx"], ]
ytest_3 <- y[ridx[ridx$spl == 3 & ridx$type == "test", "idx"], ]
ytest_4 <- y[ridx[ridx$spl == 4 & ridx$type == "test", "idx"], ]
ytest_5 <- y[ridx[ridx$spl == 5 & ridx$type == "test", "idx"], ]
ytest_6 <- y[ridx[ridx$spl == 6 & ridx$type == "test", "idx"], ]
y_true_all <- list(ytest_1, ytest_2, ytest_3, ytest_4, ytest_5, ytest_6)

yval_1 <- y[ridx[ridx$spl == 1 & ridx$type == "val", "idx"], ]
yval_2 <- y[ridx[ridx$spl == 2 & ridx$type == "val", "idx"], ]
yval_3 <- y[ridx[ridx$spl == 3 & ridx$type == "val", "idx"], ]
yval_4 <- y[ridx[ridx$spl == 4 & ridx$type == "val", "idx"], ]
yval_5 <- y[ridx[ridx$spl == 5 & ridx$type == "val", "idx"], ]
yval_6 <- y[ridx[ridx$spl == 6 & ridx$type == "val", "idx"], ]
y_true_val_all <- list(yval_1, yval_2, yval_3, yval_4, yval_5, yval_6)


# Evaluation --------------------------------------------------------------

## CI ############################################################

## NLL LOSS

met_cinll_all <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                                    type = "all", topk = TRUE, order_metric = "nll",
                                    lys_cdf_val_all = cdfval_cinll, y_true_val =  y_true_val_all,
                                    metrics = "all")

indi_cinll_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                                           metrics = "all")

write.csv(met_cinll_all, file = paste0(out_dir, "met_", fname_cinll, ".csv"), row.names = FALSE)
write.csv(indi_cinll_all, file = paste0(out_dir, "indivmet_", fname_cinll, ".csv"), row.names = FALSE)


## RPS LOSS

met_cirps_all <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                                    type = "all", topk = TRUE, order_metric = "rps",
                                    lys_cdf_val_all = cdfval_cirps, y_true_val =  y_true_val_all,
                                    metrics = "all")

indi_cirps_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                                           metrics = "all")

write.csv(met_cirps_all, file = paste0(out_dir, "met_", fname_cirps, ".csv"), row.names = FALSE)
write.csv(indi_cirps_all, file = paste0(out_dir, "indivmet_", fname_cirps, ".csv"), row.names = FALSE)
