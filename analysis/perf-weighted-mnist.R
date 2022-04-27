# Performance measures for non-equally weighted ensembles (CI) fitted to mnist data
# Andrea Goetschi
# April 2022

# Dependencies ------------------------------------------------------------

library(etram)

# Params ------------------------------------------------------------------

source("experiments/functions/functions_DE.R")

in_dir <- out_dir <- "experiments/results/DE/MNIST/"

splits <- 6
ensembles <- 5

fname_cinll <- "mnist_ci_lossnll_wsno_augno"
fname_cirps <- "mnist_ci_lossrps_wsno_augno"

# Read data ---------------------------------------------------------------

tab_dat <- load_data("mnist")$tab_dat
y <- model.matrix(~ 0 + y, data = tab_dat)

ridx <- get_ridx(in_dir, "mnist")

# Load results ------------------------------------------------------------

cdftest_cinll <- list_cdfs(in_dir, fname_cinll, splits, ensembles, "test")
cdfval_cinll <- list_cdfs(in_dir, fname_cinll, splits, ensembles, "val")

cdftest_cirps <- list_cdfs(in_dir, fname_cirps, splits, ensembles, "test")
cdfval_cirps <- list_cdfs(in_dir, fname_cirps, splits, ensembles, "val")

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

### NLL

# weights per ensemble method (for all splits)
w_cinll_l <- t(mapply(get_w, lys_cdf_val = cdfval_cinll, y_true_val = y_true_val_all, type = "linear",
                      optim_metric = "nll"))
w_cinll_ll <- t(mapply(get_w, lys_cdf_val = cdfval_cinll, y_true_val = y_true_val_all, type = "log-linear",
                       optim_metric = "nll"))
w_cinll_trf <- t(mapply(get_w, lys_cdf_val = cdfval_cinll, y_true_val = y_true_val_all, type = "trafo",
                        optim_metric = "nll"))

w_cinll <- data.frame(w = rbind(w_cinll_l,
                                w_cinll_ll,
                                w_cinll_trf),
                      spl = rep(1:splits, 3),
                      method = rep(c("linear",
                                     "log-linear",
                                     "trafo"), each = splits))

# metrics per method (for all splits)
met_cinll_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "avg",
                                     metrics = "all", topk = FALSE, weights = w_cinll_l) # linear weights
met_cinll_avgl$method <- "avgl"

met_cinll_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "avg",
                                      metrics = "all", topk = FALSE, weights = w_cinll_ll) # log-linear weights
met_cinll_avgll$method <- "avgll"

met_cinll_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "avg",
                                       metrics = "all", topk = FALSE, weights = w_cinll_trf) # trafo weights
met_cinll_avgtrf$method <- "avgtrf"

met_cinll_l <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "linear",
                                  metrics = "all", topk = FALSE, weights = w_cinll_l)

met_cinll_ll <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "log-linear",
                                   metrics = "all", topk = FALSE, weights = w_cinll_ll)
met_cinll_trf <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "trafo",
                                    metrics = "all", topk = FALSE, weights = w_cinll_trf)


met_cinll <- rbind(met_cinll_avgl,
                   met_cinll_avgll,
                   met_cinll_avgtrf,
                   met_cinll_l,
                   met_cinll_ll,
                   met_cinll_trf)

### RPS

# weights per ensemble method (for all splits)
w_cirps_l <- t(mapply(get_w, lys_cdf_val = cdfval_cirps, y_true_val = y_true_val_all, type = "linear",
                      optim_metric = "rps"))
w_cirps_ll <- t(mapply(get_w, lys_cdf_val = cdfval_cirps, y_true_val = y_true_val_all, type = "log-linear",
                       optim_metric = "rps"))
w_cirps_trf <- t(mapply(get_w, lys_cdf_val = cdfval_cirps, y_true_val = y_true_val_all, type = "trafo",
                        optim_metric = "rps"))

w_cirps <- data.frame(w = rbind(w_cirps_l,
                                w_cirps_ll,
                                w_cirps_trf),
                      spl = rep(1:splits, 3),
                      method = rep(c("linear",
                                     "log-linear",
                                     "trafo"), each = splits))

# metrics per method (for all splits)
met_cirps_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "avg",
                                     metrics = "all", topk = FALSE, weights = w_cirps_l) # linear weights
met_cirps_avgl$method <- "avgl"

met_cirps_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "avg",
                                      metrics = "all", topk = FALSE, weights = w_cirps_ll) # log-linear weights
met_cirps_avgll$method <- "avgll"

met_cirps_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "avg",
                                       metrics = "all", topk = FALSE, weights = w_cirps_trf) # trafo weights
met_cirps_avgtrf$method <- "avgtrf"

met_cirps_l <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "linear",
                                  metrics = "all", topk = FALSE, weights = w_cirps_l)
met_cirps_ll <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "log-linear",
                                   metrics = "all", topk = FALSE, weights = w_cirps_ll)
met_cirps_trf <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "trafo",
                                    metrics = "all", topk = FALSE, weights = w_cirps_trf)


met_cirps <- rbind(met_cirps_avgl,
                   met_cirps_avgll,
                   met_cirps_avgtrf,
                   met_cirps_l,
                   met_cirps_ll,
                   met_cirps_trf)

# save results
write.csv(met_cinll, file = paste0(out_dir, "met_", fname_cinll, "_weighted.csv"), row.names = FALSE)
write.csv(met_cirps, file = paste0(out_dir, "met_", fname_cirps, "_weighted.csv"), row.names = FALSE)

write.csv(w_cinll, file = paste0(out_dir, "w_", fname_cinll, ".csv"), row.names = FALSE)
write.csv(w_cirps, file = paste0(out_dir, "w_", fname_cirps, ".csv"), row.names = FALSE)
