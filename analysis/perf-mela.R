# Performance measures for equally weighted ensembles (SI-CS, CI, SI, SI-LS)
# fitted to melanoma data
# Andrea Goetschi
# April 2022

# Dependencies ------------------------------------------------------------

library(readr)
library(rhdf5)
library(ontram)
library(etram)
library(dplyr)

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

fname_si <- "mela_si"
fname_sils <- "mela_sils"

# Read data ---------------------------------------------------------------

dat <- load_data("melanoma", path = path, im_path = im_path)
tab_dat <- dat$tab_dat
y <- model.matrix(~ 0 + target, data = tab_dat)

ridx <- get_ridx(in_dir, fname = "melanoma")

# Load results ------------------------------------------------------------

## CDFs all splits

### CI-LS
cdftest_cilsnll <- list_cdfs(in_dir, fname_cilsnll, splits, ensembles, "test")
cdfval_cilsnll <- list_cdfs(in_dir, fname_cilsnll, splits, ensembles, "val")

cdftest_cilsrps <- list_cdfs(in_dir, fname_cilsrps, splits, ensembles, "test")
cdfval_cilsrps <- list_cdfs(in_dir, fname_cilsrps, splits, ensembles, "val")

### CI
cdftest_cinll <- list_cdfs(in_dir, fname_cinll, splits, ensembles, "test")
cdfval_cinll <- list_cdfs(in_dir, fname_cinll, splits, ensembles, "val")

cdftest_cirps <- list_cdfs(in_dir, fname_cirps, splits, ensembles, "test")
cdfval_cirps <- list_cdfs(in_dir, fname_cirps, splits, ensembles, "val")

### SI

cdftest_si <- list_cdfs(in_dir, fname_si, splits = splits, ensembles = NULL, "test")

### SI-LS

cdftest_sils <- list_cdfs(in_dir, fname_sils, splits = splits, ensembles = NULL, "test")

## Y true

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

## CILS ##########################################################

met_cilsnll_all <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
                                      type = "all", topk = TRUE, order_metric = "nll",
                                      lys_cdf_val_all = cdfval_cilsnll, y_true_val_all = y_true_val_all,
                                      metrics = "all", cutoff = 1)

indi_cilsnll_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
                                             metrics = "all", cutoff = 1)

write.csv(met_cilsnll_all, file = paste0(out_dir, "met_", fname_cilsnll, ".csv"), row.names = FALSE)
write.csv(indi_cilsnll_all, file = paste0(out_dir, "indivmet_", fname_cilsnll, ".csv"), row.names = FALSE)


met_cilsrps_all <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
                                      type = "all", topk = TRUE, order_metric = "rps",
                                      lys_cdf_val_all = cdfval_cilsrps, y_true_val_all = y_true_val_all,
                                      metrics = "all", cutoff = 1)

indi_cilsrps_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
                                             metrics = "all", cutoff = 1)

write.csv(met_cilsrps_all, file = paste0(out_dir, "met_", fname_cilsrps, ".csv"), row.names = FALSE)
write.csv(indi_cilsrps_all, file = paste0(out_dir, "indivmet_", fname_cilsrps, ".csv"), row.names = FALSE)


## CI ############################################################

met_cinll_all <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                                    type = "all", topk = TRUE, order_metric = "nll",
                                    lys_cdf_val_all = cdfval_cinll, y_true_val_all = y_true_val_all,
                                    metrics = "all", cutoff = 1)

indi_cinll_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                                           metrics = "all", cutoff = 1)

write.csv(met_cinll_all, file = paste0(out_dir, "met_", fname_cinll, ".csv"), row.names = FALSE)
write.csv(indi_cinll_all, file = paste0(out_dir, "indivmet_", fname_cinll, ".csv"), row.names = FALSE)


met_cirps_all <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                                    type = "all", topk = TRUE, order_metric = "rps",
                                    lys_cdf_val_all = cdfval_cirps, y_true_val_all = y_true_val_all,
                                    metrics = "all", cutoff = 1)

indi_cirps_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                                           metrics = "all", cutoff = 1)

write.csv(met_cirps_all, file = paste0(out_dir, "met_", fname_cirps, ".csv"), row.names = FALSE)
write.csv(indi_cirps_all, file = paste0(out_dir, "indivmet_", fname_cirps, ".csv"), row.names = FALSE)


## SI ###########################################################

met_si_all <- get_metrics_allspl(lys_cdf_all = cdftest_si, y_true_all = y_true_all, type = "all",
                                 topk = FALSE, metrics = "all", cutoff = 1)

write.csv(met_si_all[met_si_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_si, ".csv"), row.names = FALSE)


## SILS #########################################################

met_sils_all <- get_metrics_allspl(lys_cdf_all = cdftest_sils, y_true_all = y_true_all, type = "all",
                                   topk = FALSE, metrics = "all", cutoff = 1)

write.csv(met_sils_all[met_sils_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_sils, ".csv"), row.names = FALSE)
