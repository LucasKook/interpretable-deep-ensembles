# Performance measures for equally weighted ensembles 
# (SI-CS, CI, SI-LS-CS, CI-LS, SI, SI-LS) fitted to utkface data
# Andrea Goetschi
# April 2022

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

fname_si <- "utkface_si"
fname_sils <- "utkface_sils"


# Read data ---------------------------------------------------------------

dat <- load_data("utkface", path = path)
tab_dat <- dat$tab_dat
y <- model.matrix(~ 0 + age_group, data = tab_dat)

ridx <- get_ridx(in_dir, fname = "utkface")

# Load results ------------------------------------------------------------

## CDFs all splits

### SI-LS-CS
cdftest_silscsnll <- list_cdfs(in_dir, fname_silscsnll, splits, ensembles, "test")
cdfval_silscsnll <- list_cdfs(in_dir, fname_silscsnll, splits, ensembles, "val")

cdftest_silscsrps <- list_cdfs(in_dir, fname_silscsrps, splits, ensembles, "test")
cdfval_silscsrps <- list_cdfs(in_dir, fname_silscsrps, splits, ensembles, "val")

### CI-LS
cdftest_cilsnll <- list_cdfs(in_dir, fname_cilsnll, splits, ensembles, "test")
cdfval_cilsnll <- list_cdfs(in_dir, fname_cilsnll, splits, ensembles, "val")

cdftest_cilsrps <- list_cdfs(in_dir, fname_cilsrps, splits, ensembles, "test")
cdfval_cilsrps <- list_cdfs(in_dir, fname_cilsrps, splits, ensembles, "val")

### SI-CS
cdftest_sicsnll <- list_cdfs(in_dir, fname_sicsnll, splits, ensembles, "test")
cdfval_sicsnll <- list_cdfs(in_dir, fname_sicsnll, splits, ensembles, "val")

cdftest_sicsrps <- list_cdfs(in_dir, fname_sicsrps, splits, ensembles, "test")
cdfval_sicsrps <- list_cdfs(in_dir, fname_sicsrps, splits, ensembles, "val")

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

## SILSCS ########################################################

met_silscsnll_all <- get_metrics_allspl(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, 
                                        type = "all", topk = TRUE, order_metric = "nll",
                                        lys_cdf_val_all = cdfval_silscsnll, y_true_val_all = y_true_val_all,
                                        metrics = "all")

indi_silscsnll_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all,
                                               metrics = "all")

write.csv(met_silscsnll_all, file = paste0(out_dir, "met_", fname_silscsnll, ".csv"), row.names = FALSE)
write.csv(indi_silscsnll_all, file = paste0(out_dir, "indivmet_", fname_silscsnll, ".csv"), row.names = FALSE)


met_silscsrps_all <- get_metrics_allspl(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, 
                                        type = "all", topk = TRUE, order_metric = "rps",
                                        lys_cdf_val_all = cdfval_silscsrps, y_true_val_all = y_true_val_all,
                                        metrics = "all")

indi_silscsrps_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all,
                                               metrics = "all")

write.csv(met_silscsrps_all, file = paste0(out_dir, "met_", fname_silscsrps, ".csv"), row.names = FALSE)
write.csv(indi_silscsrps_all, file = paste0(out_dir, "indivmet_", fname_silscsrps, ".csv"), row.names = FALSE)

## CILS ##########################################################

met_cilsnll_all <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
                                      type = "all", topk = TRUE, order_metric = "nll",
                                      lys_cdf_val_all = cdfval_cilsnll, y_true_val_all = y_true_val_all,
                                      metrics = "all")

indi_cilsnll_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all,
                                             metrics = "all")

write.csv(met_cilsnll_all, file = paste0(out_dir, "met_", fname_cilsnll, ".csv"), row.names = FALSE)
write.csv(indi_cilsnll_all, file = paste0(out_dir, "indivmet_", fname_cilsnll, ".csv"), row.names = FALSE)


met_cilsrps_all <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
                                      type = "all", topk = TRUE, order_metric = "rps",
                                      lys_cdf_val_all = cdfval_cilsrps, y_true_val_all = y_true_val_all,
                                      metrics = "all")

indi_cilsrps_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all,
                                             metrics = "all")

write.csv(met_cilsrps_all, file = paste0(out_dir, "met_", fname_cilsrps, ".csv"), row.names = FALSE)
write.csv(indi_cilsrps_all, file = paste0(out_dir, "indivmet_", fname_cilsrps, ".csv"), row.names = FALSE)


## SICS ##########################################################

met_sicsnll_all <- get_metrics_allspl(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all,
                                      type = "all", topk = TRUE, order_metric = "nll",
                                      lys_cdf_val_all = cdfval_sicsnll, y_true_val_all = y_true_val_all,
                                      metrics = "all")

indi_sicsnll_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all,
                                             metrics = "all")

write.csv(met_sicsnll_all, file = paste0(out_dir, "met_", fname_sicsnll, ".csv"), row.names = FALSE)
write.csv(indi_sicsnll_all, file = paste0(out_dir, "indivmet_", fname_sicsnll, ".csv"), row.names = FALSE)


met_sicsrps_all <- get_metrics_allspl(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all,
                                      type = "all", topk = TRUE, order_metric = "rps",
                                      lys_cdf_val_all = cdfval_sicsrps, y_true_val_all = y_true_val_all,
                                      metrics = "all")

indi_sicsrps_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all,
                                             metrics = "all")

write.csv(met_sicsrps_all, file = paste0(out_dir, "met_", fname_sicsrps, ".csv"), row.names = FALSE)
write.csv(indi_sicsrps_all, file = paste0(out_dir, "indivmet_", fname_sicsrps, ".csv"), row.names = FALSE)

## CI ############################################################

met_cinll_all <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                                    type = "all", topk = TRUE, order_metric = "nll",
                                    lys_cdf_val_all = cdfval_cinll, y_true_val_all = y_true_val_all,
                                    metrics = "all")

indi_cinll_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all,
                                           metrics = "all")

write.csv(met_cinll_all, file = paste0(out_dir, "met_", fname_cinll, ".csv"), row.names = FALSE)
write.csv(indi_cinll_all, file = paste0(out_dir, "indivmet_", fname_cinll, ".csv"), row.names = FALSE)


met_cirps_all <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                                    type = "all", topk = TRUE, order_metric = "rps",
                                    lys_cdf_val_all = cdfval_cirps, y_true_val_all = y_true_val_all,
                                    metrics = "all")

indi_cirps_all <- get_indiv_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all,
                                           metrics = "all")

write.csv(met_cirps_all, file = paste0(out_dir, "met_", fname_cirps, ".csv"), row.names = FALSE)
write.csv(indi_cirps_all, file = paste0(out_dir, "indivmet_", fname_cirps, ".csv"), row.names = FALSE)


## SI ###########################################################

met_si_all <- get_metrics_allspl(lys_cdf_all = cdftest_si, y_true_all = y_true_all, type = "all",
                                 topk = FALSE, metrics = "all")

write.csv(met_si_all[met_si_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_si, ".csv"), row.names = FALSE)


## SILS #########################################################

met_sils_all <- get_metrics_allspl(lys_cdf_all = cdftest_sils, y_true_all = y_true_all, type = "all",
                                   topk = FALSE, metrics = "all")

write.csv(met_sils_all[met_sils_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_sils, ".csv"), row.names = FALSE)

