# Performance measures for equally weighted ensembles 
# (SI-CS, CI, SI-LS-CS, CI-LS, SI, SI-LS) fitted to utkface data
# Andrea Goetschi
# April 2022

# Dependencies ------------------------------------------------------------

library(etram)

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
fname_sils <- "utkface_sils"
fname_sirps <- "utkface_si_rps"
fname_silsrps <- "utkface_sils_rps"

# Load results ------------------------------------------------------------

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
cdfval_silscsnll <- load_lys_cdf_all(all_cdf, m = "silscs", K = K, l = "nll", t = "val")

cdftest_silscsrps <- load_lys_cdf_all(all_cdf, m = "silscs", K = K, l = "rps", t = "test")
cdfval_silscsrps <- load_lys_cdf_all(all_cdf, m = "silscs", K = K, l = "rps", t = "val")

### SI-CS
cdftest_sicsnll <- load_lys_cdf_all(all_cdf, m = "sics", K = K, l = "nll", t = "test")
cdfval_sicsnll <- load_lys_cdf_all(all_cdf, m = "sics", K = K, l = "nll", t = "val")

cdftest_sicsrps <- load_lys_cdf_all(all_cdf, m = "sics", K = K, l = "rps", t = "test")
cdfval_sicsrps <- load_lys_cdf_all(all_cdf, m = "sics", K = K, l = "rps", t = "val")

### CI-LS
cdftest_cilsnll <- load_lys_cdf_all(all_cdf, m = "cils", K = K, l = "nll", t = "test")
cdfval_cilsnll <- load_lys_cdf_all(all_cdf, m = "cils", K = K, l = "nll", t = "val")

cdftest_cilsrps <- load_lys_cdf_all(all_cdf, m = "cils", K = K, l = "rps", t = "test")
cdfval_cilsrps <- load_lys_cdf_all(all_cdf, m = "cils", K = K, l = "rps", t = "val")

### CI
cdftest_cinll <- load_lys_cdf_all(all_cdf, m = "ci", K = K, l = "nll", t = "test")
cdfval_cinll <- load_lys_cdf_all(all_cdf, m = "ci", K = K, l = "nll", t = "val")

cdftest_cirps <- load_lys_cdf_all(all_cdf, m = "ci", K = K, l = "rps", t = "test")
cdfval_cirps <- load_lys_cdf_all(all_cdf, m = "ci", K = K, l = "rps", t = "val")

### SI
cdftest_si <- load_lys_cdf_all(all_cdf, m = "si", K = K, l = "nll", t = "test")
cdftest_sirps <- load_lys_cdf_all(all_cdf, m = "si", K = K, l = "rps", t = "test")

### SI-LS
cdftest_sils <- load_lys_cdf_all(all_cdf, m = "sils", K = K, l = "nll", t = "test")
cdftest_silsrps <- load_lys_cdf_all(all_cdf, m = "sils", K = K, l = "rps", t = "test")

## Y true

y_true_all <- load_y_true_all(all_y, K = K, t = "test")
y_true_val_all <- load_y_true_all(all_y, K = K, t = "val")

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

met_sirps_all <- get_metrics_allspl(lys_cdf_all = cdftest_sirps, y_true_all = y_true_all, type = "all",
                                    topk = FALSE, metrics = "all")

write.csv(met_sirps_all[met_sirps_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_sirps, ".csv"), row.names = FALSE)


## SILS #########################################################

met_sils_all <- get_metrics_allspl(lys_cdf_all = cdftest_sils, y_true_all = y_true_all, type = "all",
                                   topk = FALSE, metrics = "all")

write.csv(met_sils_all[met_sils_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_sils, ".csv"), row.names = FALSE)

met_silsrps_all <- get_metrics_allspl(lys_cdf_all = cdftest_silsrps, y_true_all = y_true_all, type = "all",
                                      topk = FALSE, metrics = "all")

write.csv(met_silsrps_all[met_silsrps_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_silsrps, ".csv"), row.names = FALSE)

