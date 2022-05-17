# Performance measures for equally weighted ensembles (SI-CS, CI, SI, SI-LS)
# fitted to melanoma data
# Andrea Goetschi
# April 2022

# Dependencies ------------------------------------------------------------

library(etram)

# Directories -------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/melanoma/"

# Params -----------------------------------------------------------------

K <- 2

fname_cilsnll <- "mela_silscs_lossnll_wsyes_augno"
fname_cilsrps <- "mela_silscs_lossrps_wsyes_augno"
fname_cinll <- "mela_ci_lossnll_wsyes_augno"
fname_cirps <- "mela_ci_lossrps_wsyes_augno"

fname_si <- "mela_si"
fname_sils <- "mela_sils"
fname_sirps <- "mela_si_rps"
fname_silsrps <- "mela_sils_rps"

# Load results ------------------------------------------------------------

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

met_sirps_all <- get_metrics_allspl(lys_cdf_all = cdftest_sirps, y_true_all = y_true_all, type = "all",
                                    topk = FALSE, metrics = "all", cutoff = 1)

write.csv(met_sirps_all[met_sirps_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_sirps, ".csv"), row.names = FALSE)


## SILS #########################################################

met_sils_all <- get_metrics_allspl(lys_cdf_all = cdftest_sils, y_true_all = y_true_all, type = "all",
                                   topk = FALSE, metrics = "all", cutoff = 1)

write.csv(met_sils_all[met_sils_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_sils, ".csv"), row.names = FALSE)

met_silsrps_all <- get_metrics_allspl(lys_cdf_all = cdftest_silsrps, y_true_all = y_true_all, type = "all",
                                      topk = FALSE, metrics = "all", cutoff = 1)

write.csv(met_silsrps_all[met_silsrps_all$method == "linear", ],
          file = paste0(out_dir, "met_", fname_silsrps, ".csv"), row.names = FALSE)

