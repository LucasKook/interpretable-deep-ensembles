# Performance measures for non-equally weighted ensembles (SI-CS, CI, SI, SI-LS)
# fitted to melanoma data
# Andrea Goetschi
# April 2022

# Dependencies ------------------------------------------------------------

library(etram)

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
fname_sils <- "mela_sils"

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

## Y true

y_true_all <- load_y_true_all(all_y, K = K, t = "test")
y_true_val_all <- load_y_true_all(all_y, K = K, t = "val")

# Evaluation --------------------------------------------------------------

## CILS ##########################################################

### NLL

# weights per ensemble method (for all splits)
w_cilsnll_l <- t(mapply(get_w, lys_cdf_val = cdfval_cilsnll, y_true_val = y_true_val_all, type = "linear",
                        optim_metric = "nll"))
w_cilsnll_ll <- t(mapply(get_w, lys_cdf_val = cdfval_cilsnll, y_true_val = y_true_val_all, type = "log-linear",
                         optim_metric = "nll"))
w_cilsnll_trf <- t(mapply(get_w, lys_cdf_val = cdfval_cilsnll, y_true_val = y_true_val_all, type = "trafo",
                          optim_metric = "nll"))

# metrics per method (for all splits)
met_cilsnll_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "avg",
                                       metrics = "all", topk = FALSE, weights = w_cilsnll_l,
                                       cutoff = 1) # linear weights
met_cilsnll_avgl$method <- "avgl"
met_cilsnll_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "avg",
                                        metrics = "all", topk = FALSE, weights = w_cilsnll_ll,
                                        cutoff = 1) # log-linear weights
met_cilsnll_avgll$method <- "avgll"
met_cilsnll_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "avg",
                                         metrics = "all", topk = FALSE, weights = w_cilsnll_trf,
                                         cutoff = 1) # trafo weights
met_cilsnll_avgtrf$method <- "avgtrf"
met_cilsnll_l <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "linear",
                                    metrics = "all", topk = FALSE, weights = w_cilsnll_l, cutoff = 1)
met_cilsnll_ll <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "log-linear",
                                     metrics = "all", topk = FALSE, weights = w_cilsnll_ll, cutoff = 1)
met_cilsnll_trf <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "trafo",
                                      metrics = "all", topk = FALSE, weights = w_cilsnll_trf, cutoff = 1)
# combine results
met_cilsnll <- rbind(met_cilsnll_avgl, met_cilsnll_avgll, met_cilsnll_avgtrf,
                     met_cilsnll_l, met_cilsnll_ll, met_cilsnll_trf)

### RPS

# weights per ensemble method (for all splits)
w_cilsrps_l <- t(mapply(get_w, lys_cdf_val = cdfval_cilsrps, y_true_val = y_true_val_all, type = "linear",
                        optim_metric = "rps"))
w_cilsrps_ll <- t(mapply(get_w, lys_cdf_val = cdfval_cilsrps, y_true_val = y_true_val_all, type = "log-linear",
                         optim_metric = "rps"))
w_cilsrps_trf <- t(mapply(get_w, lys_cdf_val = cdfval_cilsrps, y_true_val = y_true_val_all, type = "trafo",
                          optim_metric = "rps"))

# metrics per method (for all splits)
met_cilsrps_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "avg",
                                       metrics = "all", topk = FALSE, weights = w_cilsrps_l,
                                       cutoff = 1) # linear weights
met_cilsrps_avgl$method <- "avgl"
met_cilsrps_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "avg",
                                        metrics = "all", topk = FALSE, weights = w_cilsrps_ll,
                                        cutoff = 1) # log-linear weights
met_cilsrps_avgll$method <- "avgll"
met_cilsrps_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "avg",
                                         metrics = "all", topk = FALSE, weights = w_cilsrps_trf,
                                         cutoff = 1) # trafo weights
met_cilsrps_avgtrf$method <- "avgtrf"
met_cilsrps_l <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "linear",
                                    metrics = "all", topk = FALSE, weights = w_cilsrps_l, cutoff = 1)
met_cilsrps_ll <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "log-linear",
                                     metrics = "all", topk = FALSE, weights = w_cilsrps_ll, cutoff = 1)
met_cilsrps_trf <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "trafo",
                                      metrics = "all", topk = FALSE, weights = w_cilsrps_trf, cutoff = 1)
# combine results
met_cilsrps <- rbind(met_cilsrps_avgl, met_cilsrps_avgll, met_cilsrps_avgtrf,
                     met_cilsrps_l, met_cilsrps_ll, met_cilsrps_trf)
# save results
write.csv(met_cilsnll, file = paste0(out_dir, "met_", fname_cilsnll, "_weighted.csv"), row.names = FALSE)
write.csv(met_cilsrps, file = paste0(out_dir, "met_", fname_cilsrps, "_weighted.csv"), row.names = FALSE)

w_cilsnll <- data.frame(w = rbind(w_cilsnll_l,
                                  w_cilsnll_ll,
                                  w_cilsnll_trf),
                        spl = rep(1:splits, 3),
                        method = rep(c("linear", 
                                       "log-linear",
                                       "trafo"), each = splits))

w_cilsrps <- data.frame(w = rbind(w_cilsrps_l,
                                  w_cilsrps_ll,
                                  w_cilsrps_trf),
                        spl = rep(1:splits, 3),
                        method = rep(c("linear", 
                                       "log-linear",
                                       "trafo"), each = splits))

write.csv(w_cilsnll, file = paste0(out_dir, "w_", fname_cilsnll, ".csv"), row.names = FALSE)
write.csv(w_cilsrps, file = paste0(out_dir, "w_", fname_cilsrps, ".csv"), row.names = FALSE)


## CI ############################################################

### NLL

# weights per ensemble method (for all splits)
w_cinll_l <- t(mapply(get_w, lys_cdf_val = cdfval_cinll, y_true_val = y_true_val_all, type = "linear",
                      optim_metric = "nll"))
w_cinll_ll <- t(mapply(get_w, lys_cdf_val = cdfval_cinll, y_true_val = y_true_val_all, type = "log-linear",
                       optim_metric = "nll"))
w_cinll_trf <- t(mapply(get_w, lys_cdf_val = cdfval_cinll, y_true_val = y_true_val_all, type = "trafo",
                        optim_metric = "nll"))

# metrics per method (for all splits)
met_cinll_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "avg",
                                     metrics = "all", topk = FALSE, weights = w_cinll_l,
                                     cutoff = 1) # linear weights
met_cinll_avgl$method <- "avgl"
met_cinll_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "avg",
                                      metrics = "all", topk = FALSE, weights = w_cinll_ll,
                                      cutoff = 1) # log-linear weights
met_cinll_avgll$method <- "avgll"
met_cinll_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "avg",
                                       metrics = "all", topk = FALSE, weights = w_cinll_trf,
                                       cutoff = 1) # trafo weights
met_cinll_avgtrf$method <- "avgtrf"
met_cinll_l <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "linear",
                                  metrics = "all", topk = FALSE, weights = w_cinll_l, cutoff = 1)
met_cinll_ll <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "log-linear",
                                   metrics = "all", topk = FALSE, weights = w_cinll_ll, cutoff = 1)
met_cinll_trf <- get_metrics_allspl(lys_cdf_all = cdftest_cinll, y_true_all = y_true_all, type = "trafo",
                                    metrics = "all", topk = FALSE, weights = w_cinll_trf, cutoff = 1)
# combine results
met_cinll <- rbind(met_cinll_avgl, met_cinll_avgll, met_cinll_avgtrf,
                   met_cinll_l, met_cinll_ll, met_cinll_trf)

### RPS

# weights per ensemble method (for all splits)
w_cirps_l <- t(mapply(get_w, lys_cdf_val = cdfval_cirps, y_true_val = y_true_val_all, type = "linear",
                      optim_metric = "rps"))
w_cirps_ll <- t(mapply(get_w, lys_cdf_val = cdfval_cirps, y_true_val = y_true_val_all, type = "log-linear",
                       optim_metric = "rps"))
w_cirps_trf <- t(mapply(get_w, lys_cdf_val = cdfval_cirps, y_true_val = y_true_val_all, type = "trafo",
                        optim_metric = "rps"))

# metrics per method (for all splits)
met_cirps_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "avg",
                                     metrics = "all", topk = FALSE, weights = w_cirps_l,
                                     cutoff = 1) # linear weights
met_cirps_avgl$method <- "avgl"
met_cirps_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "avg",
                                      metrics = "all", topk = FALSE, weights = w_cirps_ll,
                                      cutoff = 1) # log-linear weights
met_cirps_avgll$method <- "avgll"
met_cirps_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "avg",
                                       metrics = "all", topk = FALSE, weights = w_cirps_trf,
                                       cutoff = 1) # trafo weights
met_cirps_avgtrf$method <- "avgtrf"
met_cirps_l <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "linear",
                                  metrics = "all", topk = FALSE, weights = w_cirps_l, cutoff = 1)
met_cirps_ll <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "log-linear",
                                   metrics = "all", topk = FALSE, weights = w_cirps_ll, cutoff = 1)
met_cirps_trf <- get_metrics_allspl(lys_cdf_all = cdftest_cirps, y_true_all = y_true_all, type = "trafo",
                                    metrics = "all", topk = FALSE, weights = w_cirps_trf, cutoff = 1)
# combine results
met_cirps <- rbind(met_cirps_avgl, met_cirps_avgll, met_cirps_avgtrf,
                   met_cirps_l, met_cirps_ll, met_cirps_trf)
# save results
write.csv(met_cinll, file = paste0(out_dir, "met_", fname_cinll, "_weighted.csv"), row.names = FALSE)
write.csv(met_cirps, file = paste0(out_dir, "met_", fname_cirps, "_weighted.csv"), row.names = FALSE)

w_cinll <- data.frame(w = rbind(w_cinll_l,
                                w_cinll_ll,
                                w_cinll_trf),
                      spl = rep(1:splits, 3),
                      method = rep(c("linear", 
                                     "log-linear",
                                     "trafo"), each = splits))

w_cirps <- data.frame(w = rbind(w_cirps_l,
                                w_cirps_ll,
                                w_cirps_trf),
                      spl = rep(1:splits, 3),
                      method = rep(c("linear", 
                                     "log-linear",
                                     "trafo"), each = splits))

write.csv(w_cinll, file = paste0(out_dir, "w_", fname_cinll, ".csv"), row.names = FALSE)
write.csv(w_cirps, file = paste0(out_dir, "w_", fname_cirps, ".csv"), row.names = FALSE)
