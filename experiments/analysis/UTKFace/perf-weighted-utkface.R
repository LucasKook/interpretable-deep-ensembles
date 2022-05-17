# Performance measures for non-equally weighted ensembles 
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

## Y true

y_true_all <- load_y_true_all(all_y, K = K, t = "test")
y_true_val_all <- load_y_true_all(all_y, K = K, t = "val")

# Evaluation --------------------------------------------------------------

## SILSCS ########################################################

### NLL

# weights per ensemble method (for all splits)
w_silscsnll_l <- t(mapply(get_w, lys_cdf_val = cdfval_silscsnll, y_true_val = y_true_val_all, type = "linear",
                          optim_metric = "nll"))
w_silscsnll_ll <- t(mapply(get_w, lys_cdf_val = cdfval_silscsnll, y_true_val = y_true_val_all, type = "log-linear",
                           optim_metric = "nll"))
w_silscsnll_trf <- t(mapply(get_w, lys_cdf_val = cdfval_silscsnll, y_true_val = y_true_val_all, type = "trafo",
                            optim_metric = "nll"))

# metrics per method (for all splits)
met_silscsnll_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, type = "avg",
                                         metrics = "all", topk = FALSE, weights = w_silscsnll_l) # linear weights
met_silscsnll_avgl$method <- "avgl"
met_silscsnll_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, type = "avg",
                                          metrics = "all", topk = FALSE, weights = w_silscsnll_ll) # log-linear weights
met_silscsnll_avgll$method <- "avgll"
met_silscsnll_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, type = "avg",
                                           metrics = "all", topk = FALSE, weights = w_silscsnll_trf) # trafo weights
met_silscsnll_avgtrf$method <- "avgtrf"
met_silscsnll_l <- get_metrics_allspl(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, type = "linear",
                                      metrics = "all", topk = FALSE, weights = w_silscsnll_l)
met_silscsnll_ll <- get_metrics_allspl(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, type = "log-linear",
                                       metrics = "all", topk = FALSE, weights = w_silscsnll_ll)
met_silscsnll_trf <- get_metrics_allspl(lys_cdf_all = cdftest_silscsnll, y_true_all = y_true_all, type = "trafo",
                                        metrics = "all", topk = FALSE, weights = w_silscsnll_trf)
# combine results
met_silscsnll <- rbind(met_silscsnll_avgl, met_silscsnll_avgll, met_silscsnll_avgtrf,
                       met_silscsnll_l, met_silscsnll_ll, met_silscsnll_trf)

### RPS

# weights per ensemble method (for all splits)
w_silscsrps_l <- t(mapply(get_w, lys_cdf_val = cdfval_silscsrps, y_true_val = y_true_val_all, type = "linear",
                          optim_metric = "rps"))
w_silscsrps_ll <- t(mapply(get_w, lys_cdf_val = cdfval_silscsrps, y_true_val = y_true_val_all, type = "log-linear",
                           optim_metric = "rps"))
w_silscsrps_trf <- t(mapply(get_w, lys_cdf_val = cdfval_silscsrps, y_true_val = y_true_val_all, type = "trafo",
                            optim_metric = "rps"))

# metrics per method (for all splits)
met_silscsrps_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, type = "avg",
                                         metrics = "all", topk = FALSE, weights = w_silscsrps_l) # linear weights
met_silscsrps_avgl$method <- "avgl"
met_silscsrps_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, type = "avg",
                                          metrics = "all", topk = FALSE, weights = w_silscsrps_ll) # log-linear weights
met_silscsrps_avgll$method <- "avgll"
met_silscsrps_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, type = "avg",
                                           metrics = "all", topk = FALSE, weights = w_silscsrps_trf) # trafo weights
met_silscsrps_avgtrf$method <- "avgtrf"
met_silscsrps_l <- get_metrics_allspl(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, type = "linear",
                                      metrics = "all", topk = FALSE, weights = w_silscsrps_l)
met_silscsrps_ll <- get_metrics_allspl(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, type = "log-linear",
                                       metrics = "all", topk = FALSE, weights = w_silscsrps_ll)
met_silscsrps_trf <- get_metrics_allspl(lys_cdf_all = cdftest_silscsrps, y_true_all = y_true_all, type = "trafo",
                                        metrics = "all", topk = FALSE, weights = w_silscsrps_trf)
# combine results
met_silscsrps <- rbind(met_silscsrps_avgl, met_silscsrps_avgll, met_silscsrps_avgtrf,
                       met_silscsrps_l, met_silscsrps_ll, met_silscsrps_trf)
# save results
write.csv(met_silscsnll, file = paste0(out_dir, "met_", fname_silscsnll, "_weighted.csv"), row.names = FALSE)
write.csv(met_silscsrps, file = paste0(out_dir, "met_", fname_silscsrps, "_weighted.csv"), row.names = FALSE)

w_silscsnll <- data.frame(w = rbind(w_silscsnll_l,
                                    w_silscsnll_ll,
                                    w_silscsnll_trf),
                          spl = rep(1:splits, 3),
                          method = rep(c("linear", 
                                         "log-linear",
                                         "trafo"), each = splits))

w_silscsrps <- data.frame(w = rbind(w_silscsrps_l,
                                    w_silscsrps_ll,
                                    w_silscsrps_trf),
                          spl = rep(1:splits, 3),
                          method = rep(c("linear", 
                                         "log-linear",
                                         "trafo"), each = splits))

write.csv(w_silscsnll, file = paste0(out_dir, "w_", fname_silscsnll, ".csv"), row.names = FALSE)
write.csv(w_silscsrps, file = paste0(out_dir, "w_", fname_silscsrps, ".csv"), row.names = FALSE)


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
                                       metrics = "all", topk = FALSE, weights = w_cilsnll_l) # linear weights
met_cilsnll_avgl$method <- "avgl"
met_cilsnll_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "avg",
                                        metrics = "all", topk = FALSE, weights = w_cilsnll_ll) # log-linear weights
met_cilsnll_avgll$method <- "avgll"
met_cilsnll_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "avg",
                                         metrics = "all", topk = FALSE, weights = w_cilsnll_trf) # trafo weights
met_cilsnll_avgtrf$method <- "avgtrf"
met_cilsnll_l <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "linear",
                                    metrics = "all", topk = FALSE, weights = w_cilsnll_l)
met_cilsnll_ll <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "log-linear",
                                     metrics = "all", topk = FALSE, weights = w_cilsnll_ll)
met_cilsnll_trf <- get_metrics_allspl(lys_cdf_all = cdftest_cilsnll, y_true_all = y_true_all, type = "trafo",
                                      metrics = "all", topk = FALSE, weights = w_cilsnll_trf)
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
                                       metrics = "all", topk = FALSE, weights = w_cilsrps_l) # linear weights
met_cilsrps_avgl$method <- "avgl"
met_cilsrps_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "avg",
                                        metrics = "all", topk = FALSE, weights = w_cilsrps_ll) # log-linear weights
met_cilsrps_avgll$method <- "avgll"
met_cilsrps_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "avg",
                                         metrics = "all", topk = FALSE, weights = w_cilsrps_trf) # trafo weights
met_cilsrps_avgtrf$method <- "avgtrf"
met_cilsrps_l <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "linear",
                                    metrics = "all", topk = FALSE, weights = w_cilsrps_l)
met_cilsrps_ll <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "log-linear",
                                     metrics = "all", topk = FALSE, weights = w_cilsrps_ll)
met_cilsrps_trf <- get_metrics_allspl(lys_cdf_all = cdftest_cilsrps, y_true_all = y_true_all, type = "trafo",
                                      metrics = "all", topk = FALSE, weights = w_cilsrps_trf)
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

## SICS ##########################################################

### NLL

# weights per ensemble method (for all splits)
w_sicsnll_l <- t(mapply(get_w, lys_cdf_val = cdfval_sicsnll, y_true_val = y_true_val_all, type = "linear",
                        optim_metric = "nll"))
w_sicsnll_ll <- t(mapply(get_w, lys_cdf_val = cdfval_sicsnll, y_true_val = y_true_val_all, type = "log-linear",
                         optim_metric = "nll"))
w_sicsnll_trf <- t(mapply(get_w, lys_cdf_val = cdfval_sicsnll, y_true_val = y_true_val_all, type = "trafo",
                          optim_metric = "nll"))

# metrics per method (for all splits)
met_sicsnll_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, type = "avg",
                                       metrics = "all", topk = FALSE, weights = w_sicsnll_l) # linear weights
met_sicsnll_avgl$method <- "avgl"
met_sicsnll_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, type = "avg",
                                        metrics = "all", topk = FALSE, weights = w_sicsnll_ll) # log-linear weights
met_sicsnll_avgll$method <- "avgll"
met_sicsnll_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, type = "avg",
                                         metrics = "all", topk = FALSE, weights = w_sicsnll_trf) # trafo weights
met_sicsnll_avgtrf$method <- "avgtrf"
met_sicsnll_l <- get_metrics_allspl(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, type = "linear",
                                    metrics = "all", topk = FALSE, weights = w_sicsnll_l)
met_sicsnll_ll <- get_metrics_allspl(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, type = "log-linear",
                                     metrics = "all", topk = FALSE, weights = w_sicsnll_ll)
met_sicsnll_trf <- get_metrics_allspl(lys_cdf_all = cdftest_sicsnll, y_true_all = y_true_all, type = "trafo",
                                      metrics = "all", topk = FALSE, weights = w_sicsnll_trf)
# combine results
met_sicsnll <- rbind(met_sicsnll_avgl, met_sicsnll_avgll, met_sicsnll_avgtrf,
                     met_sicsnll_l, met_sicsnll_ll, met_sicsnll_trf)

### RPS

# weights per ensemble method (for all splits)
w_sicsrps_l <- t(mapply(get_w, lys_cdf_val = cdfval_sicsrps, y_true_val = y_true_val_all, type = "linear",
                        optim_metric = "rps"))
w_sicsrps_ll <- t(mapply(get_w, lys_cdf_val = cdfval_sicsrps, y_true_val = y_true_val_all, type = "log-linear",
                         optim_metric = "rps"))
w_sicsrps_trf <- t(mapply(get_w, lys_cdf_val = cdfval_sicsrps, y_true_val = y_true_val_all, type = "trafo",
                          optim_metric = "rps"))

# metrics per method (for all splits)
met_sicsrps_avgl <- get_metrics_allspl(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, type = "avg",
                                       metrics = "all", topk = FALSE, weights = w_sicsrps_l) # linear weights
met_sicsrps_avgl$method <- "avgl"
met_sicsrps_avgll <- get_metrics_allspl(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, type = "avg",
                                        metrics = "all", topk = FALSE, weights = w_sicsrps_ll) # log-linear weights
met_sicsrps_avgll$method <- "avgll"
met_sicsrps_avgtrf <- get_metrics_allspl(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, type = "avg",
                                         metrics = "all", topk = FALSE, weights = w_sicsrps_trf) # trafo weights
met_sicsrps_avgtrf$method <- "avgtrf"
met_sicsrps_l <- get_metrics_allspl(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, type = "linear",
                                    metrics = "all", topk = FALSE, weights = w_sicsrps_l)
met_sicsrps_ll <- get_metrics_allspl(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, type = "log-linear",
                                     metrics = "all", topk = FALSE, weights = w_sicsrps_ll)
met_sicsrps_trf <- get_metrics_allspl(lys_cdf_all = cdftest_sicsrps, y_true_all = y_true_all, type = "trafo",
                                      metrics = "all", topk = FALSE, weights = w_sicsrps_trf)
# combine results
met_sicsrps <- rbind(met_sicsrps_avgl, met_sicsrps_avgll, met_sicsrps_avgtrf,
                     met_sicsrps_l, met_sicsrps_ll, met_sicsrps_trf)
# save results
write.csv(met_sicsnll, file = paste0(out_dir, "met_", fname_sicsnll, "_weighted.csv"), row.names = FALSE)
write.csv(met_sicsrps, file = paste0(out_dir, "met_", fname_sicsrps, "_weighted.csv"), row.names = FALSE)

w_sicsnll <- data.frame(w = rbind(w_sicsnll_l,
                                  w_sicsnll_ll,
                                  w_sicsnll_trf),
                        spl = rep(1:splits, 3),
                        method = rep(c("linear", 
                                       "log-linear",
                                       "trafo"), each = splits))

w_sicsrps <- data.frame(w = rbind(w_sicsrps_l,
                                  w_sicsrps_ll,
                                  w_sicsrps_trf),
                        spl = rep(1:splits, 3),
                        method = rep(c("linear", 
                                       "log-linear",
                                       "trafo"), each = splits))

write.csv(w_sicsnll, file = paste0(out_dir, "w_", fname_sicsnll, ".csv"), row.names = FALSE)
write.csv(w_sicsrps, file = paste0(out_dir, "w_", fname_sicsrps, ".csv"), row.names = FALSE)


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
