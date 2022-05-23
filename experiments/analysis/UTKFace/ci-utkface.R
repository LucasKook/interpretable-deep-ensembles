# Bootstrap confidence intervals
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(etram)
library(boot)

# Directories -------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/UTKFace/"

# Params ------------------------------------------------------------------

fname_silscsnll <- "utkface_silscs_lossnll_wsyes_augno"
fname_silscsrps <- "utkface_silscs_lossrps_wsyes_augno"
fname_cilsnll <- "utkface_cils_lossnll_wsyes_augno"
fname_cilsrps <- "utkface_cils_lossrps_wsyes_augno"
fname_sicsnll <- "utkface_sics_lossnll_wsyes_augno"
fname_sicsrps <- "utkface_sics_lossrps_wsyes_augno"
fname_cinll <- "utkface_ci_lossnll_wsyes_augno"
fname_cirps <- "utkface_ci_lossrps_wsyes_augno"

fname_si <- "utkface_si"
fname_sirps <- "utkface_si_rps"
fname_sils <- "utkface_sils"
fname_silsrps <- "utkface_sils_rps"

nsemi <- 4
K <- 7

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

## Load results of reference model SI-LS
met_sils <-  read.csv(file = paste0(in_dir, "met_", fname_sils, ".csv"))
met_silsrps <- read.csv(file = paste0(in_dir, "met_", fname_silsrps, ".csv"))

# Confidence interval -----------------------------------------------------

args_nll <- data.frame(cdf_all = "all_cdf", y_true_all = "all_y",
                       met_ref = "met_sils",
                       binary = FALSE,
                       mod = c(rep(c("silscs", "cils", "sics", "ci"), each = 10),
                               rep(c("sils", "si"), each = 1)),
                       meth = c(rep(c("linear", "log-linear", "trafo", "avg",
                                      "linear", "log-linear", "trafo",
                                      "avg", "avgll", "avgtrf"), nsemi),
                                rep(NA, 2)),
                       K = K, loss = "nll",
                       weighted = c(rep(c(rep(FALSE, 4), rep(TRUE, 6)), nsemi),
                                    rep(c(FALSE), 2)),
                       avg = c(rep(c(FALSE, FALSE, FALSE, TRUE,
                                     FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), nsemi),
                               rep(FALSE, 2)),
                       fname = c(rep(c(fname_silscsnll, fname_cilsnll,
                                       fname_sicsnll, fname_cinll), each = 10),
                                 rep(c(fname_sils, fname_si), each = 1)),
                       in_dir = in_dir,
                       out_dir = out_dir)

args_rps <- args_nll %>% mutate(met_ref = "met_silsrps",
                                loss = "rps",
                                fname = c(rep(c(fname_silscsrps, fname_cilsrps,
                                                fname_sicsrps, fname_cirps), each = 10),
                                          rep(c(fname_silsrps, fname_sirps), each = 1)))


do.call(Map, c(f = boot_ci, args_nll))
do.call(Map, c(f = boot_ci, args_rps))
