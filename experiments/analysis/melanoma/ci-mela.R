# Bootstrap confidence intervals
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(etram)
library(boot)

# Directories -------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/melanoma/"

# Params ------------------------------------------------------------------

fname_cilsnll <- "mela_silscs_lossnll_wsyes_augno"
fname_cilsrps <- "mela_silscs_lossrps_wsyes_augno"
fname_cinll <- "mela_ci_lossnll_wsyes_augno"
fname_cirps <- "mela_ci_lossrps_wsyes_augno"

fname_si <- "mela_si"
fname_sirps <- "mela_si_rps"
fname_sils <- "mela_sils"
fname_silsrps <- "mela_sils_rps"

nsemi <- 2
K <- 2

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

## Load results of reference model SI-LS
met_sils <-  read.csv(file = paste0(in_dir, "met_", fname_sils, ".csv"))
met_silsrps <- read.csv(file = paste0(in_dir, "met_", fname_silsrps, ".csv"))

# Confidence interval -----------------------------------------------------

args_nll <- data.frame(cdf_all = "all_cdf", y_true_all = "all_y",
                       met_ref = "met_sils",
                       binary = TRUE,
                       mod = c(rep(c("cils", "ci"), each = 10),
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
                       fname = c(rep(c(fname_cilsnll, fname_cinll), each = 10),
                                 rep(c(fname_sils, fname_si), each = 1)),
                       in_dir = in_dir,
                       out_dir = out_dir)

args_rps <- args_nll %>% mutate(met_ref = "met_silsrps",
                                loss = "rps",
                                fname = c(rep(c(fname_cilsrps, fname_cirps), each = 10),
                                          rep(c(fname_silsrps, fname_sirps), each = 1)))


do.call(Map, c(f = boot_ci, args_nll))
do.call(Map, c(f = boot_ci, args_rps))
