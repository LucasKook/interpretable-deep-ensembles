# Bootstrap confidence intervals
# Andrea Goetschi
# May 2022

# Dependencies ------------------------------------------------------------

library(etram)
library(boot)

# Directories -------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/MNIST/"

# Params ------------------------------------------------------------------

fname_cinll <- "mnist_ci_lossnll_wsno_augno"
fname_cirps <- "mnist_ci_lossrps_wsno_augno"

nsemi <- 1
K <- 10

# Load results ------------------------------------------------------------

## all CDF, Y
cdf_nll <- read.csv(paste0(in_dir, "mnist_merged_cdf_cinll.csv"))
cdf_rps <- read.csv(paste0(in_dir, "mnist_merged_cdf_cirps.csv"))
all_cdf <- rbind(cdf_nll, cdf_rps)

all_y <- read.csv(paste0(in_dir, "mnist_merged_y.csv"))

# Confidence interval -----------------------------------------------------

args_nll <- data.frame(cdf_all = "all_cdf", y_true_all = "all_y",
                       met_ref = "NULL",
                       binary = FALSE,
                       mod = rep("ci", each = 10),
                       meth = rep(c("linear", "log-linear", "trafo", "avg",
                                    "linear", "log-linear", "trafo",
                                    "avg", "avgll", "avgtrf"), nsemi),
                       K = K, loss = "nll",
                       weighted = rep(c(rep(FALSE, 4), rep(TRUE, 6)), nsemi),
                       avg = rep(c(FALSE, FALSE, FALSE, TRUE,
                                   FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), nsemi),
                       fname = rep(fname_cinll, each = 10),
                       in_dir = in_dir,
                       out_dir = out_dir)

args_rps <- args_nll %>% mutate(loss = "rps",
                                fname = rep(fname_cirps, each = 10))

do.call(Map, c(f = boot_ci, args_nll))
do.call(Map, c(f = boot_ci, args_rps))
