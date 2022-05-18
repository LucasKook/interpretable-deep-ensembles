# Calibration mnist data
# Andrea Goetschi
# April 2022

# Deps --------------------------------------------------------------------

library(etram)
library(caret)
library(Hmisc)

# Params ------------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/MNIST/"

fname_cinll <- "mnist_ci_lossnll_wsno_augno"
fname_cirps <- "mnist_ci_lossrps_wsno_augno"

nsemi <- 1
K <- 10

# Functions ---------------------------------------------------------------

cutf <- function(x) {
  q1 <- quantile(x, 0.05)
  q2 <- quantile(x, 0.95)
  cstep <- (q2 - q1) / 9
  cuts <- seq(from = q1, to = q2, by = cstep)
  return(cuts)
}

# Load results ------------------------------------------------------------

## all CDF, Y
cdf_nll <- read.csv(paste0(in_dir, "mnist_merged_cdf_cinll.csv"))
cdf_rps <- read.csv(paste0(in_dir, "mnist_merged_cdf_cirps.csv"))
all_cdf <- rbind(cdf_nll, cdf_rps)

all_y <- read.csv(paste0(in_dir, "mnist_merged_y.csv"))

# Calibration -------------------------------------------------------------

args_nll <- data.frame(cdf_all = "all_cdf", y_true_all = "all_y",
                       mod = rep("ci", each = 10),
                       meth = rep(c("linear", "log-linear", "trafo", "avg",
                                    "linear", "log-linear", "trafo",
                                    "avg", "avgll", "avgtrf"), nsemi),
                       K = K, loss = "nll",
                       weighted = rep(c(rep(FALSE, 4), rep(TRUE, 6)), nsemi),
                       avg = rep(c(FALSE, FALSE, FALSE, TRUE,
                                   FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), nsemi),
                       cuts_fun = "cutf",
                       fname = rep(fname_cinll, each = 10),
                       in_dir = in_dir)

args_rps <- args_nll %>% mutate(loss = "rps",
                                fname = rep(fname_cirps, each = 10))

res_nll <- do.call(Map, c(f = calc_calibration, args_nll))
res_nll <- bind_rows(res_nll)

res_rps <- do.call(Map, c(f = calc_calibration, args_rps))
res_rps <- bind_rows(res_rps)

# Save results ------------------------------------------------------------

write.csv(res_nll, file = paste0(out_dir, "cumcal_avgnll.csv"))
write.csv(res_rps, file = paste0(out_dir, "cumcal_avgrps.csv"))
