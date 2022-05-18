# Calibration melanoma data
# Andrea Goetschi
# April 2022

# Deps --------------------------------------------------------------------

library(etram)
library(caret)
library(Hmisc)

# Params ------------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/melanoma/"

fname_cilsnll <- "mela_silscs_lossnll_wsyes_augno"
fname_cilsrps <- "mela_silscs_lossrps_wsyes_augno"
fname_cinll <- "mela_ci_lossnll_wsyes_augno"
fname_cirps <- "mela_ci_lossrps_wsyes_augno"
fname_sils <- "mela_sils"
fname_silsrps <- "mela_sils_rps"

nsemi <- 2
K <- 2

# Functions ---------------------------------------------------------------

cutf <- function(x) quantile(x, c(0.5, 0.9, 0.99, 0.999))

# cut points for sils model
cutf_sils <- function(x) quantile(x, c(0.5, 0.9, 0.99))

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

# Calibration -------------------------------------------------------------

args_nll <- data.frame(cdf_all = "all_cdf", y_true_all = "all_y",
                       mod = c(rep(c("cils", "ci"), each = 10),
                               rep(c("sils"), each = 2)),
                       meth = c(rep(c("linear", "log-linear", "trafo", "avg",
                                      "linear", "log-linear", "trafo",
                                      "avg", "avgll", "avgtrf"), nsemi),
                                rep(NA, 2)),
                       K = K, loss = "nll",
                       weighted = c(rep(c(rep(FALSE, 4), rep(TRUE, 6)), nsemi),
                                    rep(c(FALSE, TRUE), each = 1)),
                       avg = c(rep(c(FALSE, FALSE, FALSE, TRUE,
                                     FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), nsemi),
                               rep(FALSE, 2)),
                       cuts_fun = "cutf",
                       fname = c(rep(c(fname_cilsnll, fname_cinll), each = 10),
                                 rep(c(fname_sils), each = 2)),
                       in_dir = in_dir)

args_rps <- args_nll %>% mutate(loss = "rps",
                                fname = c(rep(c(fname_cilsrps, fname_cirps), each = 10),
                                          rep(c(fname_silsrps), each = 2)))

# use different cut function for sils models
args_nll[args_nll$mod == "sils", "cutf"] <- "cutf_sils"
args_rps[args_rps$mod == "sils", "cutf"] <- "cutf_sils"


res_nll <- do.call(Map, c(f = calc_calibration, args_nll))
res_nll <- bind_rows(res_nll)

res_rps <- do.call(Map, c(f = calc_calibration, args_rps))
res_rps <- bind_rows(res_rps)

# Save results ------------------------------------------------------------

write.csv(res_nll, file = paste0(out_dir, "cal_avgnll_emp.csv"))
write.csv(res_rps, file = paste0(out_dir, "cal_avgrps_emp.csv"))
