# Calibration utkface data
# Andrea Goetschi
# April 2022

# Deps --------------------------------------------------------------------

library(etram)
library(caret)
library(Hmisc)

# Params ------------------------------------------------------------------

source("experiments/functions/functions_DE.R")
in_dir <- out_dir <- "experiments/results/DE/UTKFace/"

fname_silscsnll <- "utkface_silscs_lossnll_wsyes_augno"
fname_silscsrps <- "utkface_silscs_lossrps_wsyes_augno"
fname_cilsnll <- "utkface_cils_lossnll_wsyes_augno"
fname_cilsrps <- "utkface_cils_lossrps_wsyes_augno"
fname_sicsnll <- "utkface_sics_lossnll_wsyes_augno"
fname_sicsrps <- "utkface_sics_lossrps_wsyes_augno"
fname_cinll <- "utkface_ci_lossnll_wsyes_augno"
fname_cirps <- "utkface_ci_lossrps_wsyes_augno"

fname_sils <- "utkface_sils"
fname_silsrps <- "utkface_sils_rps"

nsemi <- 4
K <- 7

# Functions ---------------------------------------------------------------

cutf <- function(x) {
  q1 <- quantile(x, 0.05)
  q2 <- quantile(x, 0.95)
  cstep <- (q2 - q1) / 9
  cuts <- seq(from = q1, to = q2, by = cstep)
  return(cuts)
}

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

# Calibration -------------------------------------------------------------

args_nll <- data.frame(cdf_all = "all_cdf", y_true_all = "all_y",
                       mod = c(rep(c("silscs", "cils", "sics", "ci"), each = 10),
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
                       fname = c(rep(c(fname_silscsnll, fname_cilsnll,
                                       fname_sicsnll, fname_cinll), each = 10),
                                 rep(c(fname_sils), each = 2)),
                       in_dir = in_dir)

args_rps <- args_nll %>% mutate(loss = "rps",
                                fname = c(rep(c(fname_silscsrps, fname_cilsrps,
                                                fname_sicsrps, fname_cirps), each = 10),
                                          rep(c(fname_silsrps), each = 2)))


res_nll <- do.call(Map, c(f = calc_calibration, args_nll))
res_nll <- bind_rows(res_nll)

res_rps <- do.call(Map, c(f = calc_calibration, args_rps))
res_rps <- bind_rows(res_rps)

# Save results ------------------------------------------------------------

write.csv(res_nll, file = paste0(out_dir, "cal_avgnll.csv"))
write.csv(res_rps, file = paste0(out_dir, "cal_avgrps.csv"))
