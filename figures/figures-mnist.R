# All Figures MNIST
# Andrea Goetschi
# April 2022

# Dependencies ------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(colorspace)
library(ggbeeswarm)
library(etram)

# Params -------------------------------------------------------------------

source("experiments/functions/functions_DE.R")

in_dir <- "experiments/results/DE/MNIST/"
out_dir <- "experiments/results/DE/MNIST/figures/"

fname_cinll <- "mnist_ci_lossnll_wsno_augno"
fname_cirps <- "mnist_ci_lossrps_wsno_augno"

spl <- 6
ens <- 5

# Performance equal weights -----------------------------------------------

met_cinll <- read.csv(file = paste0(in_dir, "met_", fname_cinll, ".csv"))
met_cinll <- met_cinll[met_cinll$topn == ens, -2]
met_cinll$mod <- "ci"
met_cinll$weights <- "equal"

indivmet_cinll <- read.csv(file = paste0(in_dir, "indivmet_", fname_cinll, ".csv"))
indivmet_cinll$mod <- "ci"

met_cirps <- read.csv(file = paste0(in_dir, "met_", fname_cirps, ".csv"))
met_cirps <- met_cirps[met_cirps$topn == ens, -2]
met_cirps$mod <- "ci"
met_cirps$weights <- "equal"

indivmet_cirps <- read.csv(file = paste0(in_dir, "indivmet_", fname_cirps, ".csv"))
indivmet_cirps$mod <- "ci"

## Bootstrap confidence intervals

# ci nll
f_nll_nw <- list.files(in_dir, pattern = "^boot.*nll.*\\.csv$") # ^: start, $: end, .*: any pattern,
lys_nll_nw <- lapply(f_nll_nw, function(x) read.csv(paste0(in_dir, x)))
ci_nll_nw <- do.call("rbind", lys_nll_nw)
ci_nll_nw$weights <- "equal"

# ci rps
f_rps_nw <- list.files(in_dir, pattern = "^boot.*rps.*\\.csv$")
lys_rps_nw <- lapply(f_rps_nw, function(x) read.csv(paste0(in_dir, x)))
ci_rps_nw <- do.call("rbind", lys_rps_nw)
ci_rps_nw$weights <- "equal"


# Performance weighted ----------------------------------------------------

met_cinll_w <- read.csv(file = paste0(in_dir, "met_", fname_cinll, "_weighted.csv"))
met_cinll_w$mod <- "ci"
met_cinll_w[met_cinll_w$method == "avgl", "method"] <- "avg"
met_cinll_w$weights <- "tuned"

met_cirps_w <- read.csv(file = paste0(in_dir, "met_", fname_cirps, "_weighted.csv"))
met_cirps_w$mod <- "ci"
met_cirps_w[met_cirps_w$method == "avgl", "method"] <- "avg"
met_cirps_w$weights <- "tuned"

## Bootstrap confidence intervals

# ci nll
f_nll_w <- list.files(in_dir, pattern = "^wboot.*nll.*\\.csv$")
lys_nll_w <- lapply(f_nll_w, function(x) read.csv(paste0(in_dir, x)))
ci_nll_w <- do.call("rbind", lys_nll_w)
ci_nll_w$weights <- "tuned"

# ci rps
f_rps_w <- list.files(in_dir, pattern = "^wboot.*rps.*\\.csv$")
lys_rps_w <- lapply(f_rps_w, function(x) read.csv(paste0(in_dir, x)))
ci_rps_w <- do.call("rbind", lys_rps_w)
ci_rps_w$weights <- "tuned"

# Combine -----------------------------------------------------------------

met_cinll_all <- bind_rows(met_cinll, met_cinll_w)
met_cinll_all$spl <- factor(met_cinll_all$spl)
met_cirps_all <- bind_rows(met_cirps, met_cirps_w)
met_cirps_all$spl <- factor(met_cirps_all$spl)

ci_nll <- bind_rows(ci_nll_nw, ci_nll_w)
ci_rps <- bind_rows(ci_rps_nw, ci_rps_w)

## Reorder levels

ord_met <- c("nll", "rps", "eacc")
cal_met <- c("cint", "cslope")
meths <- c("trafo", "avgtrf",
           "log-linear", "avgll",
           "linear", "avg")

met_cinll_all <- relev(met_cinll_all, "metric", ord_met)
indivmet_cinll <- relev(indivmet_cinll, "metric", ord_met)
met_cinll_all <- relev(met_cinll_all, "method", meths)
ci_nll <- relev(ci_nll, "method", meths)

met_cirps_all <- relev(met_cirps_all, "metric", ord_met)
indivmet_cirps <- relev(indivmet_cirps, "metric", ord_met)
met_cirps_all <- relev(met_cirps_all, "method", meths)
ci_rps <- relev(ci_rps, "method", meths)


# Performance plots -------------------------------------------------------

# Plots NLL

pl_ordnll <- pl_met(spl_met = met_cinll_all,
                    ci = ci_nll,
                    xlab = "",
                    metrics = ord_met)

pl_ordnll_indiv <- pl_met(spl_met = met_cinll_all, indiv_met = indivmet_cinll,
                          ci = ci_nll,
                          xlab = "",
                          metrics = ord_met)

# Plots RPS

pl_ordrps <- pl_met(spl_met = met_cirps_all,
                    ci = ci_rps,
                    xlab = "",
                    metrics = ord_met)

pl_ordrps_indiv <- pl_met(spl_met = met_cirps_all,
                          ci = ci_rps,
                          indiv_met = indivmet_cirps,
                          xlab = "",
                          metrics = ord_met)


## FIGURE E1 A, B

c_ord <- (pl_ordnll + labs(tag = "A", subtitle = "Loss: NLL")) /
         (pl_ordrps + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_ord + plot_layout(guides = "collect")
# ggsave(paste0(out_dir, "mnist_ci_wvsnw.pdf"), height = 9, width = 8)

## FIGURE E2

c_ord_indiv <- (pl_ordnll_indiv + labs(tag = "A", subtitle = "Loss: NLL")) /
               (pl_ordrps_indiv + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_ord_indiv + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mnist_ci_wvsnw_indiv.pdf"), height = 7, width = 8)


# Calibration intercept, slope --------------------------------------------

# Plots NLL

pl_calnll_indiv <- pl_met(spl_met = met_cinll_all,
                          metrics = cal_met,
                          ci = ci_nll,
                          indiv_met = indivmet_cinll,
                          ylab = "", xlab = "")

# Plots RPS

pl_calrps_indiv <- pl_met(spl_met = met_cirps_all,
                          metrics = cal_met,
                          ci = ci_rps,
                          indiv_met = indivmet_cirps,
                          ylab = "", xlab = "")


## FIGURE E3

c_cal_indiv <- (pl_calnll_indiv + labs(tag = "A", subtitle = "Loss: NLL")) /
               (pl_calrps_indiv + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_cal_indiv + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mnist_cal_indiv.pdf"), height = 7, width = 7.5)


# Calibration plots -------------------------------------------------------

## Read results

avg_nll <- read.csv(paste0(in_dir, "cumcal_avgnll.csv"))
avg_rps <- read.csv(paste0(in_dir, "cumcal_avgrps.csv"))

## Prep

avgnll <- avg_nll %>% mutate(method = factor(method, levels = c("avg", "linear",
                                                                "avgll", "log-linear",
                                                                "avgtrf", "trafo")))
avgrps <- avg_rps %>% mutate(method = factor(method, levels = c("avg", "linear",
                                                                "avgll", "log-linear",
                                                                "avgtrf", "trafo")))
nll <- pl_cal(avg = avgnll)

rps <- pl_cal(avg = avgrps)


## FIGURE E1 C, D

avg <- (nll + labs(tag = "A", subtitle = "Loss: NLL")) +
       (rps + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
avg + plot_layout(guides = "collect")
# ggsave(paste0(out_dir, "mnist_cumcalpl.pdf"), height = 6.5, width = 8.7)


# Combine performance and calibration plots -------------------------------

## FIGURE E1

prf_cal <- (pl_ordnll + labs(tag = "A", subtitle = "Loss: NLL")) +
           (nll + theme(legend.position = "none")  + labs(tag = "C", subtitle = "Loss: NLL")) +
           (pl_ordrps + labs(tag = "B", subtitle = "Loss: RPS")) +
           (rps + theme(legend.position = "none") + labs(tag = "D", subtitle = "Loss: RPS"))
prf_cal + plot_layout(guides = "collect", widths = c(4, 1))

ggsave(paste0(out_dir, "mnist_prf_cumcalpl.pdf"), height = 8, width = 12)
