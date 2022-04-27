# Test error and calibration intercept, slope melanoma
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

in_dir <- "experiments/results/DE/melanoma/"
out_dir <- "experiments/results/DE/melanoma/figures/"

fname_cilsnll <- "mela_silscs_lossnll_wsyes_augno"
fname_cilsrps <- "mela_silscs_lossrps_wsyes_augno"
fname_cinll <- "mela_ci_lossnll_wsyes_augno"
fname_cirps <- "mela_ci_lossrps_wsyes_augno"
fname_si <- "mela_si"
fname_sils <- "mela_sils"

spl <- 6
ens <- 5

# Results equal weights ---------------------------------------------------

## NLL

### CILS

met_cilsnll <- read.csv(file = paste0(in_dir, "met_", fname_cilsnll, ".csv"))
met_cilsnll <- met_cilsnll[met_cilsnll$topn == ens, -2]
met_cilsnll$mod <- "cils"
met_cilsnll$weights <- "equal"

indiv_cilsnll <- read.csv(file = paste0(in_dir, "indivmet_", fname_cilsnll, ".csv"))
indiv_cilsnll$mod <- "cils"

### CI

met_cinll <- read.csv(file = paste0(in_dir, "met_", fname_cinll, ".csv"))
met_cinll <- met_cinll[met_cinll$topn == ens, -2]
met_cinll$mod <- "ci"
met_cinll$weights <- "equal"

indiv_cinll <- read.csv(file = paste0(in_dir, "indivmet_", fname_cinll, ".csv"))
indiv_cinll$mod <- "ci"

### SI

met_si <-  read.csv(file = paste0(in_dir, "met_", fname_si, ".csv"))
met_si$mod <- "si"
met_si$weights <- "equal"

### SILS

met_sils <-  read.csv(file = paste0(in_dir, "met_", fname_sils, ".csv"))
met_sils$mod <- "sils"
met_sils$weights <- "equal"


## RPS

### CILS

met_cilsrps <- read.csv(file = paste0(in_dir, "met_", fname_cilsrps, ".csv"))
met_cilsrps <- met_cilsrps[met_cilsrps$topn == ens, -2]
met_cilsrps$mod <- "cils"
met_cilsrps$weights <- "equal"

indiv_cilsrps <- read.csv(file = paste0(in_dir, "indivmet_", fname_cilsrps, ".csv"))
indiv_cilsrps$mod <- "cils"

### CI

met_cirps <- read.csv(file = paste0(in_dir, "met_", fname_cirps, ".csv"))
met_cirps <- met_cirps[met_cirps$topn == ens, -2]
met_cirps$mod <- "ci"
met_cirps$weights <- "equal"

indiv_cirps <- read.csv(file = paste0(in_dir, "indivmet_", fname_cirps, ".csv"))
indiv_cirps$mod <- "ci"

# Results weighted --------------------------------------------------------

## NLL

### CILS

met_cilsnll_w <- read.csv(file = paste0(in_dir, "met_", fname_cilsnll, "_weighted.csv"))
met_cilsnll_w$mod <- "cils"
met_cilsnll_w$weights <- "tuned"

### CI

met_cinll_w <- read.csv(file = paste0(in_dir, "met_", fname_cinll, "_weighted.csv"))
met_cinll_w$mod <- "ci"
met_cinll_w$weights <- "tuned"

### SI

met_si_w <- met_si
met_si_w$weights <- "tuned"

### SILS

met_sils_w <- met_sils
met_sils_w$weights <- "tuned"


## RPS

### CILS

met_cilsrps_w <- read.csv(file = paste0(in_dir, "met_", fname_cilsrps, "_weighted.csv"))
met_cilsrps_w$mod <- "cils"
met_cilsrps_w$weights <- "tuned"

### CI

met_cirps_w <- read.csv(file = paste0(in_dir, "met_", fname_cirps, "_weighted.csv"))
met_cirps_w$mod <- "ci"
met_cirps_w$weights <- "tuned"

# Combine -----------------------------------------------------------------

met_negloglik <- bindr(pat1 = "met", pat2 = "nll")
met_negloglik[met_negloglik$method == "avgl", "method"] <- "avg"
met_negloglik <- bind_rows(met_negloglik, met_si, met_si_w, met_sils, met_sils_w)
met_negloglik[met_negloglik$mod %in% c("si", "sils"), "method"] <- "ref"
met_negloglik$spl <- factor(met_negloglik$spl)

ind_negloglik <- bindr(pat1 = "indiv", pat2 = "nll")

met_ranked <- bindr(pat1 = "met", pat2 = "rps")
met_ranked[met_ranked$method == "avgl", "method"] <- "avg"
met_ranked <- bind_rows(met_ranked,met_si, met_si_w, met_sils, met_sils_w)
met_ranked[met_ranked$mod %in% c("si", "sils"), "method"] <- "ref"
met_ranked$spl <- factor(met_ranked$spl)

ind_ranked <- bindr(pat1 = "indiv", pat2 = "rps")

# Reorder levels ----------------------------------------------------------

prf_metrics <- c("nll", "brier", "eauc", "ebinacc")
cal_metrics <- c("cint", "cslope")
meths <- c("trafo", "avgtrf", 
           "log-linear", "avgll",
           "linear", "avg")
mods <- c("si", "sils", "ci", "cils")

met_negloglik <- relev(met_negloglik, "metric", c(prf_metrics))
met_negloglik <- relev(met_negloglik, "method", meths)
met_negloglik <- relev(met_negloglik, "mod", mods)
ind_negloglik <- relev(ind_negloglik, "metric", prf_metrics)

met_ranked <- relev(met_ranked, "method", meths)
met_ranked <- relev(met_ranked, "metric", c(prf_metrics))
met_ranked <- relev(met_ranked, "mod", mods)
ind_ranked <- relev(ind_ranked, "metric", prf_metrics)


# Absolute performance plots ----------------------------------------------

pl_nll <- pl_met(spl_met = met_negloglik, 
                    metrics = prf_metrics,
                    ref = c("si", "sils"))

# with indiv
pl_nll_indiv <- pl_met(spl_met = met_negloglik,
                          indiv_met = ind_negloglik,
                          metrics = prf_metrics,
                          ref = c("si", "sils"))

pl_rps <-  pl_met(spl_met = met_ranked,
                     metrics = prf_metrics,
                     ref = c("si", "sils"))

# with indiv
pl_rps_indiv <- pl_met(spl_met = met_ranked,
                          indiv_met = ind_ranked,
                          metrics = prf_metrics,
                          ref = c("si", "sils"))


## FIGURE 6

c_prf <- (pl_nll + labs(tag = "A", subtitle = "Loss: NLL")) / 
         (pl_rps + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_prf + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mela_wvsnw.pdf"), height = 13, width = 13.5)

## FIGURE E4

c_prf_indiv <- (pl_nll_indiv + labs(tag = "A", subtitle = "Loss: NLL")) / 
               (pl_rps_indiv + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_prf_indiv + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mela_wvsnw_indiv.pdf"), height = 13, width = 13.5)


# Relative performance plots ----------------------------------------------

pl_nll_rel <- pl_met(spl_met = met_negloglik,
                        metrics = prf_metrics,
                        ref = "si",
                        rel = TRUE)

pl_rps_rel <-  pl_met(spl_met = met_ranked,
                         metrics = prf_metrics,
                         ref = "si",
                         rel = TRUE)

## FIGURE E3

c_prf_rel <- (pl_nll_rel + labs(tag = "A", subtitle = "Loss: NLL")) / 
             (pl_rps_rel + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_prf_rel + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mela_wvsnw_rel.pdf"), height = 12.5, width = 13.5)


# Calibration intercept, slope --------------------------------------------

pl_calnll_indiv <- pl_met(spl_met = met_negloglik,
                          metrics = cal_metrics,
                          indiv_met = ind_negloglik,
                          ref = c("si", "sils"),
                          ylab = "")

pl_calrps_indiv <- pl_met(spl_met = met_ranked,
                          metrics = cal_metrics,
                          indiv_met = ind_ranked,
                          ref = c("si", "sils"),
                          ylab = "")


## FIGURE E5
c_cal_indiv <- (pl_calnll_indiv + labs(tag = "A", subtitle = "Loss: NLL")) / 
               (pl_calrps_indiv + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_cal_indiv + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mela_cal_indiv.pdf"), height = 12, width = 7.5)

