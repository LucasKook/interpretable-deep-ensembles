# Test error and calibration intercept, slope utkface
# Andrea Goetschi
# April 2022

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(colorspace)
library(ggbeeswarm)
library(patchwork)

# Params -------------------------------------------------------------------

source("experiments/functions/functions_DE.R")

in_dir <- "experiments/results/DE/UTKFace/"
out_dir <- "experiments/results/DE/UTKFace/figures/"

fname_silscsnll <- "utkface_silscs_lossnll_wsyes_augno"
fname_silscsrps <- "utkface_silscs_lossrps_wsyes_augno"
fname_cilsnll <- "utkface_cils_lossnll_wsyes_augno"
fname_cilsrps <- "utkface_cils_lossrps_wsyes_augno"
fname_sicsnll <- "utkface_sics_lossnll_wsyes_augno"
fname_sicsrps <- "utkface_sics_lossrps_wsyes_augno"
fname_cinll <- "utkface_ci_lossnll_wsyes_augno"
fname_cirps <- "utkface_ci_lossrps_wsyes_augno"

fname_sils <- "utkface_sils"
fname_si <- "utkface_si"
fname_silsrps <- "utkface_sils_rps"
fname_sirps <- "utkface_si_rps"

spl <- 6
ens <- 5

# Performance results equal weights ---------------------------------------

## NLL

### SILSCS

met_silscsnll <- read.csv(file = paste0(in_dir, "met_", fname_silscsnll, ".csv"))
met_silscsnll <- met_silscsnll[met_silscsnll$topn == ens, -2]
met_silscsnll$mod <- "silscs"
met_silscsnll$weights <- "equal"

indiv_silscsnll <- read.csv(file = paste0(in_dir, "indivmet_", fname_silscsnll, ".csv"))
indiv_silscsnll$mod <- "silscs"

### CILS

met_cilsnll <- read.csv(file = paste0(in_dir, "met_", fname_cilsnll, ".csv"))
met_cilsnll <- met_cilsnll[met_cilsnll$topn == ens, -2]
met_cilsnll$mod <- "cils"
met_cilsnll$weights <- "equal"

indiv_cilsnll <- read.csv(file = paste0(in_dir, "indivmet_", fname_cilsnll, ".csv"))
indiv_cilsnll$mod <- "cils"

### SICS

met_sicsnll <- read.csv(file = paste0(in_dir, "met_", fname_sicsnll, ".csv"))
met_sicsnll <- met_sicsnll[met_sicsnll$topn == ens, -2]
met_sicsnll$mod <- "sics"
met_sicsnll$weights <- "equal"

indiv_sicsnll <- read.csv(file = paste0(in_dir, "indivmet_", fname_sicsnll, ".csv"))
indiv_sicsnll$mod <- "sics"

### CI

met_cinll <- read.csv(file = paste0(in_dir, "met_", fname_cinll, ".csv"))
met_cinll <- met_cinll[met_cinll$topn == ens, -2]
met_cinll$mod <- "ci"
met_cinll$weights <- "equal"

indiv_cinll <- read.csv(file = paste0(in_dir, "indivmet_", fname_cinll, ".csv"))
indiv_cinll$mod <- "ci"

### SI

met_sinll <-  read.csv(file = paste0(in_dir, "met_", fname_si, ".csv"))
met_sinll$mod <- "si"
met_sinll$weights <- "equal"

### SILS

met_silsnll <-  read.csv(file = paste0(in_dir, "met_", fname_sils, ".csv"))
met_silsnll$mod <- "sils"
met_silsnll$weights <- "equal"


## RPS

### SILSCS

met_silscsrps <- read.csv(file = paste0(in_dir, "met_", fname_silscsrps, ".csv"))
met_silscsrps <- met_silscsrps[met_silscsrps$topn == ens, -2]
met_silscsrps$mod <- "silscs"
met_silscsrps$weights <- "equal"

indiv_silscsrps <- read.csv(file = paste0(in_dir, "indivmet_", fname_silscsrps, ".csv"))
indiv_silscsrps$mod <- "silscs"

### CILS

met_cilsrps <- read.csv(file = paste0(in_dir, "met_", fname_cilsrps, ".csv"))
met_cilsrps <- met_cilsrps[met_cilsrps$topn == ens, -2]
met_cilsrps$mod <- "cils"
met_cilsrps$weights <- "equal"

indiv_cilsrps <- read.csv(file = paste0(in_dir, "indivmet_", fname_cilsrps, ".csv"))
indiv_cilsrps$mod <- "cils"

### SICS

met_sicsrps <- read.csv(file = paste0(in_dir, "met_", fname_sicsrps, ".csv"))
met_sicsrps <- met_sicsrps[met_sicsrps$topn == ens, -2]
met_sicsrps$mod <- "sics"
met_sicsrps$weights <- "equal"

indiv_sicsrps <- read.csv(file = paste0(in_dir, "indivmet_", fname_sicsrps, ".csv"))
indiv_sicsrps$mod <- "sics"

### CI

met_cirps <- read.csv(file = paste0(in_dir, "met_", fname_cirps, ".csv"))
met_cirps <- met_cirps[met_cirps$topn == ens, -2]
met_cirps$mod <- "ci"
met_cirps$weights <- "equal"

indiv_cirps <- read.csv(file = paste0(in_dir, "indivmet_", fname_cirps, ".csv"))
indiv_cirps$mod <- "ci"

### SI

met_sirps <-  read.csv(file = paste0(in_dir, "met_", fname_sirps, ".csv"))
met_sirps$mod <- "si"
met_sirps$weights <- "equal"

### SILS

met_silsrps <-  read.csv(file = paste0(in_dir, "met_", fname_silsrps, ".csv"))
met_silsrps$mod <- "sils"
met_silsrps$weights <- "equal"


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


# Performance results weighted --------------------------------------------

## NLL

### SILSCS

met_silscsnll_w <- read.csv(file = paste0(in_dir, "met_", fname_silscsnll, "_weighted.csv"))
met_silscsnll_w$mod <- "silscs"
met_silscsnll_w$weights <- "tuned"

### CILS

met_cilsnll_w <- read.csv(file = paste0(in_dir, "met_", fname_cilsnll, "_weighted.csv"))
met_cilsnll_w$mod <- "cils"
met_cilsnll_w$weights <- "tuned"

### SICS

met_sicsnll_w <- read.csv(file = paste0(in_dir, "met_", fname_sicsnll, "_weighted.csv"))
met_sicsnll_w$mod <- "sics"
met_sicsnll_w$weights <- "tuned"

### CI

met_cinll_w <- read.csv(file = paste0(in_dir, "met_", fname_cinll, "_weighted.csv"))
met_cinll_w$mod <- "ci"
met_cinll_w$weights <- "tuned"

### SI

met_sinll_w <- met_sinll
met_sinll_w$weights <- "tuned"

### SILS

met_silsnll_w <- met_silsnll
met_silsnll_w$weights <- "tuned"


## RPS

### SILSCS

met_silscsrps_w <- read.csv(file = paste0(in_dir, "met_", fname_silscsrps, "_weighted.csv"))
met_silscsrps_w$mod <- "silscs"
met_silscsrps_w$weights <- "tuned"

### CILS

met_cilsrps_w <- read.csv(file = paste0(in_dir, "met_", fname_cilsrps, "_weighted.csv"))
met_cilsrps_w$mod <- "cils"
met_cilsrps_w$weights <- "tuned"

### SICS

met_sicsrps_w <- read.csv(file = paste0(in_dir, "met_", fname_sicsrps, "_weighted.csv"))
met_sicsrps_w$mod <- "sics"
met_sicsrps_w$weights <- "tuned"

### CI

met_cirps_w <- read.csv(file = paste0(in_dir, "met_", fname_cirps, "_weighted.csv"))
met_cirps_w$mod <- "ci"
met_cirps_w$weights <- "tuned"

### SI

met_sirps_w <- met_sirps
met_sirps_w$weights <- "tuned"

### SILS

met_silsrps_w <- met_silsrps
met_silsrps_w$weights <- "tuned"


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

met_negloglik <- bindr(pat1 = "met", pat2 = "nll")
met_negloglik[met_negloglik$method == "avgl", "method"] <- "avg"
met_negloglik[met_negloglik$mod %in% c("si", "sils"), "method"] <- "ref"
met_negloglik$spl <- factor(met_negloglik$spl)

ind_negloglik <- bindr(pat1 = "indiv", pat2 = "nll")

met_ranked <- bindr(pat1 = "met", pat2 = "rps")
met_ranked[met_ranked$method == "avgl", "method"] <- "avg"
met_ranked[met_ranked$mod %in% c("si", "sils"), "method"] <- "ref"
met_ranked$spl <- factor(met_ranked$spl)

ind_ranked <- bindr(pat1 = "indiv", pat2 = "rps")

# confidence intervals
ci_nll <- bind_rows(ci_nll_nw, ci_nll_w, ci_nll_nw %>% 
                      filter(mod %in% c("si", "sils")) %>%
                      mutate(weights = "tuned"))
ci_rps <- bind_rows(ci_rps_nw, ci_rps_w, ci_rps_nw %>% 
                      filter(mod %in% c("si", "sils")) %>%
                      mutate(weights = "tuned"))

# Reorder levels ----------------------------------------------------------

ord_metrics <- c("nll", "rps", "eqwk")
cal_metrics <- c("cint", "cslope")
meths <- c("trafo", "avgtrf", 
           "log-linear", "avgll",
           "linear", "avg")
mods <- c("si", "sils", "sics", "silscs", "ci", "cils")

met_negloglik <- relev(met_negloglik, "metric", c(ord_metrics))
met_negloglik <- relev(met_negloglik, "method", meths)
ci_nll <- relev(ci_nll, "method", meths)
met_negloglik <- relev(met_negloglik, "mod", mods)
ind_negloglik <- relev(ind_negloglik, "metric", c(ord_metrics))

met_ranked <- relev(met_ranked, "metric", c(ord_metrics))
met_ranked <- relev(met_ranked, "method", meths)
ci_rps <- relev(ci_rps, "method", meths)
met_ranked <- relev(met_ranked, "mod", mods)
ind_ranked <- relev(ind_ranked, "metric", c(ord_metrics))


# Absolute plots ----------------------------------------------------------

# without reference models
met_negloglik_noref <- met_negloglik[!(met_negloglik$mod %in% c("si", "sils")), ]
met_negloglik_noref$mod <- factor(met_negloglik_noref$mod, levels = c("sics", "silscs", "ci", "cils"))
met_negloglik_noref$method <- factor(met_negloglik_noref$method, levels = c("trafo", "avgtrf", 
                                                                            "log-linear", "avgll",
                                                                            "linear", "avg"))
pl_ordnll_noref <- pl_met(spl_met = met_negloglik_noref,
                          ci = ci_nll,
                          metrics = ord_metrics)

# with indiv
pl_ordnll_indiv <- pl_met(spl_met = met_negloglik,
                          indiv_met = ind_negloglik,
                          ci = ci_nll,
                          metrics = ord_metrics,
                          ref = c("si", "sils"))

# without reference models
met_ranked_noref <- met_ranked[!(met_ranked$mod %in% c("si", "sils")), ]
met_ranked_noref$mod <- factor(met_ranked_noref$mod, levels = c("sics", "silscs", "ci", "cils"))
met_ranked_noref$method <- factor(met_ranked_noref$method, levels = c("trafo", "avgtrf", 
                                                                      "log-linear", "avgll",
                                                                      "linear", "avg"))
pl_ordrps_noref <- pl_met(spl_met = met_ranked_noref[is.finite(met_ranked_noref$val), ], 
                          ci = ci_rps,
                          metrics = ord_metrics)

# with indiv
pl_ordrps_indiv <- pl_met(spl_met = met_ranked[is.finite(met_ranked$val), ],
                          indiv_met = ind_ranked[is.finite(ind_ranked$val), ],
                          ci = ci_rps,
                          metrics = ord_metrics,
                          ref = c("si", "sils"))

## FIGURE 9

c_ord <- (pl_ordnll_noref + labs(tag = "A", subtitle = "Loss: NLL")) / 
         (pl_ordrps_noref + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_ord + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "utkface_wvsnw.pdf"), height = 13.5, width = 11.5)

## FIGURE E7

c_ord_indiv <- (pl_ordnll_indiv + labs(tag = "A", subtitle = "Loss: NLL")) / 
               (pl_ordrps_indiv + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_ord_indiv + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "utkface_wvsnw_indiv.pdf"), height = 13.5, width = 11.5)


# Relative plots ----------------------------------------------------------

pl_ordnll_rel <- pl_met(spl_met = met_negloglik,
                        metrics = ord_metrics,
                        ci = ci_nll,
                        ref = "si",
                        rel = TRUE)

pl_ordrps_rel <-  pl_met(spl_met = met_ranked[is.finite(met_ranked$val), ],
                         metrics = ord_metrics,
                         ci = ci_rps,
                         ref = "si",
                         rel = TRUE)

## FIGURE E6

c_ord_rel <- (pl_ordnll_rel + labs(tag = "A", subtitle = "Loss: NLL")) / 
             (pl_ordrps_rel + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_ord_rel + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "utkface_wvsnw_rel.pdf"), height = 13, width = 11.5)


# Calibration int, slope --------------------------------------------------

pl_calnll_indiv <- pl_met(spl_met = met_negloglik,
                          metrics = cal_metrics,
                          ci = ci_nll,
                          indiv_met = ind_negloglik,
                          ref = c("si", "sils"),
                          ylab = "")

pl_calrps_indiv <- pl_met(spl_met = met_ranked,
                          metrics = cal_metrics,
                          ci = ci_rps,
                          indiv_met = ind_ranked,
                          ref = c("si", "sils"),
                          ylab = "")

## FIGURE E8

c_cal_indiv <- (pl_calnll_indiv + labs(tag = "A", subtitle = "Loss: NLL")) / 
               (pl_calrps_indiv + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
c_cal_indiv + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "utkface_cal_indiv.pdf"), height = 12.1, width = 7.5)
