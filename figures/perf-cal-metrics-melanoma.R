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

source("../experiments/functions/functions_DE.R")

in_dir <- "../experiments/results/DE/melanoma/"
out_dir <- "./"

fname_cilsnll <- "mela_silscs_lossnll_wsyes_augno"
fname_cilsrps <- "mela_silscs_lossrps_wsyes_augno"
fname_cinll <- "mela_ci_lossnll_wsyes_augno"
fname_cirps <- "mela_ci_lossrps_wsyes_augno"
fname_sinll <- "mela_si"
fname_silsnll <- "mela_sils"
fname_sirps <- "mela_si_rps"
fname_silsrps <- "mela_sils_rps"

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

met_sinll <-  read.csv(file = paste0(in_dir, "met_", fname_sinll, ".csv"))
met_sinll$mod <- "si"
met_sinll$weights <- "equal"

### SILS

met_silsnll <-  read.csv(file = paste0(in_dir, "met_", fname_silsnll, ".csv"))
met_silsnll$mod <- "sils"
met_silsnll$weights <- "equal"


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

met_sinll_w <- met_sinll
met_sinll_w$weights <- "tuned"

### SILS

met_silsnll_w <- met_silsnll
met_silsnll_w$weights <- "tuned"


## RPS

### CILS

met_cilsrps_w <- read.csv(file = paste0(in_dir, "met_", fname_cilsrps, "_weighted.csv"))
met_cilsrps_w$mod <- "cils"
met_cilsrps_w$weights <- "tuned"

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

prf_metrics <- c("nll", "brier", "eauc", "eacc")
cal_metrics <- c("cint", "cslope")
meths <- c("trafo", "avgtrf",
           "log-linear", "avgll",
           "linear", "avg")
mods <- c("si", "sils", "ci", "cils")

met_negloglik <- relev(met_negloglik, "metric", c(prf_metrics))
met_negloglik <- relev(met_negloglik, "method", meths)
ci_nll <- relev(ci_nll, "method", meths)
met_negloglik <- relev(met_negloglik, "mod", mods)
ind_negloglik <- relev(ind_negloglik, "metric", prf_metrics)

met_ranked <- relev(met_ranked, "metric", c(prf_metrics))
met_ranked <- relev(met_ranked, "method", meths)
ci_rps <- relev(ci_rps, "method", meths)
met_ranked <- relev(met_ranked, "mod", mods)
ind_ranked <- relev(ind_ranked, "metric", prf_metrics)


# Absolute performance plots ----------------------------------------------

pl_nll <- pl_met(spl_met = met_negloglik,
                 metrics = prf_metrics,
                 ci = ci_nll,
                 ref = c("si", "sils"))

# with indiv
pl_nll_indiv <- pl_met(spl_met = met_negloglik,
                       indiv_met = ind_negloglik,
                       ci = ci_nll,
                       metrics = prf_metrics,
                       ref = c("si", "sils"))

pl_rps <-  pl_met(spl_met = met_ranked,
                  metrics = prf_metrics,
                  ci = ci_rps,
                  ref = c("si", "sils"))

# with indiv
pl_rps_indiv <- pl_met(spl_met = met_ranked,
                       indiv_met = ind_ranked,
                       ci = ci_rps,
                       metrics = prf_metrics,
                       ref = c("si", "sils"))


## FIGURE 5

c_prf <- (pl_nll + labs(tag = "A", subtitle = "Loss: NLL")) /
         (pl_rps + labs(tag = "B", subtitle = "Loss: Brier score")) &
  theme(legend.position = "right", text = element_text(size = 13))
c_prf + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mela_wvsnw.pdf"), height = 13, width = 13.5)

## FIGURE E5

c_prf_indiv <- (pl_nll_indiv + labs(tag = "A", subtitle = "Loss: NLL")) /
               (pl_rps_indiv + labs(tag = "B", subtitle = "Loss: Brier score")) &
  theme(legend.position = "right", text = element_text(size = 13))
c_prf_indiv + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mela_wvsnw_indiv.pdf"), height = 13, width = 13.5)


# Relative performance plots ----------------------------------------------

pl_nll_rel <- pl_met(spl_met = met_negloglik,
                     metrics = prf_metrics,
                     ci = ci_nll,
                     ref = "si",
                     rel = TRUE)

pl_rps_rel <-  pl_met(spl_met = met_ranked,
                      metrics = prf_metrics,
                      ci = ci_rps,
                      ref = "si",
                      rel = TRUE)

## FIGURE E4

c_prf_rel <- (pl_nll_rel + labs(tag = "A", subtitle = "Loss: NLL")) /
             (pl_rps_rel + labs(tag = "B", subtitle = "Loss: Brier score")) &
  theme(legend.position = "right", text = element_text(size = 13))
c_prf_rel + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mela_wvsnw_rel.pdf"), height = 12.5, width = 13.5)


# Calibration intercept, slope --------------------------------------------

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


## FIGURE E6
c_cal_indiv <- (pl_calnll_indiv + labs(tag = "A", subtitle = "Loss: NLL")) /
               (pl_calrps_indiv + labs(tag = "B", subtitle = "Loss: Brier score")) &
  theme(legend.position = "right", text = element_text(size = 13))
c_cal_indiv + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "mela_cal_indiv.pdf"), height = 12, width = 7.5)

