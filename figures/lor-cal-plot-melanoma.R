# Log odds-ratios and calibration plots melanoma
# Andrea Goetschi
# April 2022

# Dep ---------------------------------------------------------------------

library(colorspace)
library(tidyverse)
library(ggbeeswarm)
library(dplyr)
library(patchwork)
library(ggpp) # position dodgenudge

# Params ------------------------------------------------------------------

source("experiments/functions/functions_DE.R")

in_dir <- "experiments/results/DE/melanoma/"
out_dir <- "experiments/results/DE/melanoma/figures/"

splits <- 6
ensembles <- 5
nvars <- 1

fname_cilsnll <- "mela_silscs_lossnll_wsyes_augno"
fname_cilsrps <- "mela_silscs_lossrps_wsyes_augno"

fname_sils <- "mela_sils"
fname_silsrps <- "mela_sils_rps"

# Load results calibration ------------------------------------------------

avg_nll <- read.csv(paste0(in_dir, "cal_avgnll_emp.csv"))
avg_rps <- read.csv(paste0(in_dir, "cal_avgrps_emp.csv"))

# Prep --------------------------------------------------------------------

avgnll <- avg_nll %>% filter(!(mod %in% c("si", "sils"))) %>%
  mutate(mod = factor(mod, levels = c("sics", "silscs", "ci", "cils")),
         method = factor(method, levels = c("avg", "linear",
                                            "avgll", "log-linear",
                                            "avgtrf", "trafo")))

avgrps <- avg_rps %>% filter(!(mod %in% c("si", "sils"))) %>%
  mutate(mod = factor(mod, levels = c("sics", "silscs", "ci", "cils")),
         method = factor(method, levels = c("avg", "linear",
                                            "avgll", "log-linear",
                                            "avgtrf", "trafo")))

avgref <- avg_nll %>% filter(mod %in% c("sils")) %>%
  mutate(mod = factor(mod, levels = c("sils"))) %>%
  dplyr::select(-c("method"))


# Calibration plots -------------------------------------------------------

nll <- pl_cal(avg = avgnll, avg_ref = avgref)
rps <- pl_cal(avg = avgrps, avg_ref = avgref)

## FIGURE 7 A, B

avg <- (nll + labs(tag = "A", subtitle = "Loss: NLL")) /
  (rps + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
avg + plot_layout(guides = "collect")


# Load results log OR -----------------------------------------------------

# sils
sils_files <- list.files(path = in_dir,
                         pattern = paste0(fname_sils, "_lor"))

sils_lor <- lapply(sils_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})

sils_lor <- do.call("rbind", sils_lor)
sils_lor <- sils_lor %>% gather(key = "var", value = "lor")
sils_lor$var <- factor("s_age")
sils_lor$mod <- "sils"
sils_lor$spl <- factor(rep(1:splits, nvars))
silsnll_indivlor <- sils_lor %>% mutate(w = 1)

# sils rps
silsrps_files <- list.files(path = in_dir,
                            pattern = paste0(fname_silsrps, "_lor"))

silsrps_lor <- lapply(silsrps_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})

silsrps_lor <- do.call("rbind", silsrps_lor)
silsrps_lor <- silsrps_lor %>% gather(key = "var", value = "lor")
silsrps_lor$var <- factor("s_age")
silsrps_lor$mod <- "sils"
silsrps_lor$spl <- factor(rep(1:splits, nvars))
silsrps_indivlor <- silsrps_lor %>% mutate(w = 1)


# cils nll
cilsnll_files <- list.files(path = in_dir,
                            pattern = paste0(fname_cilsnll, "_lor"))

cilsnll_lor <- lapply(cilsnll_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})

# Trafo weights
w_cilsnll <- read.csv(paste0(in_dir, "w_", fname_cilsnll, ".csv"))
w_cilsnll_trf <- extract_w(w_cilsnll, meth = "trafo")
w_cilsnll_trf <- w_cilsnll_trf %>% do.call("rbind", .)

cilsnll_indivlor <- do.call("rbind", cilsnll_lor)
cilsnll_indivlor <- cilsnll_indivlor %>% gather(key = "var", value = "lor")
cilsnll_indivlor$var <- factor(cilsnll_indivlor$var)
cilsnll_indivlor$mod <- factor("cils")
cilsnll_indivlor$spl <- factor(rep(1:splits, each = ensembles))

cilsnll_indivlor$w <- NA
for (s in seq_len(splits)) {
  for (v in unique(cilsnll_indivlor$var)) {
    cilsnll_indivlor[cilsnll_indivlor$spl == s & cilsnll_indivlor$var == v, "w"] <- w_cilsnll_trf[s, ]
  }
}

# cils rps
cilsrps_files <- list.files(path = in_dir,
                            pattern = paste0(fname_cilsrps, "_lor"))

cilsrps_lor <- lapply(cilsrps_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})

# Trafo weights
w_cilsrps <- read.csv(paste0(in_dir, "w_", fname_cilsrps, ".csv"))
w_cilsrps_trf <- extract_w(w_cilsrps, meth = "trafo")
w_cilsrps_trf <- w_cilsrps_trf %>% do.call("rbind", .)

cilsrps_indivlor <- do.call("rbind", cilsrps_lor)
cilsrps_indivlor <- cilsrps_indivlor %>% gather(key = "var", value = "lor")
cilsrps_indivlor$var <- factor(cilsrps_indivlor$var)
cilsrps_indivlor$mod <- factor("cils")
cilsrps_indivlor$spl <- factor(rep(1:splits, each = ensembles))

cilsrps_indivlor$w <- NA
for (s in seq_len(splits)) {
  for (v in unique(cilsrps_indivlor$var)) {
    cilsrps_indivlor[cilsrps_indivlor$spl == s & cilsrps_indivlor$var == v, "w"] <- w_cilsrps_trf[s, ]
  }
}


indivnll <- bindr(pat1 = "indivlor", pat2 = "nll")
indivnll$mod <- factor(indivnll$mod, levels = c("sils", "cils"))

indivrps <- bindr(pat1 = "indivlor", pat2 = "rps")
indivrps$mod <- factor(indivrps$mod, levels = c("sils", "cils"))

# Plot --------------------------------------------------------------------

## FIGURE 6 C

ornll <- pl_or(indiv = indivnll, var_labs = c("s_age" = "age"),
               width = 0.2, ylim = c(0.55, 1.05))
## FIGURE 6 D

orrps <- pl_or(indiv = indivrps, var_labs = c("s_age" = "age"),
               width = 0.2, ylim = c(0.55, 1.05))

## FIGURE 6

pl <- (nll + labs(tag = "A", subtitle = "Loss: NLL")) +
      (ornll + labs(tag = "C", subtitle = "Loss: NLL")) +
      (rps + labs(tag = "B", subtitle = "Loss: Brier score")) +
      (orrps + labs(tag = "D", subtitle = "Loss: Brier score")) & theme(legend.position = "right")
pl + plot_layout(guides = "collect", widths = c(4, 3))

ggsave(paste0(out_dir, "mela_lor_calpl_emp.pdf"), height = 8, width = 11.6)

