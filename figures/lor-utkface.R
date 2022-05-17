# log odds-ratios utkface
# Andrea Goetschi
# April 2022

# Dep ---------------------------------------------------------------------

library(tidyverse)
library(ggbeeswarm)
library(colorspace)
library(patchwork)
library(ggpp) # position dodgenudge

# Params ------------------------------------------------------------------

source("experiments/functions/functions_DE.R")

in_dir <- "experiments/results/DE/UTKFace/"
out_dir <- "experiments/results/DE/UTKFace/figures/"

splits <- 6
ensembles <- 5
nvars <- 11

fname_silscsnll <- "utkface_silscs_lossnll_wsyes_augno"
fname_silscsrps <- "utkface_silscs_lossrps_wsyes_augno"
fname_cilsnll <- "utkface_cils_lossnll_wsyes_augno"
fname_cilsrps <- "utkface_cils_lossrps_wsyes_augno"

fname_sils <- "utkface_sils"
fname_silsrps <- "utkface_sils_rps"

# Load results ------------------------------------------------------------

# sils
sils_files <- list.files(path = in_dir,
                         pattern = paste0(fname_sils, "_lor"))

sils_lor <- lapply(sils_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})

sils_lor <- do.call("rbind", sils_lor)
sils_lor <- sils_lor %>% gather(key = "var", value = "lor")
sils_lor$var <- factor(sils_lor$var)
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
silsrps_lor$var <- factor(silsrps_lor$var)
silsrps_lor$mod <- "sils"
silsrps_lor$spl <- factor(rep(1:splits, nvars))
silsrps_indivlor <- silsrps_lor %>% mutate(w = 1)

# silscs nll
silscsnll_files <- list.files(path = in_dir,
                              pattern = paste0(fname_silscsnll, "_lor"))

silscsnll_lor <- lapply(silscsnll_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})

# Trafo weights
w_silscsnll <- read.csv(paste0(in_dir, "w_", fname_silscsnll, ".csv"))
w_silscsnll_trf <- extract_w(w_silscsnll, meth = "trafo")
w_silscsnll_trf <- w_silscsnll_trf %>% do.call("rbind", .)

silscsnll_indivlor <- do.call("rbind", silscsnll_lor)
silscsnll_indivlor <- silscsnll_indivlor %>% gather(key = "var", value = "lor")
silscsnll_indivlor$var <- factor(silscsnll_indivlor$var)
silscsnll_indivlor$mod <- factor("silscs")
silscsnll_indivlor$spl <- factor(rep(1:splits, each = ensembles))

silscsnll_indivlor$w <- NA
for (s in seq_len(splits)) {
  for (v in unique(silscsnll_indivlor$var)) {
    silscsnll_indivlor[silscsnll_indivlor$spl == s & silscsnll_indivlor$var == v, "w"] <- w_silscsnll_trf[s, ]
  }
}

# silscs rps
silscsrps_files <- list.files(path = in_dir,
                              pattern = paste0(fname_silscsrps, "_lor"))

silscsrps_lor <- lapply(silscsrps_files, function(fname) {
  read.csv(paste0(in_dir, fname))
})

# Trafo weights
w_silscsrps <- read.csv(paste0(in_dir, "w_", fname_silscsrps, ".csv"))
w_silscsrps_trf <- extract_w(w_silscsrps, meth = "trafo")
w_silscsrps_trf <- w_silscsrps_trf %>% do.call("rbind", .)

silscsrps_indivlor <- do.call("rbind", silscsrps_lor)
silscsrps_indivlor <- silscsrps_indivlor %>% gather(key = "var", value = "lor")
silscsrps_indivlor$var <- factor(silscsrps_indivlor$var)
silscsrps_indivlor$mod <- factor("silscs")
silscsrps_indivlor$spl <- factor(rep(1:splits, each = ensembles))

silscsrps_indivlor$w <- NA
for (s in seq_len(splits)) {
  for (v in unique(silscsrps_indivlor$var)) {
    silscsrps_indivlor[silscsrps_indivlor$spl == s & silscsrps_indivlor$var == v, "w"] <- w_silscsrps_trf[s, ]
  }
}

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

# Combine -----------------------------------------------------------------

indivnll <- bindr(pat1 = "indivlor", pat2 = "nll")
indivnll$var <- factor(indivnll$var, levels = c("gender",
                                                "x_1", "x_2", "x_3", "x_4", "x_5",
                                                "x_6", "x_7", "x_8", "x_9", "x_10"))
indivnll$mod <- factor(indivnll$mod, levels = c("sils", "silscs", "cils"))


indivrps <- bindr(pat1 = "indivlor", pat2 = "rps")
indivrps$var <- factor(indivrps$var, levels = c("gender",
                                                "x_1", "x_2", "x_3", "x_4", "x_5",
                                                "x_6", "x_7", "x_8", "x_9", "x_10"))
indivrps$mod <- factor(indivrps$mod, levels = c("sils", "silscs", "cils"))

# Plot --------------------------------------------------------------------

ornll <- pl_or(indiv = indivnll, width = 1, ylim = c(-0.53, 0.5),
               order_vars = F, lbetvar = T)
orrps <- pl_or(indiv = indivrps, width = 1, ylim = c(-0.53, 0.5),
               order_vars = F, lbetvar = T)

## FIGURE E10
lor <- (ornll + labs(tag = "A", subtitle = "Loss: NLL")) +
       (orrps + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
lor + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "utkface_lor.pdf"), height = 13, width = 10)

