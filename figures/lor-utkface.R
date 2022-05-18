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

source("../experiments/functions/functions_DE.R")

in_dir <- "../experiments/results/DE/UTKFace/"
out_dir <- "./"

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

indivnll <- read.csv(paste0(in_dir, "lor_nll.csv"))
indivrps <- read.csv(paste0(in_dir, "lor_rps.csv"))

indivnll$var <- factor(indivnll$var, levels = c("gender",
                                                "x_1", "x_2", "x_3", "x_4", "x_5",
                                                "x_6", "x_7", "x_8", "x_9", "x_10"))
indivnll$mod <- factor(indivnll$mod, levels = c("sils", "silscs", "cils"))

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

