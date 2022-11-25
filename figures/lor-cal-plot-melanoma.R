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

source("../experiments/functions/functions_DE.R")

in_dir <- "../experiments/results/DE/melanoma/"
out_dir <- "./"

splits <- 6
ensembles <- 5

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
  (rps + labs(tag = "B", subtitle = "Loss: RPS")) &
  theme(legend.position = "right", text = element_text(size = 13))
avg + plot_layout(guides = "collect")


# Load results log OR -----------------------------------------------------

indivnll <- read.csv(paste0(in_dir, "lor_nll.csv"))
indivrps <- read.csv(paste0(in_dir, "lor_rps.csv"))

indivnll$mod <- factor(indivnll$mod, levels = c("sils", "cils"))
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
      (orrps + labs(tag = "D", subtitle = "Loss: Brier score")) &
  theme(legend.position = "right", text = element_text(size = 13))
pl + plot_layout(guides = "collect", widths = c(4, 3))

ggsave(paste0(out_dir, "mela_lor_calpl_emp.pdf"), height = 8, width = 11.6)

