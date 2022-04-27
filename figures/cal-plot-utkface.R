# Calibration plots utkface
# Andrea Goetschi
# April 2022

# Dep ---------------------------------------------------------------------

library(colorspace)
library(tidyverse)
library(patchwork)

# Paths -------------------------------------------------------------------

source("experiments/functions/functions_DE.R")

in_dir <- "experiments/results/DE/UTKFace/"
out_dir <- "experiments/results/DE/UTKFace/figures/"

# Load results ------------------------------------------------------------

avg_nll <- read.csv(paste0(in_dir, "cal_avgnll.csv"))
avg_rps <- read.csv(paste0(in_dir, "cal_avgrps.csv"))

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
  select(-c("method"))

# Plots -------------------------------------------------------------------

# mean cal across all 6 splits
nll <- pl_cal(avg = avgnll, avg_ref = avgref)

# mean cal across all 6 splits
rps <- pl_cal(avg = avgrps, avg_ref = avgref)

## FIGURE 10

avg <- (nll + labs(tag = "A", subtitle = "Loss: NLL")) / 
       (rps + labs(tag = "B", subtitle = "Loss: RPS")) & theme(legend.position = "right")
avg + plot_layout(guides = "collect")
ggsave(paste0(out_dir, "utkface_calpl.pdf"), height = 10.5, width = 12.5)

