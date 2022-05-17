# Fit unconditional (SI) and linear shift model (SI-LSx) by minimizing RPS
# Andrea Goetschi
# May 2022

# Reproducibility ---------------------------------------------------------

set.seed(738593)

# Deps --------------------------------------------------------------------

library(rhdf5)
library(tram)
library(etram)
library(readr)

# Paths ------------------------------------------------------------------

path <- "~/../data/UTKFace/UTKFace.h5"
out_dir <- "experiments/results/DE/UTKFace/"

# Params ------------------------------------------------------------------

fml_cond <- age_group ~ gender + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10
fml_uncond <- age_group ~ 1
decay <- 1e-4
lr <- 0.01
epochs <- 60
spl <- 6

# Load data ---------------------------------------------------------------

dat <- load_data("utkface", path = path, im_path = im_path)
tab_dat <- dat$tab_dat

## get folds
ridx <- get_ridx(out_dir, "utkface")

# SI (RPS loss) -----------------------------------------------------------

fname <- "utkface_si_rps"

fit_ref_lossrps(mod = "si", fml = fml_uncond, tab_dat = tab_dat,
                lr = lr, decay = decay, epochs = epochs,
                ridx = ridx, splits = spl, out_dir = out_dir, fname = fname)


## SILS (RPS loss) --------------------------------------------------------

fname <- "utkface_sils_rps"

fit_ref_lossrps(mod = "sils", fml = fml_cond, tab_dat = tab_dat,
                lr = lr, decay = decay, epochs = epochs,
                ridx = ridx, splits = spl, out_dir = out_dir, fname = fname)


