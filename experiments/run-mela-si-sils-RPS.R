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

im_path <- "~/../data/mela_all/data/train_mela_images/"
path <- "~/../data/mela_all/data/train.csv"
out_dir <- "~/git-repos/MasterThesis/experiments/results/DE/melanoma/"

# Params ------------------------------------------------------------------

fml_cond <- target ~ age_s
fml_uncond <- target ~ 1
decay <- 1e-4
lr <- 0.01
epochs <- 50
spl <- 6

# Load data ---------------------------------------------------------------

dat <- load_data("melanoma", path = path, im_path = im_path)
tab_dat <- dat$tab_dat

## get folds
ridx <- get_ridx(out_dir, "melanoma")

# SI (RPS loss) -----------------------------------------------------------

fname <- "mela_si_rps"

fit_ref_lossrps(mod = "si", fml = fml_uncond, tab_dat = tab_dat,
                lr = lr, decay = decay, epochs = epochs,
                ridx = ridx, splits = spl, out_dir = out_dir, fname = fname)


## SILS (RPS loss) --------------------------------------------------------

fname <- "mela_sils_rps"

fit_ref_lossrps(mod = "sils", fml = fml_cond, tab_dat = tab_dat,
                lr = lr, decay = decay, epochs = epochs,
                ridx = ridx, splits = spl, out_dir = out_dir, fname = fname)


