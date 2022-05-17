# Run experiments (utkface)
# Andrea Goetschi
# April 2022

# Command line args -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (identical(args, character(0))) {
  loss <- "rps"
  fml <- age_group ~ 1
  # fml <- age_group ~ gender + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 +
  #          x_9 + x_10
  mod <- "ci"
  fname <- paste0("utkface_", mod, "_loss", loss, "_wsyes_augno")
} else {
  loss <- args[1]
  fml <- as.formula(args[2])
  mod <- args[3]
  fname <- paste0("utkface_", mod, "_loss", loss, "_wsyes_augno")
}

# Reproducibility ---------------------------------------------------------

set.seed(738593)

# Deps --------------------------------------------------------------------

library(rhdf5)
library(tram)
library(etram)

# Paths ------------------------------------------------------------------

path <- "~/../data/UTKFace/UTKFace.h5"
out_dir <- "experiments/results/DE/UTKFace/"

# Params ------------------------------------------------------------------

bs <- 32
lr <- 10^-4
epochs <- 150
spl <- 6
ens <- 5

# Data --------------------------------------------------------------------

dat <- load_data("utkface", path = path)
tab_dat <- dat$tab_dat
im <- dat$im

# fml_cond <- age_group ~ gender + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10
# save_ridx(nsplts = spl, prptest = 0.2, prpval = 0.2, tab_dat = tab_dat,
#           fml = fml_cond, out_dir = out_dir, fname = "utkface")

## get folds
ridx <- get_ridx(out_dir, "utkface")

# Run experiments ---------------------------------------------------------

ensemble(mod = mod,
         fml = fml, tab_dat = tab_dat, im = im, ridx = ridx,
         splits = spl, ensembles = ens,
         nn = cnn_utkface, input_shape = dim(im)[2:4],
         bs = bs, lr = lr, epochs = epochs,
         loss = loss,
         ws = TRUE, augment = FALSE,
         train_batchwise = FALSE,
         out_dir = out_dir, fname = fname)
