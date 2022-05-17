# Run experiments (melanoma)
# Andrea Goetschi
# April 2022

# Command line args -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (identical(args, character(0))) {
  loss <- "rps"
  fml <- "target ~ 1"
  mod <- "ci"
  fname <- paste0("mela_", mod, "_loss", loss, "_wsyes_augno")
} else {
  loss <- args[1]
  fml <- as.formula(args[2])
  mod <- args[3]
  fname <- paste0("mela_", mod, "_loss", loss, "_wsyes_augno")
}

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
out_dir <- "experiments/results/DE/melanoma/"

# Params ------------------------------------------------------------------

bs <- 64
lr <- 10^-4
epochs <- 90
spl <- 6
ens <- 5

# Load data ---------------------------------------------------------------

dat <- load_data("melanoma", path = path, im_path = im_path)
tab_dat <- dat$tab_dat
im <- dat$im

# fml_cond <- target ~ age_s
# save_ridx(nsplts = spl, prptest = 0.2, prpval = 0.2, tab_dat = tab_dat,
#           fml = fml_cond, out_dir = out_dir, fname = "melanoma")

## get folds
ridx <- get_ridx(out_dir, "melanoma")

# Run experiments ---------------------------------------------------------

ensemble(mod = mod,
         fml = fml, tab_dat = tab_dat, im = im, ridx = ridx,
         splits = spl, ensembles = ens,
         nn = cnn_melanoma, input_shape = dim(im)[2:4],
         bs = bs, lr = lr, epochs = epochs,
         loss = loss,
         ws = TRUE, augment = FALSE,
         train_batchwise = FALSE,
         out_dir = out_dir, fname = fname)

