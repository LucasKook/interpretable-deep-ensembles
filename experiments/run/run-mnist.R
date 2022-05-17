# Run experiments (mnist)
# Andrea Goetschi
# April 2022

# Reproducibility ---------------------------------------------------------

set.seed(93746)

# Deps --------------------------------------------------------------------

library(etram)
library(tram)

# Paths ------------------------------------------------------------------

out_dir <- "experiments/results/DE/MNIST/"

# Params ------------------------------------------------------------------

bs <- 512
lr <- 10^-5
epochs <- 30
spl <- 6
ens <- 5

fml <- y ~ 1

# Data --------------------------------------------------------------------

dat <- load_data("mnist")
tab_dat <- dat$tab_dat
im <- dat$im

# save_ridx(nsplts = spl, prptest = 0.1, prpval = 0.1, tab_dat = tab_dat,
#           fml = fml, out_dir = out_dir, fname = "mnist")

ridx <- get_ridx(out_dir, "mnist")

# Run experiments ---------------------------------------------------------

fname <- "mnist_ci_lossnll_wsno_augno"

ensemble(mod = "ci",
         fml = fml, tab_dat = tab_dat, im = im, ridx = ridx,
         splits = spl, ensembles = ens, nn = cnn_mnist, input_shape = dim(im)[2:4],
         bs = bs, lr = lr, epochs = epochs,
         loss = "nll",
         ws = FALSE, augment = FALSE, aug_params = NULL,
         out_dir = out_dir, fname = fname)

fname <- "mnist_ci_lossrps_wsno_augno"

ensemble(mod = "ci",
         fml = fml, tab_dat = tab_dat, im = im, ridx = ridx,
         splits = spl, ensembles = ens, nn = cnn_mnist, input_shape = dim(im)[2:4],
         bs = bs, lr = lr, epochs = epochs,
         loss = "rps",
         ws = FALSE, augment = FALSE, aug_params = NULL,
         out_dir = out_dir, fname = fname)

