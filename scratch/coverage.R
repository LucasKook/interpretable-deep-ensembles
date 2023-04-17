### Algorithmic uncertainty and coverage
### Lucas Kook, March 2023

set.seed(24101968)

# Dependencies ------------------------------------------------------------

library("deeptrafo")
library("tidyverse")

# Params ------------------------------------------------------------------

nep <- 0
nens <- 5
pat <- 50

# Run ---------------------------------------------------------------------

dd <- read.table("benchmark/benchmark_data/yacht.data")
df <- list(Y = dd[, ncol(dd)], X = as.matrix(dd[, -ncol(dd)]))

nn <- \(x) x %>%
  layer_dense(input_shape = ncol(df$X), units = 16L, activation = "relu") %>%
  layer_dense(units = 16L, activation = "relu") %>%
  layer_dense(units = 8L, activation = "relu") %>%
  layer_dense(1L)

ens <- trafoensemble(
  Y ~ 0 + nn(X), data = df, list_of_deep_models = list(nn = nn), epochs = nep,
  callbacks = list(callback_early_stopping(patience = pat, restore_best_weights = TRUE)),
  n_ensemble = nens, seed = sample.int(1e6, nens)
)

ngr <- 1e3
nd <- list(
  Y = seq(min(df$Y), max(df$Y), length.out = ngr),
  X = matrix(0, ncol = ncol(df$X), nrow = ngr)
)

preds <- do.call("cbind", predict(ens, type = "trafo", newdata = nd))
ens <- apply(preds, 1, mean)
sd <- apply(preds, 1, sd)

qalp <- qnorm(0.975) # qt(0.975, df = nens)
lwr <- ens - qalp * sd
upr <- ens + qalp * sd

plot(nd$Y, ens, type = "l", lwd = 3)
lines(nd$Y, lwr, type = "l", lwd = 3, lty = 2)
lines(nd$Y, upr, type = "l", lwd = 3, lty = 2)

nens <- trafoensemble(
  Y ~ 0 + nn(X), data = df, list_of_deep_models = list(nn = nn), epochs = nep,
  callbacks = list(callback_early_stopping(patience = pat, restore_best_weights = TRUE)),
  n_ensemble = nens, seed = sample.int(1e6, nens)
)

cover <- function(x) {
  mean(lwr <= x & x <= upr)
}

npreds <- predict(nens, type = "trafo", newdata = nd)
mean(unlist(lapply(npreds, cover)))

matlines(nd$Y, do.call("cbind", npreds))
