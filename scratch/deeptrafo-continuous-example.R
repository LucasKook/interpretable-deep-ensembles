
library("deeptrafo")

missp <- function(n = 100) {
  x <- runif(n, min = -4, max = 4)
  y <- 3 * plogis(3 * x) - rnorm(n, sd = 0.6^2)
  y <- y - mean(y)
  data.frame(y = y, x = x)
}

d <- missp(1e3)

m <- BoxCoxNN(y | x ~ 1, data = d, optimizer = optimizer_adam(1e-1, decay = 1e-4))
(ens <- ensemble(m, epochs = 1, verbose = TRUE, validation_split = 0.1))
plot(ens, type = "trafo", newdata = data.frame(x = 0))
