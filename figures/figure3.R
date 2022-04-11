# Produce Figure 3
# Lucas Kook
# Feb 2022

set.seed(1)

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(patchwork)
theme_set(theme_bw() + theme(legend.position = "top"))

# FUNs --------------------------------------------------------------------

h <- function(y, theta) y %*% theta
pZ <- function(z) plogis(z)
dZ <- function(z) dlogis(z)
hp <- function(dy, theta) dy %*% theta
cols <- colorspace::qualitative_hcl(n = 3, l = 40)
names(cols) <- c("LIN", "LOG", "TRF")

# Data --------------------------------------------------------------------

res <- 1e3
ys <- seq(-8, 8, length.out = res)

thetas <- matrix(rnorm(2 * 5), ncol = 2)
thetas[, 2] <- log(1 + exp(thetas[, 2]))

y <- cbind(1, ys)
dy <- cbind(0, rep(1, length(ys)))

probs <- pZ(h(y, t(thetas)))
dens <- dZ(h(y, t(thetas))) * hp(dy, t(thetas))

# Linear ensemble
pcensemble <- apply(probs, 1, mean)
censemble <- apply(dens, 1, mean)

# Transformation ensemble
mth <- apply(thetas, 2, mean)
iensemble <- c(dZ(h(y, mth)) * hp(dy, mth))
piensemble <- c(pZ(h(y, mth)))

# Mean NLL
mNLL <- exp(rowMeans(log(dens)))

# log-linear ensemble
lensemble <- mNLL / sum(mNLL * diff(ys)[1])
plensemble <- cumsum(mNLL) / max(cumsum(mNLL))

edat <- data.frame(
  y = ys,
  LIN = censemble,
  TRF = iensemble,
  LOG = lensemble,
  AVG = mNLL
) %>%
  gather("method", "value", LIN, TRF, LOG)

adat <- data.frame(y = ys, value = mNLL)

mdat <- data.frame(
  y = ys,
  dens = dens
) %>%
  gather("member", "value", dens.1:dens.5)

p1 <- ggplot(edat,
       aes(x = y, y = value, color = method, group = method)) +
  geom_line(aes(x = y, y = value, group = member), data = mdat,
            inherit.aes = FALSE, alpha = 0.3, linetype = 2) +
  geom_line() +
  labs(x = "y", y = "density", color = "Ensemble") +
  scale_color_manual(values = cols)

p2 <- ggplot(edat, aes(x = y, y = -log(value), color = method, group = method)) +
  geom_line(aes(x = y, y = -log(value), group = member), data = mdat,
            inherit.aes = FALSE, alpha = 0.3, linetype = 2) +
  geom_line(aes(x = y, y = -log(value), lwd = "AVG"), data = adat,
            inherit.aes = FALSE, alpha = 1, linetype = 1) +
  geom_line(show.legend = FALSE) +
  labs(x = "y", y = "NLL", color = "Ensemble", size = "") +
  scale_size_manual(values = c("AVG" = 0.7)) +
  scale_color_manual(values = cols)

(p1 + labs(tag = "A")) + (p2 + labs(tag = "B"))
ggsave("figure3.pdf", height = 4, width = 9)
