# Produce Figure 4
# Lucas Kook
# Feb 2022

set.seed(123)

# Dependencies ------------------------------------------------------------

library(tram)
library(etram)
library(tidyverse)
library(patchwork)
theme_set(theme_bw())

# Params ------------------------------------------------------------------

tn <- 150
n_mods <- 5
cols <- colorspace::qualitative_hcl(n = 3, l = 40)
names(cols) <- c("LIN", "LOG", "TRF")

# FUNs --------------------------------------------------------------------

missp <- function(n = tn) {
  x <- rnorm(n, sd = 2)
  h <- rnorm(n, sd = 1^2)
  y <- 2 * plogis(x) - 0 * h + rnorm(n, sd = 0.8^2)
  data.frame(y = y, x = x)
}

get_trafos <- function(dgp) {
  tdat <- dgp()
  mods <- lapply(1:n_mods, function(idx) as.mlt(BoxCox(y ~ x, data = dgp(),
                                                       support = range(tdat$y))))

  # Classical ensemble trafo
  lcdf <- lapply(mods, predict, newdata = tdat, type = "distribution")
  CE <- qnorm(get_ensemble(lcdf, "linear"))
  nd <- data.frame(y = seq(min(tdat$y), max(tdat$y), length.out = 1e3), x = 0)
  trafos <- do.call("cbind", lapply(mods, predict, newdata = nd))

  ndd <- data.frame(x = seq(min(tdat$x), max(tdat$x), length.out = 1e3))
  cpreds <- t(do.call("rbind", lapply(mods, predict, newdata = ndd, type = "quantile", p = 0.5)))

  # Trafo ensemble
  mparm <- colMeans(do.call("rbind", lapply(mods, coef)))
  TE <- do.call("cbind", lapply(mods, predict, type = "trafo", newdata = tdat))

  list(
    td = tdat,
    df = data.frame(obs = 1:nrow(tdat), CE = CE, TE = TE) %>%
      gather("member", "trafo", TE.1:TE.5),
    tf = data.frame(obs = 1:nrow(trafos), y = nd$y, trafo = trafos),
    cp = data.frame(obs = 1:nrow(cpreds), x = ndd$x, cps = unname(cpreds))
  )
}

pdat <- get_trafos(missp)

udat <- pdat$tf %>%
  gather("mem", "val", trafo.1:trafo.5) %>%
  group_by(obs) %>%
  mutate(avg = mean(val), sd = sd(val))

p1 <- ggplot(pdat$df, aes(x = CE, y = trafo, col = member)) +
  geom_point(alpha = 0.3, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  labs(
    x = expression(F[Z]^{-1}*(bar(F)[M]^c*(y[i]*'|'*x[i]))),
    y = expression(h[m](y[i]*'|'*x[i]))
  )

p2 <- ggplot(udat, aes(x = y, y = pnorm(val), group = mem)) +
  geom_ribbon(aes(x = y, ymin = pnorm(avg - 2 * sd), ymax = pnorm(avg + 2 * sd)),
              data = udat, inherit.aes = FALSE, alpha = 0.3) +
  geom_line(aes(linetype = "Individual"), alpha = 0.7) +
  geom_line(aes(x = y, y = pnorm(avg), linetype = "Ensemble"), data = udat,
            inherit.aes = FALSE) +
  labs(x = "y", y = "CDF of Y | x = 0", linetype = "") +
  theme(legend.position = c(0.25, 0.8),
        legend.background = element_rect(fill = "transparent"))

cdat <- pdat$cp %>%
  gather("mem", "val", cps.1:cps.5) %>%
  group_by(obs) %>%
  mutate(avg = mean(val), sd = sd(val))

p3 <- ggplot(cdat, aes(x = x, y = val, group = mem)) +
  geom_ribbon(aes(ymin = avg - 2 * sd, ymax = avg + 2 * sd), alpha = 0.3, fill = "gray") +
  geom_line(aes(linetype = "Individual"), alpha = 0.7) +
  geom_line(aes(y = avg, linetype = "Ensemble")) +
  geom_point(aes(x = x, y = y), data = pdat$td, inherit.aes = FALSE, alpha = 0.3) +
  labs(y = "y", color = "", linetype = "Conditional median", fill = "") +
  theme(legend.position = c(0.25, 0.8),
        legend.background = element_rect(fill = "transparent")) +
  scale_fill_manual(values = "gray")

(p1 + labs(tag = "A")) + (p3 + labs(tag = "B")) + (p2 + labs(tag = "C"))

ggsave("figure4.pdf", height = 4, width = 12)
