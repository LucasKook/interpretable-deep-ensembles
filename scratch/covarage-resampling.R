# Produce Figure 4
# Lucas Kook
# Feb 2022

set.seed(123)

# Dependencies ------------------------------------------------------------

library(tram)
library(etram)
library(tidyverse)
library(patchwork)
theme_set(
  theme_bw() +
    theme(legend.position = "top", text = element_text(size = 13))
)

# Params ------------------------------------------------------------------

tn <- 150
n_mods <- 100
xvals <- c(0, 2)
cols <- colorspace::qualitative_hcl(n = 3, l = 40)
names(cols) <- c("LIN-Ens", "LOG-Ens", "TRF-Ens")
col2 <- colorspace::diverge_hcl(n = 2, l = 40)

# FUNs --------------------------------------------------------------------

missp <- function(n = tn) {
  x <- runif(n, min = -4, max = 4)
  y <- 3 * plogis(3 * x) - rnorm(n, sd = 0.6^2)
  y <- y - mean(y)
  data.frame(y = y, x = x)
}

get_trafos <- function(dgp) {
  tdat <- dgp()
  mods <- lapply(1:n_mods, function(idx) as.mlt(BoxCox(y ~ x, data = dgp(),
                                                       support = range(tdat$y))))

  # Classical ensemble trafo
  lcdf <- lapply(mods, predict, newdata = tdat, type = "distribution")
  CE <- qnorm(get_ensemble(lcdf, "linear"))
  nd <- data.frame(expand.grid(y = seq(min(tdat$y), max(tdat$y), length.out = 1e3), x = xvals))
  trafos <- do.call("cbind", lapply(lapply(mods, predict, newdata = nd), as.double))

  ndd <- data.frame(x = seq(min(tdat$x), max(tdat$x), length.out = 1e3))
  cpreds <- t(do.call("rbind", lapply(lapply(mods, predict, newdata = ndd,
                                             type = "quantile", p = 0.5), as.double)))

  # Trafo ensemble
  mparm <- colMeans(do.call("rbind", lapply(mods, coef)))
  TE <- do.call("cbind", lapply(lapply(mods, predict, type = "trafo", newdata = tdat), as.double))

  list(
    td = tdat,
    df = data.frame(obs = 1:nrow(tdat), CE = CE, TE = TE) %>%
      gather("member", "trafo", paste0("TE.", 1:n_mods)),
    tf = data.frame(obs = 1:nrow(trafos), y = nd$y, x= nd$x, trafo = trafos),
    cp = data.frame(obs = 1:nrow(cpreds), x = ndd$x, cps = unname(cpreds))
  )
}

pdat <- get_trafos(missp)

udat <- pdat$tf %>%
  gather("mem", "val", paste0("trafo.", 1:n_mods)) %>%
  group_by(obs, x) %>%
  mutate(avg = mean(val), sd = sd(val))

cdat <- udat %>% ungroup() %>% mutate(cover = (val > avg - 1.96 * sd & val < avg + 1.96 * sd)) %>%
  group_by(x, y) %>% summarize(cover = mean(cover))

cdat %>% summarize(mean(cover))

ggplot(cdat, aes(x = y, y = cover, color = factor(x))) +
  geom_step()

p1 <- ggplot(pdat$df, aes(x = CE, y = trafo, col = member)) +
  geom_point(alpha = 0.3, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  labs(
    x = expression(F[Z]^{-1}*(bar(F)[M]^c*(y[i]*'|'*x[i]))),
    # y = expression(F[Z]^{-1}*(F[m](y[i]*'|'*x[i])))
    y = expression(h[m](y[i]*'|'*x[i]))
  )

p2 <- ggplot(udat, aes(x = y, y = pnorm(val), group = interaction(mem, x))) +
  geom_ribbon(aes(x = y, ymin = pnorm(avg - 1.96 * sd), ymax = pnorm(avg + 1.96 * sd),
                  group = x, fill = factor(x)), show.legend = FALSE,
              data = udat, inherit.aes = FALSE, alpha = 0.3) +
  geom_line(aes(linetype = "Member", color = factor(x)), alpha = 0.7,
            show.legend = FALSE) +
  geom_line(aes(x = y, y = pnorm(avg), linetype = "TRF-Ens", group = x,
                color = factor(x)), data = udat, show.legend = FALSE,
            inherit.aes = FALSE, lwd = 0.8) +
  labs(x = "y", y = "CDF of Y | X = x", linetype = "Conditional CDF") +
  theme(legend.position = c(0.25, 0.8),
        legend.background = element_rect(fill = "transparent")) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_color_manual(values = col2) +
  scale_fill_manual(values = col2)

cdat <- pdat$cp %>%
  gather("mem", "val", paste0("cps.", 1:n_mods)) %>%
  group_by(obs) %>%
  mutate(avg = mean(val), sd = sd(val))

p3 <- ggplot(cdat, aes(x = x, y = val, group = mem)) +
  geom_ribbon(aes(ymin = avg - 1.96 * sd, ymax = avg + 1.96 * sd), alpha = 0.3, fill = "gray") +
  geom_line(aes(linetype = "Member"), alpha = 0.7) +
  geom_line(aes(y = avg, linetype = "TRF-Ens")) +
  geom_point(aes(x = x, y = y), data = pdat$td, inherit.aes = FALSE, alpha = 0.3) +
  labs(y = "y", color = "", linetype = "Conditional median", fill = "") +
  theme(legend.position = c(0.3, 0.8),
        legend.background = element_rect(fill = "transparent")) +
  scale_fill_manual(values = "gray") +
  scale_linetype_manual(values = c(2, 1)) +
  geom_vline(xintercept = xvals, color = col2, lty = 3)

(p1 + labs(tag = "A")) + (p3 + labs(tag = "B")) + (p2 + labs(tag = "C"))
