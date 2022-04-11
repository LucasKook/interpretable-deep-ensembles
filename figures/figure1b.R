# Transformation model illustration
# Ordered outcomes
# LK Apr 2022

# Params, FUNs
oy <- seq(0.1, 0.9, length.out = 7)
ncl <- length(oy)
y <- seq(0 + 1e-2, 1 - 0.15, length.out = res <- 1e3)
alp <- pi
bet <- pi
pY <- function(y) plogis(1.1 * qlogis(y) - 0.3) # pbeta(y, alp, bet)
dY <- function(y) dbeta(y, alp, bet)
pZ <- function(z) 1 - exp(-exp(z)) # plogis
dZ <- dlogis
qY <- function(p) qbeta(p, alp, bet)
qZ <- qlogis

h <- function(y) qZ(pY(y))
hp <- function(y) {
 tmp <- (1 / dZ(qZ(pY(y)))) * dY(y)
 tmp <- ifelse(is.nan(tmp), 0, tmp)
}

opar <- par(no.readonly = TRUE)

par(opar)

# Plot
pdf("figure1b.pdf", width = 6, height = 5)
par(mar = c(0.1, 0.1, 0.1, 0.1) + 0.1)
layout(matrix(c(1,2,2,2,
                3,4,4,4,
                3,4,4,4,
                3,4,4,4), nrow = 4, ncol = 4, byrow = TRUE))

cols <- colorspace::diverge_hcl(length(oy))
tcx <- 1
tcxa <- 1.4

plot.new()
plot(oy, c(pZ(h(oy))[-ncl], 1), type = "s", axes = FALSE, xlim = c(1, 0),
     ylim = c(0, 1), col = rgb(.1, .1, .1, .5), pch = 20, cex = 2)
mtext(expression(F[Y](y~'|'~x)), 2, line = 3, cex = tcx, adj = 1)
axis(1, at = oy, labels = parse(text = paste0("y[" , seq_along(oy), "]")),
     cex.axis = tcxa)
axis(2, las = 1)
arrows(oy, 0, oy, c(pZ(h(oy))[-ncl], 1), length = 0, col = cols, lwd = 2)
plot(pZ(h(y)), h(y), type = "l", axes = FALSE, xlim = c(1, 0),
     ylim = range(h(oy)) + c(-0.5, 0), col = "white")
axis(4, las = 1, labels = -3:2, at = -3:2)
axis(3)
mtext(expression(F[Z](h(y~'|'~x))), 3, line = 2, cex = tcx, adj = 0)
mtext(expression(h(y~'|'~x)), 4, line = 3, cex = tcx, las = 1)
arrows(0, h(oy)[-ncl], pZ(h(oy))[-ncl], h(oy)[-ncl], length = 0, col = cols,
       lwd = 2)
lines(pZ(h(y)), h(y), lwd = 1.2, col = rgb(.1, .1, .1, .5))
plot(oy, c(h(oy)[-ncl], Inf), type = "p", axes = FALSE, xlim = c(1, 0),
     col = cols, pch = 20, cex = 2, ylim = range(h(oy)) + c(-0.5, 0))
par(opar)
dev.off()
