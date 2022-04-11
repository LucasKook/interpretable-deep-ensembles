# Transformation model illustration
# Continuous outcomes
# LK Apr 2022

# Params, FUNs
y <- seq(0.1, 11, length.out = res <- 1e4)
pY <- function(y) pchisq(y, df = 4)
dY <- function(y) dchisq(y, df = 4)
pZ <- pnorm
dZ <- dnorm
qY <- function(p) qchisq(p, df = 4)
qZ <- qnorm

h <- function(y) qZ(pY(y))
hp <- function(y) {
 tmp <- (1 / dZ(qZ(pY(y)))) * dY(y)
 tmp <- ifelse(is.nan(tmp), 0, tmp)
}

plot(y, pY(y), type = "l")
plot(y, h(y), type = "l")
plot(h(y), pZ(h(y)), type = "l")

opar <- par(no.readonly = TRUE)

# Plot
pdf("figure1a.pdf", width = 6, height = 5)
par(mar = c(0.1, 0.1, 0.1, 0.1) + 0.1)
layout(matrix(c(1,2,2,2,
                3,4,4,4,
                3,4,4,4,
                3,4,4,4), nrow = 4, ncol = 4, byrow = TRUE))

cols <- colorRampPalette(c("cornflowerblue", "orange"))(res)
tcx <- 1

plot.new()
plot(y, pZ(h(y)), type = "l", axes = FALSE, xlim = rev(range(y)),
     ylim = c(0, 1))
mtext(expression(F[Y](y~'|'~x)), 2, line = 3, cex = tcx, adj = 1)
mtext(expression(y), 1, line = 3, cex = tcx)
axis(1)
axis(2, las = 1)
arrows(y, 0, y, pZ(h(y)), length = 0, col = cols)
plot(pZ(h(y)), h(y), type = "l", axes = FALSE, xlim = c(1, 0))
axis(4, las = 1)
axis(3)
mtext(expression(F[Z](h(y~'|'~x))), 3, line = 2, cex = tcx, adj = 0)
mtext(expression(h(y~'|'~x)), 4, line = 3, cex = tcx, las = 1)
arrows(0, h(y), pZ(h(y)), h(y), length = 0, col = cols)
plot(y, h(y), type = "p", axes = FALSE, xlim = rev(range(y)), col = cols, pch = 20, cex = 0.1)
par(opar)
dev.off()
