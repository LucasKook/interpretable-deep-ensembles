
#' Plot histories
#' @param hists output of \code{\link{load_hists}}.
#' @export
plot_hists <- function(hists, ylab = "", title = "",
                       nrow = NULL, ncol = NULL,
                       ylim = NULL, ...) {
  ggplot(hists, aes(x = epoch, y = loss,
                    color = Ens, linetype = Type, group = interaction(Type, Ens))) +
    facet_wrap(~ Split, labeller = label_both, nrow = nrow, ncol = ncol) +
    geom_line(size = 0.4) +
    xlab("Epochs") +
    ylab(ylab) +
    ggtitle(title) +
    coord_cartesian(ylim = ylim) +
    scale_linetype_manual(values = c("Train" = "dashed",
                                     "Val" = "solid")) +
    # scale_colour_discrete(guide = "none") +
    theme_bw() +
    theme(...)
}
