
#' Plot evaluation metrics
#' @param met output of \code{\link{get_metrics_allspl}}.
#' @param indiv_met output of \code{\link{get_indiv_metrics_allspl}}.
#' If not \code{NULL} rugs of performance of single ensemble members will be added to the plot.
#' @param metrics which metrics should be plotted.
#' @param facet_labs labels used for facets.
#' @export
plot_metrics <- function(met, indiv_met = NULL,
                         metrics = c("all", "acc", "eacc", "nll", "binnll", "rps", "cint",
                                     "cslope", "brier", "auc", "eauc", "qwk"),
                         facet_labs = NULL, ncol = NULL, nrow = NULL, ...) {
  all_metrics <- c("nll", "rps", "acc", "eacc", "cint",
                   "cslope", "brier", "binnll", "auc", "eauc", "qwk")
  met$metric <- factor(met$metric, levels = all_metrics[all_metrics %in% metrics])
  if (!is.null(indiv_met)) {
    indiv_met$metric <- factor(indiv_met$metric, levels = all_metrics[all_metrics %in% metrics])
  }
  if (!("all" %in% metrics)) {
    met <- met[met$metric %in% metrics, ]
    if (!is.null(indiv_met)) {
      indiv_met <- indiv_met[indiv_met$metric %in% metrics, ]
    }
  }
  if (!is.null(facet_labs)) {
    levels(met$metric)[levels(met$metric) %in% metrics] <- facet_labs
    if (!is.null(indiv_met)) {
      levels(indiv_met$metric)[levels(indiv_met$metric) %in% metrics] <- facet_labs
    }
  }
  avgs <- met %>% group_by(metric, method, topn) %>% summarise(val = mean(val))
  getPalette <- colorRampPalette(rev(brewer.pal(8, name = "Dark2")))
  p <- ggplot(met, aes(x = topn, y = val, color = method, group = interaction(spl, method))) +
    facet_wrap(~ metric, scales = "free", nrow = nrow, ncol = ncol) +
    geom_line(alpha = 0.15, size = 0.5) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_line(aes(x = topn, y = val, color = method), data = avgs, inherit.aes = FALSE) +
    geom_point(aes(x = topn, y = val, color = method), data = avgs, inherit.aes = FALSE) +
    labs(x = "Top-K ensemble", y = "Test performance") +
    scale_color_manual(name = "Ensemble", labels = c("avg" = "AVG",
                                                     "linear" = "LIN",
                                                     "log-linear" = "LOG",
                                                     "trafo" = "TRF"),
                       values = getPalette(4)) +
    theme_bw() +
    theme(...)
  if (!is.null(indiv_met)) {
    p <- p + geom_rug(aes(y = val), data = indiv_met, inherit.aes = FALSE)
  }
  p
}

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
