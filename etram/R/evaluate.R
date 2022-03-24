
#' Evaluation of ensemble
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @param y_true one-hot encoded true \code{y}.
#' @param type ensemble technique.
#' @param metrics metrics to evaluate.
#' @param topk logical. Whether to evaluate top-k ensembles from \code{1:K} based on validation loss.
#' @param order_metric metric used to evaluate and order validation CDFs.
#' @param lys_cdf_val list of CDFs.
#' @param y_true_val one-hot encoded true \code{y} of validation set.
#' @param cutoff cutoff used to evaluate binary metrics.
#' @param p weighting scheme used for qwk.
#' @param weights vector of weights for ensemble members.
#' @export
get_metrics <- function(lys_cdf, y_true, type = c("all", "linear", "log-linear", "trafo", "avg"),
                        metrics = c("all", "acc", "binacc", "eacc", "ebinacc", "nll", "binnll", "rps", "cint",
                                    "cslope", "brier", "auc", "eauc", "qwk", "eqwk"),
                        topk = TRUE, order_metric = c("nll", "rps"),
                        lys_cdf_val = NULL, y_true_val = NULL, cutoff = 3, p = 2,
                        weights = rep(1, length(lys_cdf))) {
  all_metrics <- c("acc", "binacc", "eacc", "ebinacc", "nll", "binnll", "rps", "cint",
                   "cslope", "brier", "auc", "eauc", "qwk", "eqwk")
  all_methods <- c("linear", "log-linear", "trafo", "avg")
  nmetrics <- length(all_metrics)
  nmethods <- length(all_methods)
  type <- match.arg(type)
  # metrics <- match.arg(metrics)
  order_metric <- match.arg(order_metric)
  ens <- length(lys_cdf)
  lys_cdf <- lapply(lys_cdf, as.matrix)
  if (topk) {
    ranks <- get_order(lys_cdf_val, y_true_val, order_metric)
    lys_cdf <- lys_cdf[ranks]
    ret <- data.frame(metric = rep(rep(all_metrics, ens), nmethods),
                      topn = rep(rep(1:ens, each = nmetrics), nmethods),
                      val = rep(numeric(ens), nmethods),
                      method = rep(all_methods, each = nmetrics * ens))
    for (meth in all_methods) {
      if (!(meth %in% "avg")) {
        for (n in 1:ens) {
          sub_lys_cdf <- lys_cdf[1:n]
          cdf <- get_ensemble(lys_cdf = sub_lys_cdf, type = meth, weights = weights[1:n])
          ret[ret$method == meth & ret$topn == n, "val"] <- c(get_acc(cdf, y_true),
                                                              get_binacc(cdf, y_true, cutoff = cutoff),
                                                              1 - get_acc(cdf, y_true),
                                                              1 - get_binacc(cdf, y_true, cutoff = cutoff),
                                                              get_nll(cdf, y_true),
                                                              get_binnll(cdf, y_true, cutoff = cutoff),
                                                              get_rps(cdf, y_true),
                                                              get_cal(cdf, y_true)$cint,
                                                              get_cal(cdf, y_true)$cslope,
                                                              get_brier(cdf, y_true, cutoff = cutoff),
                                                              get_auc(cdf, y_true, cutoff = cutoff),
                                                              1 - get_auc(cdf, y_true, cutoff = cutoff),
                                                              get_qwk(cdf, y_true, p = p),
                                                              1 - get_qwk(cdf, y_true, p = p))
        }
      } else if (meth %in% "avg") {
        for (n in 1:ens) {
          sub_lys_cdf <- lys_cdf[1:n]
          ret[ret$method == meth & ret$topn == n, "val"] <- c(get_avg_acc(sub_lys_cdf, y_true),
                                                              get_avg_binacc(sub_lys_cdf, y_true, cutoff = cutoff),
                                                              1 - get_avg_acc(sub_lys_cdf, y_true),
                                                              1 - get_avg_binacc(sub_lys_cdf, y_true, cutoff = cutoff),
                                                              get_avg_nll(sub_lys_cdf, y_true),
                                                              get_avg_binnll(sub_lys_cdf, y_true, cutoff = cutoff),
                                                              get_avg_rps(sub_lys_cdf, y_true),
                                                              get_avg_cal(sub_lys_cdf, y_true)$cint,
                                                              get_avg_cal(sub_lys_cdf, y_true)$cslope,
                                                              get_avg_brier(sub_lys_cdf, y_true),
                                                              get_avg_auc(sub_lys_cdf, y_true, cutoff = cutoff),
                                                              1 - get_avg_auc(sub_lys_cdf, y_true, cutoff = cutoff),
                                                              get_avg_qwk(sub_lys_cdf, y_true, p = p),
                                                              1 - get_avg_qwk(sub_lys_cdf, y_true, p = p))
        }
      }
    }
  } else if (!topk) {
    ret <- data.frame(metric = rep(all_metrics, nmethods),
                      val = numeric(nmetrics * nmethods),
                      method = rep(all_methods, each = nmetrics))
    for (meth in all_methods) {
      if (!(meth %in% "avg")) {
        cdf <- get_ensemble(lys_cdf = lys_cdf, type = meth, weights = weights)
        ret[ret$method == meth, "val"] <- c(get_acc(cdf, y_true),
                                            get_binacc(cdf, y_true, cutoff = cutoff),
                                            1 - get_acc(cdf, y_true),
                                            1 - get_binacc(cdf, y_true, cutoff = cutoff),
                                            get_nll(cdf, y_true),
                                            get_binnll(cdf, y_true, cutoff = cutoff),
                                            get_rps(cdf, y_true),
                                            get_cal(cdf, y_true)$cint,
                                            get_cal(cdf, y_true)$cslope,
                                            get_brier(cdf, y_true, cutoff = cutoff),
                                            get_auc(cdf, y_true, cutoff = cutoff),
                                            1 - get_auc(cdf, y_true, cutoff = cutoff),
                                            get_qwk(cdf, y_true, p = p),
                                            1 - get_qwk(cdf, y_true, p = p))
      } else if (meth %in% "avg") {
        ret[ret$method == meth, "val"] <- c(get_avg_acc(lys_cdf, y_true),
                                            get_avg_binacc(lys_cdf, y_true, cutoff = cutoff),
                                            1 - get_avg_acc(lys_cdf, y_true),
                                            1 - get_avg_binacc(lys_cdf, y_true, cutoff = cutoff),
                                            get_avg_nll(lys_cdf, y_true),
                                            get_avg_binnll(lys_cdf, y_true, cutoff = cutoff),
                                            get_avg_rps(lys_cdf, y_true),
                                            get_avg_cal(lys_cdf, y_true)$cint,
                                            get_avg_cal(lys_cdf, y_true)$cslope,
                                            get_avg_brier(lys_cdf, y_true),
                                            get_avg_auc(lys_cdf, y_true, cutoff = cutoff),
                                            1 - get_avg_auc(lys_cdf, y_true, cutoff = cutoff),
                                            get_avg_qwk(lys_cdf, y_true, p = p),
                                            1 - get_avg_qwk(lys_cdf, y_true, p = p))
      }
    }
  }
  switch(
    type,
    "all" = (if ("all" %in% metrics) ret
             else ret[ret$metric %in% metrics, ]),
    "linear" = (if ("all" %in% metrics) ret[ret$method == "linear", ]
                else ret[ret$metric %in% metrics & ret$method == "linear", ]),
    "log-linear" = (if("all" %in% metrics) ret[ret$method == "log-linear", ]
                    else ret[ret$metric %in% metrics & ret$method == "log-linear", ]),
    "trafo" = (if ("all" %in% metrics) ret[ret$method == "trafo", ]
               else ret[ret$metric %in% metrics & ret$method == "trafo", ]),
    "avg" = (if ("all" %in% metrics) ret[ret$method == "avg", ]
             else ret[ret$metric %in% metrics & ret$method == "avg", ])
  )
}

#' Evaluation of all splits
#' @param lys_cdf_all list of all CDFs (e.g. output of \code{\link{list_cdfs}})
#' @param y_true_all list of one-hot encoded true \code{y} for all splits.
#' @param type ensemble technique.
#' @param metrics metrics to evaluate.
#' @param topk logical. Whether to evaluate top-k ensembles from \code{1:K} based on validation loss.
#' @param order_metric metric used to evaluate and order validation CDFs.
#' @param lys_cdf_val_all list of all CDFs.
#' @param y_true_val_all list of one-hot encoded true \code{y} of validation set for all splits.
#' @param cutoff cutoff used to evaluate binary metrics.
#' @param p weighting scheme used for qwk.
#' @param weights matrix of weights (cols) per split (rows).
#' @export
get_metrics_allspl <- function(lys_cdf_all, y_true_all, type = c("all", "linear", "log-linear", "trafo", "avg"),
                               metrics = c("all", "acc", "binacc", "eacc", "ebinacc", "nll", "binnll", "rps", "cint",
                                           "cslope", "brier", "auc", "eauc", "qwk", "eqwk"),
                               topk = TRUE, order_metric = c("nll", "rps"),
                               lys_cdf_val_all = NULL, y_true_val_all = NULL, cutoff = 3, p = 2,
                               weights = matrix(1, nrow = length(lys_cdf_all), ncol = length(lys_cdf_all[[1]]))) {
  spl <- length(lys_cdf_all)
  ret <- list()
  for (s in seq_len(spl)) {
    lys_cdf <- lys_cdf_all[[s]]
    y_true <- y_true_all[[s]]
    if (!is.null(lys_cdf_val_all)) {
      lys_cdf_val <- lys_cdf_val_all[[s]]
      y_true_val <- y_true_val_all[[s]]
    } else {
      lys_cdf_val <- NULL
      y_true_val <- NULL
    }
    w <- weights[s, ]
    tmp <- get_metrics(lys_cdf = lys_cdf, y_true = y_true, type = type,
                       metrics = metrics,
                       topk = topk, order_metric = order_metric,
                       lys_cdf_val = lys_cdf_val, y_true_val = y_true_val, cutoff = cutoff, p = p,
                       weights = w)
    tmp$spl <- s
    ret[[s]] <- tmp
  }
  do.call("rbind", ret)
}

#' Evaluation of single ensemble members
#' @param lys_cdf list of CDFs (e.g. CDFs of ensemble members).
#' @param y_true one-hot encoded true \code{y}.
#' @param metrics metrics to evaluate.
#' @param cutoff cutoff used to evaluate binary metrics.
#' @param p weighting scheme used for qwk.
#' @export
get_indiv_metrics <- function(lys_cdf, y_true,
                              metrics = c("all", "acc", "binacc", "eacc", "ebinacc", "nll", "binnll", "rps", "cint",
                                          "cslope", "brier", "auc", "eauc", "qwk", "eqwk"),
                              cutoff = 3, p = 2) {
  # metrics <- match.arg(metrics)
  all_metrics <- c("acc", "binacc", "eacc", "ebinacc", "nll", "binnll", "rps", "cint",
                   "cslope", "brier", "auc", "eauc", "qwk", "eqwk")
  nmetrics <- length(all_metrics)
  ens <- length(lys_cdf)
  ret <- data.frame(metric = rep(all_metrics, each = ens),
                    ens = rep(1:ens, nmetrics),
                    val = c(do.call("rbind", lapply(lys_cdf, get_acc, y_true = y_true)),
                            do.call("rbind", lapply(lys_cdf, get_binacc, y_true = y_true, cutoff = cutoff)),
                            do.call("rbind", lapply(lys_cdf, function(x) 1 - get_acc(cdf = x, y_true = y_true))),
                            do.call("rbind", lapply(lys_cdf, function(x) 1 - get_binacc(cdf = x, y_true = y_true, cutoff = cutoff))),
                            do.call("rbind", lapply(lys_cdf, get_nll, y_true = y_true)),
                            do.call("rbind", lapply(lys_cdf, get_binnll, y_true = y_true, cutoff = cutoff)),
                            do.call("rbind", lapply(lys_cdf, get_rps, y_true = y_true)),
                            do.call("rbind", lapply(lys_cdf, function(x) get_cal(cdf = x, y_true = y_true)$cint)),
                            do.call("rbind", lapply(lys_cdf, function(x) get_cal(cdf = x, y_true = y_true)$cslope)),
                            do.call("rbind", lapply(lys_cdf, get_brier, y_true = y_true, cutoff = cutoff)),
                            do.call("rbind", lapply(lys_cdf, get_auc, y_true = y_true, cutoff = cutoff)),
                            do.call("rbind", lapply(lys_cdf, function(x) 1 - get_auc(cdf = x, y_true = y_true, cutoff = cutoff))),
                            do.call("rbind", lapply(lys_cdf, get_qwk, y_true = y_true, p = p)),
                            do.call("rbind", lapply(lys_cdf, function(x) 1 - get_qwk(cdf = x, y_true = y_true, p = p)))
                    ))
  ret$metric <- factor(ret$metric, levels = c("nll", "rps", "qwk", "eqwk", "acc", "eacc",
                                              "binnll", "binacc", "ebinacc", "auc", "eauc", "brier",
                                              "cint", "cslope"))
  if ("all" %in% metrics)
    return(ret)
  else
    return(ret[ret$metric %in% metrics, ])
}

#' Evaluation of single ensemble members of all splits
#' @param lys_cdf_all list of all CDFs (e.g. output of \code{\link{list_cdfs}})
#' @param y_true_all list of one-hot encoded true \code{y} for all splits.
#' @param metrics metrics to evaluate.
#' @param cutoff cutoff used to evaluate binary metrics.
#' @param p weighting scheme used for qwk.
#' @export
get_indiv_metrics_allspl <- function(lys_cdf_all, y_true_all,
                                     metrics = c("all", "acc", "binacc", "eacc", "ebinacc", "nll", "binnll", "rps", "cint",
                                                 "cslope", "brier", "auc", "eauc", "qwk", "eqwk"),
                                     cutoff = 3, p = 2) {
  spl <- length(lys_cdf_all)
  ret <- list()
  for (s in seq_len(spl)) {
    lys_cdf <- lys_cdf_all[[s]]
    y_true <- y_true_all[[s]]
    tmp <- get_indiv_metrics(lys_cdf = lys_cdf, y_true = y_true,
                             metrics = metrics,
                             cutoff = cutoff, p = p)
    tmp$spl <- s
    ret[[s]] <- tmp
  }
  do.call("rbind", ret)
}
