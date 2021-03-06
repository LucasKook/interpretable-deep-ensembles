% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate.R
\name{get_indiv_metrics_allspl}
\alias{get_indiv_metrics_allspl}
\title{Evaluation of single ensemble members of all splits}
\usage{
get_indiv_metrics_allspl(
  lys_cdf_all,
  y_true_all,
  metrics = c("all", "acc", "binacc", "eacc", "ebinacc", "nll", "binnll", "rps",
    "cint", "cslope", "brier", "auc", "eauc", "qwk", "eqwk"),
  cutoff = 3,
  p = 2
)
}
\arguments{
\item{lys_cdf_all}{list of all CDFs (e.g. output of \code{\link{list_cdfs}})}

\item{y_true_all}{list of one-hot encoded true \code{y} for all splits.}

\item{metrics}{metrics to evaluate.}

\item{cutoff}{cutoff used to evaluate binary metrics.}

\item{p}{weighting scheme used for qwk.}
}
\description{
Evaluation of single ensemble members of all splits
}
