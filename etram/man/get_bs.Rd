% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{get_bs}
\alias{get_bs}
\title{Bootstrap confidence intervals across all splits}
\usage{
get_bs(
  lys_cdf_all,
  y_true_all,
  met_ref = NULL,
  R = 1000,
  ncpus = 10,
  weights = rep(list(rep(1, length(lys_cdf_all[[1]]))), length(lys_cdf_all)),
  binary = FALSE,
  stratified = FALSE
)
}
\arguments{
\item{lys_cdf_all}{list of lists. Each sublist contains the CDFs of all ensemble members
or the ensemble CDF per split.}

\item{y_true_all}{list of all observed responses (one-hot encoded).}

\item{met_ref}{optional. Test performance of the reference model (e.g. simple linear shift model)
as a data frame. Must contain the variables \code{spl} and \code{metric} with the following levels:
\code{c('nll', 'brier', 'eauc', 'eacc', 'cint', 'cslope')} or
\code{c('nll', 'rps', 'eqwk', 'eacc', 'cint', 'cslope')}.}

\item{R}{number of bootstrap samples.}

\item{ncpus}{Number of processes to be used in parallel operation.}

\item{weights}{list of optimized weights. Each element contains the weights for each
ensemble member as numeric vector.}

\item{binary}{logical. Whether binary or ordinal metrics should be computed
(binary metrics: NLL, Brier score, 1-AUC, 1-ACC, CITL, C Slope,
ordinal metrics: NLL, RPS, 1-QWK, 1-ACC, CITL, C Slope)}

\item{stratified}{logical. Whether samples should be taken stratified by split and member.}
}
\description{
Calculates bootstrap median and confidence intervals across all splits
for binary (NLL, Brier score, 1-AUC, 1-accuracy, calibration-in-the-large, calibration slope)
or ordinal (NLL, RPS, 1-QWK, 1-accuracy, calibration-in-the-large, calibration slope) metrics.
If desired the confidence interval is additionally calculated for the test error relative to a
reference model (performance treated as fixed).
}
