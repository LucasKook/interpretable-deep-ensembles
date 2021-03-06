% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnn.R
\name{cnn_stroke}
\alias{cnn_stroke}
\title{CNN for stroke data}
\usage{
cnn_stroke(
  output_shape = NULL,
  mbl = FALSE,
  ll_activation = "linear",
  dropout_rate = 0.3,
  ll_bias = FALSE,
  input_shape = dim(im)[2:5],
  ...
)
}
\description{
CNN for stroke data
}
