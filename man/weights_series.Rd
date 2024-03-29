% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{weights_series}
\alias{weights_series}
\title{Generate weights series}
\usage{
weights_series(t_max, serial_parameters)
}
\arguments{
\item{t_max}{maximum time to generate weights until}

\item{serial_parameters}{named list of 'mean' and 'sd' of gamma distribution
characterising the serial interval distribution}
}
\value{
a vector of normalised weights of length t_max
}
\description{
Generate weights series
}
