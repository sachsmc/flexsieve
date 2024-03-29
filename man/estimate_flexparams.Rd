% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter-estimation.R
\name{estimate_flexparams}
\alias{estimate_flexparams}
\title{Estimate parameters in the flexsurv model}
\usage{
estimate_flexparams(datain, ID.var, form, K, ...)
}
\arguments{
\item{datain}{List of data frames as returned by \link{augment_data}}

\item{ID.var}{characted indicating id variable in data}

\item{form}{Formula specifying the model. The outcome should be a \link{Surv} object}

\item{K}{the number of knots}

\item{...}{Other arguments passed to \link{optim}}
}
\description{
Estimate parameters in the flexsurv model
}
