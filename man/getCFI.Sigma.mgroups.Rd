% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effectSizes.R
\name{getCFI.Sigma.mgroups}
\alias{getCFI.Sigma.mgroups}
\title{getCFI.Sigma.mgroups}
\usage{
getCFI.Sigma.mgroups(
  SigmaHat,
  S,
  muHat = NULL,
  mu = NULL,
  N,
  fittingFunction = "ML"
)
}
\arguments{
\item{SigmaHat}{a list of model implied covariance matrix}

\item{S}{a list of observed (or population) covariance matrix}

\item{muHat}{model implied mean vector}

\item{mu}{observed (or population) mean vector}

\item{N}{a list of group weights}

\item{fittingFunction}{whether to use \code{ML} or \code{WLS}}
}
\value{
Returns CFI
}
\description{
Computes CFI given the model-implied and the observed (or population) covariance matrix for multiple group models.
\code{CFI = (F_null - F_hyp) / F_null} applying multiple group sampling weights to \code{F_hyp} and \code{F_null}.
}
