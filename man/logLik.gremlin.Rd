% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.gremlin.R
\name{logLik.gremlin}
\alias{logLik.gremlin}
\alias{npar.gremlin}
\alias{AIC.gremlin}
\title{Methods to extract log-likelihood and information criterion of a gremlin
model.}
\usage{
\method{logLik}{gremlin}(object, ...)

npar.gremlin(object)

\method{AIC}{gremlin}(object, ..., k = 2, fxdDf = FALSE)
}
\arguments{
\item{object}{An object of \code{class} \sQuote{gremlin}.}

\item{\dots}{Additional arguments.}

\item{k}{A numeric value for the penalty per parameter. Default is 2, as in
classic AIC.}

\item{fxdDf}{A logical indicating whether to penalize according to the number
of fixed effect parameters. Since only models fit by REML can be compared,
these must always be the same and so become a constant. Hence, the default
is \code{FALSE}.}
}
\value{
\code{numeric} values for the log-likelihood, the number of
  parameters estimated by the model (sum of fixed effects and random effect
  (co)variance components), and Akaike's Information Criterion.
}
\description{
Extracts the log-likelihood or AIC from a gremlin model fit.
}
\details{
Function \code{npar.gremlin} returns an object with attributes \code{n.fxd}
and \code{n.bndry} which give additional information about the parameters
estimated and contributing to the overall \code{df} of the model. \code{n.fxd}
returns the total number of parameters (No. fixed effects + No. (co)variance
components) minus the number of parameters constrained to a certain value. Thus,
\code{n.fxd} represents the number of parameters that can vary and, as a 
consequence, affect the log-likelihood.

The attribute \code{n.bndry} reports the number of parameters that were
restrained to stay inside the boundaries of allowable parameter space (e.g.,
a variance that was not allowed to be negative).
}
\examples{
grS <- gremlin(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
logLik(grS)
AIC(grS)
}
\author{
\email{matthewwolak@gmail.com}
}
