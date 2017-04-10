\name{logLik.gremlin}
\alias{logLik.gremlin}

\title{
log-likelihood of a gremlin model
}
\description{
Extracts the log-likelihood from a gremlin model fit
}

\usage{
  \method{logLik}{gremlin}(object, \dots)
}

\arguments{
  \item{object }{An object of \code{class} \sQuote{gremlin}.}
  \item{\dots }{Additional arguments to be passed to control the model fitting.}
}

\value{A \code{numeric} value for the log-likelihood.}

\author{\email{matthewwolak@gmail.com}}

\examples{
    mod11 <- gremlinR(WWG11 ~ sex - 1,
	random = ~ calf,
	data = Mrode11,
	Gstart = matrix(0.1), Rstart = matrix(0.4),
	maxit = 50, v = 2, algit = "EM")
    logLik(mod11)
}
