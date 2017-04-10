\name{summary.gremlin}
  \alias{summary.gremlin}
  \alias{print.summary.gremlin}

\title{
Gremlin model summary
}

\description{
Summarize and print results of linear mixed model fitted with gremlin.
}

\usage{
  \method{summary}{gremlin}(object, \dots)
  \method{print}{summary.gremlin}(x, digits, \dots)
}

\arguments{
  \item{object, x }{An object of \code{class} \sQuote{gremlin} or \sQuote{summary.gremlin}.}
  \item{digits }{An \code{integer} used for number formatting with \sQuote{signif()}.} 
  \item{\dots }{Additional arguments to be passed to control the output.}
}


\value{A \code{list} of class \code{summary.gremlin} or a printed value to the screen with no return values.

Returned values include the log-likelihood, model formulae, table of variance components and approximate standard errors (calculated from the inverse of the average information matrix), and fixed effects and standard errors (calculated from the corresponding diagonal elements of the inverse of the coefficient matrix).

}


\author{\email{matthewwolak@gmail.com}}

\seealso{
  \code{\link{gremlin}}
}

