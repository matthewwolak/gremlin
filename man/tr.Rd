\name{tr}
  \alias{tr}
  \alias{tr.default}
  \alias{tr.dgCMatrix}
  \alias{tr.dsCMatrix}

\title{
Matrix trace methods
}

\description{
Methods to efficiently calculate a matrix trace depending on the class of matrix
}

\usage{
  tr(X, \dots)

  \method{tr}{default}(X, \dots)
  \method{tr}{dgCMatrix}(X, \dots)
  \method{tr}{dsCMatrix}(X, \dots)
}

\arguments{
  \item{X }{A matrix.}
  \item{\dots }{Additional arguments.}
}

\value{A \code{numeric} value for the sum of the diagonal elements.}

\author{\email{matthewwolak@gmail.com}}



