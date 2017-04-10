\name{vech2matlist}
\alias{vech2matlist}

\title{
Vector to list of matrices}
\description{
Converts a vector of (co)variance parameters to a list of covariance matrices
}
\usage{
vech2matlist(vech, skeleton)
}
\arguments{
  \item{vech }{A \code{vector} of (co)variance parameters.}
  \item{skeleton }{An example structure to map \code{vech} onto.}
}

\value{
A list of matrices of the same structure as \code{skeleton}
}

\author{\email{matthewwolak@gmail.com}}

\examples{
#TODO
}

