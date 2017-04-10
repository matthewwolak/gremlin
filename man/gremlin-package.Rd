\name{gremlin-package}
\alias{gremlin-package}
\docType{package}
\title{
\packageTitle{gremlin}
}
\description{
\packageDescription{gremlin}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{gremlin}
\packageIndices{gremlin}

\code{gremlin} is a package for fitting linear mixed-effects models to Gaussian data using Restricted Maximum Likelihood (REML) for Inference. Special variance structures can be fitted through the use of generalized inverse covariance or correlation matrices for random effects (e.g., quantitative genetic models commonly known as 'animal models').

The package also can implement the average information algorithm to efficiently maximizing the log-likelihood (Thompson & Johnson 1995; Gilmour et al. 1995; Meyer & Smith 1996). The average information algorithm combined with sparse matrix techniques can potentially make model fitting very efficient.

}

\author{
\packageAuthor{gremlin}

Maintainer: \packageMaintainer{gremlin}
}

\references{
Mrode. 2005.
Meyer & Smith. 1996.
Gilmour et al. 1995.
Thompson & Johnson. 1995.
}

\keyword{ package }
\seealso{
  \code{\link[MCMCglmm:MCMCglmm-package]{MCMCglmm}}
}
\examples{
#TODO: simple examples of the most important functions
  \dontrun{
    library(nadiv)
    Ainv <- makeAinv(Mrode3[-c(1:2), 1:3])$Ainv
    mod11 <- gremlinR(WWG11 ~ sex - 1,
	random = ~ calf,
	data = Mrode11,
	ginverse = list(calf = Ainv),
	Gstart = matrix(0.2), Rstart = matrix(0.4),
	maxit = 10, v = 2, algit = "AI")
  }

}
