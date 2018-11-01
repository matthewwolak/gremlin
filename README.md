# gremlin
[![](http://www.r-pkg.org/badges/version/gremlin)](https://cran.r-project.org/package=gremlin)
[![](http://cranlogs.r-pkg.org/badges/grand-total/gremlin)](http://cranlogs.r-pkg.org/badges/grand-total/gremlin)
[![DOI](https://zenodo.org/badge/87194564.svg)](https://zenodo.org/badge/latestdoi/87194564)





[`R`](https://cran.r-project.org/) package for mixed-effects model **REML** incorporating **G**eneralized **In**verses (so, with some mental gymnastics: **GREMLIN**).


## Cite as:
If your use of the `gremlin` package contributes to a publication, please cite the package as (and feel free to let me know!):

>Wolak, M.E. 2018. gremlin: R package for mixed-effects model REML incorporating generalized Inverses. Version *Insert version here*. Zenodo [http://doi.org/10.5281/zenodo.1476565](http://doi.org/10.5281/zenodo.1476565).


## See the latest developments:
 - gremlin [NEWS page](https://github.com/matthewwolak/gremlin/blob/master/NEWS.md)


## Overview of main branches:
  - `master` branch is the most recent production version (often the same as what is available from the [R CRAN mirrors](https://cran.r-project.org/))
  - `devel` branch is a preview of the next release which *should* be functional and error/bug free, but proceed with caution

## To install gremlin:
<!--
 - From [R](https://CRAN.R-project.org/):
   - see the package page for the latest release of [gremlin on CRAN](https://CRAN.R-project.org/package=gremlin) where you can download the source.
   - install the latest release of the package directly in R:
   ```R
   install.packages("gremlin")
   ```
   then select your favorite [CRAN mirror](https://CRAN.R-project.org/)
-->   
 - From GitHub:
   - install the latest versions directly in R using the `devtools` package [https://github.com/hadley/devtools](https://github.com/hadley/devtools):
   ```
   library(devtools)

   # Install `master` branch
   install_github("matthewwolak/gremlin")
   
   # Install `devel` branch
   install_github("matthewwolak/gremlin", ref = "devel")
   ```

## Examples

### `gremlinR`
  - Estimating autosomal additive and dominance genetic variances
```
library(gremlin)
library(nadiv)  #<-- needed for creating inverse relatedness matrices

# Set up a subset of data for the example
tmpw <- warcolak[1:(72*10), ]  #<-- only use 10 units for poor young gremlin
tmpw$ID <- as.factor(as.character(tmpw$ID))
tmpw$IDD <- tmpw$ID

# Create generalized inverse matrices
Ainv <- makeAinv(tmpw[, 1:3])$Ainv
Dinv <- makeD(tmpw[, 1:3])$Dinv
  Dinv@Dimnames[[1]] <- as.character(tmpw$ID)

# Basic model structure is as follows:
## Fixed effects of sex
## ID  = autosomal additive genetic variance term
## IDD = autosomal dominance genetic variance term
system.time(grAD <- gremlinR(trait1 ~ sex-1,
	random = ~ ID + IDD,
	ginverse = list(ID = Ainv, IDD = Dinv),
	data = tmpw,
	maxit = 25, v = 2, vit = 1))

nrow(tmpw)
summary(grAD)
```
