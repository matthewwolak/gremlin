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
  - From [R](https://CRAN.R-project.org/):
    - see the package page for the latest release of [gremlin on CRAN](https://CRAN.R-project.org/package=gremlin) where you can download the source.
    - install the latest release of the package directly in R:
   ```R
   install.packages("gremlin")
   ```
   then select your favorite [CRAN mirror](https://CRAN.R-project.org/)

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

  - Estimating autosomal additive and dominance genetic variances
```
library(gremlin)
library(nadiv)  #<-- needed for creating inverse relatedness matrices

# Add unique term for including individual effects of additive and dominance
warcolak$IDD <- warcolak$ID

# Create generalized inverse matrices
Ainv <- makeAinv(warcolak[, 1:3])$Ainv
Dinv <- makeD(warcolak[, 1:3])$Dinv

# Basic model structure is as follows:
## Fixed effects of sex
## ID  = autosomal additive genetic variance term
## IDD = autosomal dominance genetic variance term
grAD <- gremlin(trait1 ~ sex-1,
	random = ~ ID + IDD,
	ginverse = list(ID = Ainv, IDD = Dinv),
	data = warcolak)

# Summary
nrow(warcolak)
summary(grAD)

# Calculate proportions of phenotypic variances (and Std. Error)
deltaSE(h2 ~ V1 / (V1 + V2 + V3), grAD)
deltaSE(d2 ~ V2 / (V1 + V2 + V3), grAD)

# Likelihood Ratio Test: Hypothesis test domimance variance=0
## Do this 2 alternative ways - both use `update()`:

### Either fix dominance variance to *almost* zero
grA_Dfxd <- update(grAD, Gstart = list(0.1, 1e-8), Gcon = list("P", "F"))

### Or drop dominance variance from the model
grA <- update(grAD, random = ~ ID)

## Compare log-likelihoods
logLik(grA_Dfxd)
logLik(grA)

## Do the Hypothesis test:
anova(grA, grAD)

```
