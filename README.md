# gremlin <!-- [![](http://www.r-pkg.org/badges/version/nadiv)](https://cran.r-project.org/package=nadiv) [![](http://cranlogs.r-pkg.org/badges/grand-total/nadiv)](http://cranlogs.r-pkg.org/badges/grand-total/nadiv) -->

[`R`](https://cran.r-project.org/) package for mixed-effects model **REML** incorporating **G**eneralized **In**verses (so, with some mental gymnastics: **GREMLIN**).


## See the latest developments:
 - gremlin [NEWS page](https://github.com/matthewwolak/gremlin/blob/master/NEWS.md)


## Overview of main branches:
  - `master` branch is the most recent production version (often the same as what is available from the [R CRAN mirrors](https://cran.r-project.org/))
  - `devel` branch is a preview of the next release which *should* be functional and error/bug free, but proceed with caution
  - `gremlinR` branch is a **pure R** implementation. This serves as a testing ground and as a learning tool to see what `gremlin` is doing 'under the hood' using just the R language. This may lag behind `devel` and `master`, but see the [NEWS page](https://github.com/matthewwolak/gremlin/blob/master/NEWS.md) for an overview of the most recent changes!

## To obtain gremlin:
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
   - clone or download the latest development version here
   - install the latest development version directly in R using the `devtools` package [https://github.com/hadley/devtools](https://github.com/hadley/devtools):
   ```R
   library(devtools); install_github("matthewwolak/gremlin")
   ```


