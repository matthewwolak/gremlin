# 1.0.0.1
## NEW
  - `update()` function
    - can now continue a model where it left off or change the structure (e.g., drop a single variance component for likelihood ratio test)

## Mnior Changes
  - Implement more efficient algorithms in the c++ code, that were developed in the R code for version 1.0.0.0.
  - Add `gremlinControl()` function for _advanced_ changes to the way gremlin runs

# 1.0.0.0
## NEW
  - Completely revised way models are built and called
    - made a "modular" series of functions for setting up the model and optimizing the REML likelihood
    - __new__ `grMod` and `gremlinR` classes.
        - `grMod` is the model structure for which a log-likelihood can be calculated
        - `gremlinR` class distinguishes from `gremlin` class in that `gremlinR` objects will only use `R` code written by the package in order to run the model. Class `gremlin` will execute underlying c++ code written in the package.

  - __Average Information__ algorithm has been _vastly_ improved
    - `ai()` efficiently calculates the AI matrix without directly computing several matrix inverses (as previously coded)

  - `lambda` and alternative parameterizations now possible and executed by the same code
    - `lambda` parameterization is the REML likelihood of the variance ratios after _factoring out a residual variance from the Mixed Model Equations_.
    - the _alternative_ does not have a special name, this is just a model of all (co)variance parameters as (co)variance parameters (as opposed to ratios, as in the `lambda` models).
    - instead of completely separate functions for these two parameterizations, there is an argument that runs alternative lines of code, wherever the calculations differ for these two different parameterizations

## Minor Changes
  - No long construct Mixed Model Array (`M`) matrix from which the Cholesky factorization (and `logDetC` and `tyPy` calculations are made)
    - Changed to directly construct coefficient matrix of mixed model equations (`C`) and obtain `tyPy` and `logDetC` using this
    - Previously had to store Cholesky factorizations of both `M` and `C`, now do a `solve` with Cholesky of `C` (`sLc`/`Lc` in `R`/`c++` code) to calculate `tyPy` based off Boldman and Van Vleck


# 0.1.0.0
## NEW
  - methods for `gremlin` objects
    - notably, `AIC`, `residuals`, `anova`, and `nobs`
    - updated the `summary`, `print`, and `logLik` methods as well

# 0.0.2.0
Improved algorithm that reduces computational resources and time! Also implemented c++ code in `gremlin()`, while keeping `gremlinR()` purely the R implementation (at least from the package writing standpoint).

# 0.0.1.0
## NEW
Documentation has switched from filling out the `.Rd` files manually to providing
documentation next to the function code in the `.R` files using `roxygen2`


# 0.0.0.1 April 2017 `gremlin` is born!

Congratulations, its a gremlin!

