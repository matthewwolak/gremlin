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

