# 1.1.0
## NEW
  - added c++ (and R) code to use the Takahashi et al. algorithm for obtaining the partial inverse of the coefficient matrix
    - This switch means almost a 20-fold speedup (for c++ routines) when using analytical first derivatives (i.e., need to calculate "trace" terms)
    - Speedup comes from only calculating elements of inverse matrix (__C-inverse__) following the non-zero pattern determined for the Cholesky decomposition of the __C__ matrix.
    - follows the SuiteSparse Matlab_Tools `sparseinv` by Tim Davis, but sparseinv works on __LDL'__ factorization of __C__ whereas I changed this to work on __LL'__ factorization of __C__.
        - Takahashi, Fagan, & Chin. 1973. Formation of a sparse bus impedance matrix and its application to short circuit study. 8th PICA Conference Proceedings, Minneapolis, MN.
    - `gremlinR()` now uses far less RAM per iteration of the model (previously was forming the entire __C-inverse__)

  - add finite difference algorithm to obtain first derivatives of likelihood function
    - introduced a parameter (`h`) inside `gremlinControl()` to set the "difference" or amount to alter parameters to calculate change in log-likelihood.
    
  - created a REML function inside the c++ code (`reml`) to calculate log-likelihood
    - moved log-likelihood calculation out of main program and reduced amount of repeated code.
    - also facilitated finite difference functions (which are just repeated log-likelihood evaluations)
    

## Minor Changes
  - Removed `error()` in c++
    - now issues with matrix singularities etc. do not stop code without returning model so far
    - should now be possible to use `update()` to get "through" trouble spots
    - also allows for user interruptions to c++ code from terminal
    
  - Changed default parameterization so _lambda_ transformation is __no longer the default__
  
  - Changed default convergence check criteria (`cctol`)
    - Models using previous values tended to only improve precision of estimates well beyond what was meaningful.
        

# 1.0.1 Released to CRAN 2020 June 25

## NEW
  - `deltaSE()` function to calculate approximate standard errors for functions of (co)variance parameters (e.g., h<sup>2</sup>, standard deviations of variances, or correlations)
    - this can take a formula for the function or a character expression
    - also allows for a list of formulas or character expressions
        e.g., calculate all variance components as proportions of total variance

  - Introduce `Gcon` and `Rcon` arguments to `gremlin()` for constraining parameters
    - enables parameters to be fixed or otherwise constrained
    - works in conjunction with the `Gstart` and `Rstart` arguments
    - For example in a simple `sire` model, we could restrain the `sire` variance `=0.38`.
```
grSf <- gremlin(WWG11 ~ sex,
	random= ~ sire,
	data = Mrode11,
	Gstart = list(matrix(0.38)),
	Gcon = list("F"),
	control = gremlinControl(lambda = FALSE))

```

  - Similar to above change (`Gcon`/`Rcon`), introduced steps to deal with parameters outside of the boundaries of their parameter space (e.g., variance < 0).
    - restrain these parameters to near their boundaries (after trying step-reduction calculation)
    - re-calculate Average Information, conditional on restrained parameters
        - See Gilmour. 2019. J. Anim. Breed. Genet. for specifics

  - change version numbering to just 3 numbers (instead of 4)
    - just dropping last number

## Minor Changes
  - create new c++ function to handle quasi Newton-Rhapson algorithm
    - allows secondary checks of appropriateness/naughtiness for proposed parameters based on a conditional AI algorithm (conditional on parameters restrained to boundary condition)


# 1.0.0.1
## NEW
  - `update()` function
    - can now continue a model where it left off or change the structure (e.g., drop a single variance component for likelihood ratio test)

  - Implement "step-halving" algorithm for AI updates
    - restricts parameter updates if AI algorithm proposes a change of >80% of original parameter value
    - amount by which a parameter change is restricted can be set in `gremlinControl()` using the `step` argument


## Minor Changes
  - Implement more efficient algorithms in the c++ code, that were developed in the R code for version 1.0.0.0.
  - Add `gremlinControl()` function for _advanced_ changes to the way gremlin runs
  - Begin major improvements to speed of gradient calculation function
    - changes to be incorporated in `em`, `ai`, and elsewhere (where relevant) in next version
    - implements calculations that take advantage of sparsity (i.e., don't calculate values where there are zeroes)


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


# 0.1.0.0 Released to CRAN 2018 October 30
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

