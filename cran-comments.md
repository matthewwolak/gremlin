# Resubmission
This package was recently archived on CRAN owing to my inability to make changes to the package in time. These changes were necessary after an update to a function in `base` R fundamentally changed the outcome of the function evaluation which caused my package to no longer work. I have now fixed the ERRORs that came about after the change to `base`.


# Test environments
  -  Ubuntu 18.10
    - R 3.5.3 (2019-03-11), platform: x86_64-pc-linux-gnu (64-bit)
  -  win-builder (devel and release): http://win-builder.r-project.org/
    - R version 3.5.3 (2019-03-11), platform: x86_64-w64-mingw32 (64-bit) 
    - R version 3.6.0 alpha (2019-04-05 r76327), platform: x86_64-w64-mingw32 (64-bit)


# R CMD check results
There were no ERRORs or WARNINGs.

  - There was 1 NOTE:
    - Maintainer: 'Matthew Wolak <matthewwolak@gmail.com>'
    - New submission
    - Package was archived on CRAN

This simply notes (correctly) information about this **existing** package.

Further, the NOTE indicates **possibly**  mis-spelled words in the DESCRIPTION. However, these words have been checked and are correctly spelled.

# Downstream dependencies
I have checked and there are currently no downstream dependencies of `gremlin`
