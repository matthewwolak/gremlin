# Resubmission
This package was recently archived on CRAN owing to my inability to make changes to the package in time. These changes were necessary after an update to a function in `base` R fundamentally changed the outcome of the function evaluation which caused my package to no longer work. I have now fixed the ERRORs that came about after the change to `base`.


# Test environments
  - Ubuntu 19.04
    - R 3.6.1 (2019-07-05) x86_64-pc-linux-gnu (64-bit)

  - win-builder (devel and release): http://win-builder.r-project.org/
    - R version 4.0.2 (2020-06-22), platform: x86_64-w64-mingw32 (64-bit) 
    - R Under development (unstable) (2020-06-23 r78741), platform: x86_64-w64-mingw32 (64-bit)

  - R-hub (`devtools::check_rhub(".", interactive = FALSE)`)
    - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
    - Ubuntu Linux 16.04 LTS, R-release, GCC
    - Fedora Linux, R-devel, clang, gfortran
    - Debian Linux, R-devel, GCC ASAN/UBSAN

# R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

Further, the NOTE indicates **possibly**  mis-spelled words in the DESCRIPTION. However, these words have been checked and are correctly spelled.

# Downstream dependencies
I have checked and there are no current downstream dependencies of `gremlin`


