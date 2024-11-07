# Test environments
  - Ubuntu 20.04 (`devtools::check(".", remote = TRUE)`)
    - R version 4.4.1 (2024-06-14 r86737) x86_64-pc-linux-gnu

  - win-builder (devel and release): http://win-builder.r-project.org/
    - R version 4.4.1 (2024-06-14 ucrt) 
    - R Under development (unstable) (2024-11-03 r87286 ucrt)

  - MacOS (`devtools::check_mac_release(".")`)
    - R version 4.4.0 (2024-04-24)
    - Apple clang version 14.0.0 (clang-1400.0.29.202)
    - GNU Fortran (GCC) 12.2.0
    - macOS Ventura 13.3.1


# R CMD check results
There were no ERRORs or NOTEs.

The MacOS check produced the following irrelevant WARNING about a dependency `Warning: package ‘Matrix’ was built under R version 4.4.1`.

# Downstream dependencies
I have checked and there are no current downstream dependencies of `gremlin`

