__Sequence Robust Association Test (SRAT)__

## Overview

SRAT (Dai et al. 2017) is a fully rank-based and flexible approach to test for association between a set of genetic variants and an outcome, while accounting for within-family correlation and adjusting for covariates. SRAT includes the well-known Wilcoxon rank sum test as a special case. Comparing to existing methods, SRAT has the advantages of allowing for unknown correlation structures and weaker assumptions about the outcome distribution. 

## Installation

If `devtools` is not installed, uncomment the code below and install it from CRAN.

``` r
# install.packages("devtools")
```

Install development version from GitHub:

``` r
devtools::install_github("celehs/SRAT")
```

## References

W Dai, M Yang, C Wang, T Cai. Sequence robust association test for familial data. _Biometrics_, 2017, 73(3):876–884. <https://doi.org/10.1111/biom.12643>
