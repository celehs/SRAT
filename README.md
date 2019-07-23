__Sequence Robust Association Test (SRAT)__

## Overview

SRAT (Dai et al. 2017) <[doi:10.1111/biom.12643](https://doi.org/10.1111/biom.12643)> is a fully rank-based and flexible approach to test for association between a set of genetic variants and an outcome, while accounting for within-family correlation and adjusting for covariates. SRAT includes the well-known Wilcoxon rank sum test as a special case. Comparing to existing methods, SRAT has the advantages of allowing for unknown correlation structures and weaker assumptions about the outcome distribution. 

## Installation

If `devtools` is not installed, uncomment the code below and install it from CRAN.

``` r
# install.packages("devtools")
```

Install development version from GitHub:

``` r
devtools::install_github("celehs/SRAT")
```

Load package `SRAT` into R:

``` r
library(SRAT)
```

## References

Dai, W., Yang, M., Wang, C., & Cai, T. (2017). Sequence robust association test for familial data. Biometrics, 73(3), 876–884. <[doi:10.1111/biom.12643](https://doi.org/10.1111/biom.12643)>
