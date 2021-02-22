
# SRAT: Sequence Robust Association Test

## Overview

Implements a fully rank-based and flexible approach to test for
association between a set of genetic variants and an outcome, while
accounting for within-family correlation and adjusting for covariates.
SRAT includes the well-known Wilcoxon rank sum test as a special case.
Comparing to existing methods, SRAT has the advantages of allowing for
unknown correlation structures and weaker assumptions about the outcome
distribution.

## Installation

If `remotes` is not installed, uncomment the code below and install it
from CRAN.

``` r
install.packages("remotes")
```

Install development version from GitHub:

``` r
remotes::install_github("celehs/SRAT")
```

## References

W Dai, M Yang, C Wang, T Cai. Sequence robust association test for
familial data. *Biometrics*, 2017, 73(3):876–884.
<https://doi.org/10.1111/biom.12643>
