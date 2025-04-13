# SRAT: Sequence Robust Association Test

## Overview

SRAT implements a fully rank-based and flexible approach to test for association between a set of genetic variants and an outcome, while accounting for within-family correlation and adjusting for covariates. The method is designed for genome-wide association studies (GWAS) and next-generation sequencing studies (NGSS) performed on family data.

### Key Features

- Fully rank-based, robust statistical method
- Accounts for within-family correlation in genetic studies
- Allows for unknown correlation structures
- Requires weaker assumptions about the outcome distribution
- Includes the well-known Wilcoxon rank sum test as a special case
- Provides better protection against type I error rate inflation
- More powerful for settings with skewed outcome distributions compared to existing methods

## Installation

Install development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("celehs/SRAT")
```

## Functions

- `srat()`: Main function for the Sequence Robust Association Test
- `srat.null()`: Fits the null model for SRAT
- `srat.test()`: Tests for genetic association using the fitted null model
- `obj.min()`: Helper function for objective function minimization

## Citation

Dai W, Yang M, Wang C, Cai T. Sequence robust association test for familial data. Biometrics. 2017 Sep;73(3):876-884. doi: 10.1111/biom.12643. Epub 2017 Mar 8. PMID: 28273695; PMCID: PMC11141465.
