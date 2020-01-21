## ------------------------------------------------------------------------
library(SRAT)
library(SKAT)

## ------------------------------------------------------------------------
data(SKAT.example)
attach(SKAT.example)

## ------------------------------------------------------------------------
obj.skat <- SKAT_Null_Model(y.c ~ X, out_type = "C")
SKAT(Z, obj.skat, weights.beta = c(1, 1))$p.value

## ------------------------------------------------------------------------
set.seed(123)
obj.srat <- srat.null(y.c, X)
srat.test(Z, obj.srat)

