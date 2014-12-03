nordpred
========

Fit power5 and poisson Age-Period-Cohort models for prediction of cancer incidence

# Note
This package is work in progress!
Until otherwise stated, you should prefer the original source files from:
http://www.kreftregisteret.no/software/nordpred/


# Getting started

Run the following from R
```{r}
# Install devtools if not already done
install.packages("devtools")

# Install nordpred
devtools::install_github("cancercentrum/nordpred")

# List all exported functions (and data sets) from the package
library(nordpred)
help(package = "nordpred")

# Read the package vignette 
vignette("nordpred")
```
