
gbdcd
=====

[![Travis Build Status](https://travis-ci.org/leandromineti/gbdcd.svg?branch=master)](https://travis-ci.org/leandromineti/gbdcd) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/7cv1gidywu5sxila?svg=true)](https://ci.appveyor.com/project/leandromineti/gbdcd) [![Coverage Status](https://codecov.io/gh/leandromineti/gbdcd/branch/master/graph/badge.svg)](https://codecov.io/gh/leandromineti/gbdcd) [![CRAN Status Badge](http://www.r-pkg.org/badges/version/gbdcd)](https://cran.r-project.org/package=gbdcd)

Overview
--------

An R package implementing the Bayesian Detection of Clusters and Discontinuities.

Installation
------------

``` r
library(devtools)

devtools::install_github("leandromineti/gbdcd")
```

Usage
-----

``` r
library(gbdcd)

data("aneeldata", package = "gbdcd")
data("aneelshape", package = "gbdcd")

target_variable <- aneelshape$z_Precipitation
neighbors <- aneeldata$connections

res <- gaussianBDCD(y = target_variable, 
                    viz = neighbors, 
                    c = 0.35, 
                    n_iterations = 10000, 
                    burn_in = 5000, 
                    mu0 = 0, 
                    sigma0 = sqrt(2))
```

### Todo list

-   \[x\] basic unit testing for all functions.
-   \[ \] turn all variable and function names into English.
-   \[ \] improve README with results.
-   \[ \] improve function documentation.
-   \[ \] publish package on CRAN.
