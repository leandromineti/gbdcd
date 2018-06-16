
gbdcd
=====

[![Travis Build Status](https://travis-ci.org/leandromineti/gbdcd.svg?branch=master)](https://travis-ci.org/leandromineti/gbdcd) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/7cv1gidywu5sxila?svg=true)](https://ci.appveyor.com/project/leandromineti/gbdcd) [![Coverage Status](https://codecov.io/gh/leandromineti/gbdcd/branch/master/graph/badge.svg)](https://codecov.io/gh/leandromineti/gbdcd) [![DOI](https://zenodo.org/badge/102520939.svg)](https://zenodo.org/badge/latestdoi/102520939)

Overview
--------

An R package implementing the Bayesian Detection of Clusters and Discontinuities.

Installation
------------

``` r
library(devtools)

devtools::install_github("leandromineti/gbdcd")
```

If you are a Windows user, make sure you have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed.

Usage
-----

``` r
library(gbdcd)

data("aneeldata", package = "gbdcd")
data("aneelshape", package = "gbdcd")

target_variable <- aneelshape$z_Precipitation
neighbors <- aneeldata$connections

out <- gaussianBDCD(y = target_variable, 
                    neigh = neighbors, 
                    c = 0.35, 
                    n_iterations = 100000, 
                    burn_in = 50000, 
                    mu0 = 0, 
                    sigma0 = sqrt(2))
```

Results
-------

From the results obtained in the previous function we can look at many variables like: posterior distribution, proximity dendogram and the map regionalization.

### Posterior distribution

``` r
barplot(table(out$k.MCMC)/sum(table(out$k.MCMC)), col="light blue",
        xlab="number of groups", ylab="estimated probability")
grid()
barplot(table(out$k.MCMC)/sum(table(out$k.MCMC)), add=T, col="light blue")
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

### Proximity dendogram

``` r
n_cluster <- as.numeric(names(which.max(table(out$k.MCMC))))

plot(out$cluster.info, ylab="d(i~j)", main="Proximity Dendrogram", xlab="", 
     labels=aneelshape$DMU, cex=0.6)

groups <- cutree(out$cluster.info, k = n_cluster)
rect.hclust(out$cluster.info, k=n_cluster, border="red")
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

### Map regionalization

``` r
library(RColorBrewer)
library(spdep)

colrs <- brewer.pal(n_cluster, "Set3")

par(mfrow = c(1,2))
plot(aneelshape, col=colrs[groups])
boxplot(aneelshape$z_Precipitation ~ groups, col=colrs,
        ylab="Precipitation", xlab="Geographical clusters")
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

### Todo list

-   \[x\] basic unit testing for all functions.
-   \[x\] turn all variable and function names into English.
-   \[x\] improve README with results.
-   \[x\] improve function documentation.
-   \[ \] publish package on CRAN.
