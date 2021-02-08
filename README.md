
# regchange

The goal of `regchange` is to fit a linear regression model with
change-points in the regression coefficients by a regularized least
squares criterion. The algorithm computes the exact solution by using a
dynamic programming approach or a local optimum by a binary segmentation
method.

## Installation

You can install this version of `regchange` with:

``` r
library(devtools)
install_github("florencialeonardi/regchange")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(regchange)
x <- matrix(runif(150), ncol=3, nrow=50)
y <- x%*%c(1,0,0)+rnorm(50, 0.01)
fit <- regchange(x,y)
#> Warning: from glmnet Fortran code (error code -75); Convergence for 75th lambda
#> value not reached after maxit=100000 iterations; solutions for larger lambdas
#> returned
```
