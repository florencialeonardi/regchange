
# regchange

The goal of `regchange` is to fit a linear regression model with
change-points in the regression coefficients by a regularized least
squares criterion. The algorithm computes the exact solution by using a
dynamic programming approach.

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
x <- matrix(runif(300), ncol=3, nrow=100)
beta <- c(1,1,0)
y <- x%*%beta+rnorm(100,0,0.1)
fit <- regchange(x,y)
```
