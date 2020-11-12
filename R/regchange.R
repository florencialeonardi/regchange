#' Inhomogeneous linear regression model
#'
#' This function fits a linear regression model with change-points in the regression coefficients by a regularized
#' least squares criterion. The algorithm computes the exact solution by using a dynamic programming approach.
#'
#' @param x matrix of covariates of dimension nxp
#' @param y response vector of length n
#' @param kmax maximum number of segments in the model
#' @param gamma penalized factor for the number of segments in the model
#' @param lambda vector of penalized parameters for the Lasso estimator on each segment
#' @param delta minimum proportion of sample points on each segment
#' @param intercept should the model have an intercept
#' @param train vector with train points
#' @param type either "glmnet" or "lars"
#'
#' @return a list containing the change point vector, a matrix with the regression parameters
#' on each segment and a matrix with the residual sum of squares
#' \describe{
#'  \item{\code{alpha}}{an array of dimension kmax x  kmax x length(lambda). For each value of lambda and number of segments j, the vector
#'     alpha\[, 1:j, lambda\] contains the change-point proportions with respect to the sample size n. }
#'  \item{\code{beta}}{a list of dimension kmax x length(lambda) containing in each entry the matrix of regression coefficients}
#'  \item{\code{rss}}{matrix of dimension kmax x length(lambda) with the sum of squared residuals for each model}
#'  }
#'
#' @export
#'
#' @examples
#' x <- matrix(runif(300), ncol=3, nrow=100)
#' y <- runif(100)
#' fit <- regchange(x,y)
#'
#' @references
#' Leonardi, F  and Bühlmann, P.  Computationally efficient change point detection for high-dimensional regression.
#' \href{https://arxiv.org/abs/1601.03704}{arXiv:1601.03704}
#'
#'
regchange <- function(x,y,kmax=5,gamma=0.01,lambda=seq(1,0.01,-0.01),delta=0.0,intercept=TRUE,train=seq(1,length(y)),type="glmnet")
{
  if(length(train) == nrow(x)) n <- nrow(x)
  else n <- length(train)
  p <- ncol(x)
  xt <- x[train,,drop=FALSE]
  yt <- y[train]
  nmin <- max(as.integer(n*delta),2)
  if (kmax > n) {
    print("error: kmax can not be bigger than the sample size length(y)")
    return()
    }
  else segmax <- kmax
  ## compute rss array for the lambda sequence
  dev_fit <- array(Inf,dim=c(n,n,length(lambda)))
  for(i in 1:(n-1)){
    for (j in (i+1):n){
      if( (j-i+1) >= nmin ){
        if (type=="glmnet") {
          fit <- glmnet::glmnet(xt[i:j,,drop=FALSE],yt[i:j],intercept=intercept)
          dev_fit[i,j,] <- apply((yt[i:j]-glmnet::predict.glmnet(fit,xt[i:j,,drop=FALSE],s=lambda))^2/n,2,sum)
        }
        else if(type=="lars") {
          fit <- lars::lars(xt[i:j,,drop=FALSE],yt[i:j],intercept=intercept)
          dev_fit[i,j,] <- apply(as.matrix((yt[i:j]-lars::predict.lars(fit,xt[i:j,,drop=FALSE],s=lambda,mode="lambda")$fit)^2/n),2,sum)
        }
        else {
          print("No valid method selected")
          return()
        }
      }
    }
  }
  ## create output objects
  rss_mat <- c()
  alpha_mat <- array(dim=c(segmax,segmax,length(lambda)))
  beta_mat <-  vector(length=segmax*length(lambda))
  beta_mat <- as.list(beta_mat)
  dim(beta_mat) <- c(segmax,length(lambda))
  ## compute rss, alpha and beta for each value of lambda
  zaux <- matrix(nrow=segmax,ncol=n)
  raux <- matrix(nrow=segmax-1,ncol=n)
  for ( lm in 1:length(lambda) ){
    zaux[1,] <- dev_fit[1,,lm]
    if ( segmax > 1){
      for(k in 2:segmax){
        for(j in k:n){
          vec <- c()
          for(i in k:j){
            vec <- c(vec,zaux[k-1,i-1]+dev_fit[i,j,lm])
          }
          zaux[k,j] <- min(vec)
          raux[k-1,j] <- which.min(vec) + k - 2
        }
      }
    }
    rss_mat <- cbind(rss_mat, zaux[,n])
    for (k in segmax:1){
      alpha_mat[k,k,lm] <- n
      beta_aux <- c()
      if ( k > 1 ) {
        for(j in (k-1):1){
          alpha_mat[k,j,lm] <- raux[j,alpha_mat[k,j+1,lm]]
          segi <- alpha_mat[k,j,lm]+1
          segf <- alpha_mat[k,j+1,lm]
          if (type=="glmnet") {
            fit <- glmnet::glmnet(xt[segi:segf,,drop=FALSE],yt[segi:segf],intercept=intercept)
            beta_aux <- cbind(glmnet::coef.glmnet(fit,s=lambda[lm]),beta_aux)
          }
          else if(type=="lars") {
            fit <- lars::lars(xt[segi:segf,,drop=FALSE],yt[segi:segf],intercept=intercept)
            a0 <- lars::predict.lars(fit,t(rep(0,p)), s=lambda[lm],mode="lambda")$fit
            beta_aux <- cbind(c(a0,lars::coef.lars(fit,s=lambda[lm],mode="lambda")),beta_aux)
          }
        }
      }
      segi <- 1
      segf <- alpha_mat[k,1,lm]
      if (type=="glmnet") {
        fit <- glmnet::glmnet(xt[segi:segf,,drop=FALSE],yt[segi:segf],intercept=intercept)
        beta_aux <- cbind(glmnet::coef.glmnet(fit,s=lambda[lm]),beta_aux)
      }
      else {
        fit <- lars::lars(xt[segi:segf,,drop=FALSE],yt[segi:segf],intercept=intercept)
        a0 <- lars::predict.lars(fit,t(rep(0,p)), s=lambda[lm],mode="lambda")$fit
        beta_aux <- cbind(c(a0,lars::coef.lars(fit,s=lambda[lm],mode="lambda")),beta_aux)
      }
      beta_mat[[k,lm]] <- beta_aux
    }
  }
  #  rescale alpha if length(train)<nrow(x)
  if(n < nrow(x)){
    for(lm in 1:length(lambda)){
      alpha_mat[1,1,lm] <- nrow(x)
      for(k in 2:segmax){
        alpha0 <- c(train[alpha_mat[k,1:(k-1),lm]],nrow(x))
        alpha_mat[k,1:k,lm] <- alpha0
      }
    }
  }
  outlist <- list(alpha=alpha_mat/nrow(x), beta=beta_mat, rss=rss_mat)
  outlist
}
