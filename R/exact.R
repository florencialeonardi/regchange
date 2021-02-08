
dynamic_segmentation <- function(xt,yt,segmax,gamma,lambda,delta,intercept,subs,type)
{
  ## compute rss array for the lambda sequence
  nmin <- max(as.integer(length(yt)*delta),2)
  p <- ncol(xt)
  n <- nrow(xt)
  m <- length(subs)
  dev_fit <- array(Inf,dim=c(m,m,length(lambda)))
  for(i in 1:(m-1)){
    for (j in (i+1):m){
      if( (subs[j]-subs[i]+1) >= nmin ){
        if (type=="glmnet") {
          fit <- glmnet::glmnet(xt[subs[i]:subs[j],,drop=FALSE],yt[subs[i]:subs[j]],intercept=intercept)
          dev_fit[i,j,] <- apply((yt[subs[i]:subs[j]]-glmnet::predict.glmnet(fit,xt[subs[i]:subs[j],,drop=FALSE],s=lambda))^2/n,2,sum)
        }
        else if(type=="lars") {
          fit <- lars::lars(xt[subs[i]:subs[j],,drop=FALSE],yt[subs[i]:subs[j]],intercept=intercept)
          dev_fit[i,j,] <- apply(as.matrix((yt[subs[i]:subs[j]]-lars::predict.lars(fit,xt[subs[i]:subs[j],,drop=FALSE],s=lambda,mode="lambda")$fit)^2/n),2,sum)
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
  zaux <- matrix(nrow=segmax,ncol=m)
  raux <- matrix(nrow=segmax-1,ncol=m)
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
    rss_mat <- cbind(rss_mat, zaux[,m])
    for (k in segmax:1){
      alpha_mat[k,k,lm] <- subs[m]
      beta_aux <- c()
      if ( k > 1 ) {
        for(j in (k-1):1){
          alpha_mat[k,j,lm] <- subs[raux[j,alpha_mat[k,j+1,lm]]]
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
  kest <- apply(rss_mat + (1:segmax)*gamma, 2, which.min)
  outlist <- list(alpha=alpha_mat, beta=beta_mat, rss=rss_mat, kest=kest)
  outlist
}



