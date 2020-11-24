utils::globalVariables(c("ite"))

lasso.like <- function(x,y,n,lambda,intercept)
{
  fit <- glmnet::glmnet(x,y,lambda=lambda,intercept=intercept)
  mat <- glmnet::predict.glmnet(fit,newx=x,s=lambda)
  lasso.like <- sum((y-mat)^2)/n
}

hieralg.tree <- function(x,y,n,lambda,findex,tree,threshold,max.depth)
{
  m <- nrow(x)
  print(paste0(c("[",findex," , ",m+findex-1,"]"),sep="",collapse=""))
  if (m > 3){
    z <- foreach::foreach(ite=2:(m-2), .combine=c, .export=c("lasso.like"), .packages=c("glmnet") ) %dopar% {
      lasso.like(x[1:ite,,drop=FALSE],y[1:ite],n,lambda)+lasso.like(x[(ite+1):m,,drop=FALSE],y[(ite+1):m],n,lambda)
    }
    k <- which.min(z)  # i = k + 1
    gamma2 <- z[k]-tree$loglike
    if ( -gamma2 > threshold && tree$depth < max.depth ){
      tree$gamma <- -gamma2
      val1 <- lasso.like(x[1:(k+1),,drop=FALSE],y[1:(k+1)],n,lambda)
      val2 <- lasso.like(x[(k+2):m,,drop=FALSE],y[(k+2):m],n,lambda)
      child.left <- tree$AddChild(paste0(c("[",findex,",",k+findex,"]"),collapse=""),gamma=0,loglike=val1,seg.left=findex,seg.right=k+findex,depth=tree$depth+1)
      child.right <- tree$AddChild(paste0(c("[",k+1+findex,",",m+findex-1,"]"),collapse=""),gamma=0,loglike=val2,seg.left=k+1+findex,seg.right=m+findex-1,depth=tree$depth+1)
      hieralg.tree(x[1:(k+1),,drop=FALSE],y[1:(k+1)],n,lambda,findex,child.left,threshold,max.depth)
      hieralg.tree(x[(k+2):m,,drop=FALSE],y[(k+2):m],n,lambda,findex+k+1,child.right,threshold,max.depth)
    }
  }
  tree
}

binalg <- function(x,y,lambda,gamma,max.depth,intercept=TRUE)
{
  ncores <- parallel::detectCores()
  cl <- parallel::makeCluster(ncores, outfile="")
  doParallel::registerDoParallel(cl)
  n <- length(y)
  tree=data.tree::Node$new("binalg.root")
  val <- lasso.like(x,y,n,lambda)
  tree$loglike <- val
  tree$seg.left <- 1
  tree$seg.right <- n
  tree$depth <- 0
  hieralg.tree(x,y,n,lambda,1,tree,gamma,max.depth)
  df.tree <- data.tree::ToDataFrameTree(tree,"loglike","gamma","seg.left","seg.right","depth")
  parallel::stopCluster(cl)
}



