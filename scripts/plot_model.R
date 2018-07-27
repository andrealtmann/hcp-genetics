
vector2matrix <- function(x){
  N <- length(x)
  xdim <- sqrt(2 * N + 0.25) + 0.5

  result <- matrix(NA, nrow=xdim, ncol=xdim)
  #diag(result) <- 0

  i <- 1
  j <- i+1
  for(k in 1:N){
    result[i,j] <- x[k]
    result[j,i] <- x[k]
    j <- j+1
    if (j > xdim){
      i <- i + 1
      j <- i + 1
    }
  }
  return(result)
}


glmnet2matrix <- function(x, mod="min", lidx=NA){

  if (is.na(lidx)){
    if (mod == "min"){
      slambda <- x$lambda.min
    }
    if (mod == "1se"){
      slambda <- x$lambda.1se
    }
    lidx <- which(x$lambda == slambda)
  }
  return( vector2matrix(x$glmnet.fit$beta[,lidx]))

}

regionCnt <- function(resMat){

  cs <-colSums(resMat!=0, na.rm=T)
  plot(cs, type="h")

}
