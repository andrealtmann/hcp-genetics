
getCVfold <- function(x, nfold){
    avail <- 1:nrow(x)
    nsize <- ceiling(nrow(x)/nfold)
    cv.idx <- rep(nfold, nrow(x))
    for (i in 1:(nfold-1)){
        tmp <- sample(avail, nsize)
        cv.idx[tmp] <- i
        avail <- setdiff(avail, tmp)
    }
    return(cv.idx)
}

#outer loop for nested CV
doubleCV <- function(y, x, nfold, fam="binomial", measure="auc", alpha=0.5, lop="1se", nla = 100, myfold=c()){

    if (length(myfold) == 0){
      cv.idx <- getCVfold(x, nfold)
    } else {
        cv.idx <- myfold
        nfold  <- max(myfold)
    }

    my2CV.predict <- list()
    my2CV.labels  <- list()
    lambdas <- c()
    df      <- c()
    #run CV
    for(i in 1:nfold){
        message("outer fold: ",i)
        test.y <- y[cv.idx==i]
        test.x <- x[cv.idx==i,]
        train.y <- y[cv.idx!=i]
        train.x <- x[cv.idx!=i,]

        #inner cv.loop
        icv <- cv.glmnet(as.matrix(train.x), train.y, family=fam, nfold=5, nlambda=nla, type.measure=measure, alpha=alpha)
        #tmp <- glmnet(as.matrix(train.x), train.y, family=fam, lambda=icv$lambda, alpha=alpha)
        if (lop == "1se") {
          laop <- icv$lambda.1se
        } else {
          laop <- icv$lambda.min
        }
        iopt <- which(icv$lambda == laop)

        my2CV.predict[[i]] <- predict(icv$glmnet.fit, newx=as.matrix(test.x), s=laop, type="response")
        #my2CV.predict[[i]] <- predict(tmp, newx=as.matrix(test.x), type="response")[,iopt]
        my2CV.labels[[i]] <- test.y

        #print(cor(my2CV.predict[[i]], my2CV.labels[[i]]))
        #print(sqrt(mean((my2CV.predict[[i]] - my2CV.labels[[i]])^2)))

        lambdas <- c(lambdas, icv$lambda[iopt])
        df      <- c(df, icv$glmnet.fit$df[iopt])
    }
    res <- rep(0, nfold)
    if (measure=="auc"){
      res <- unlist(performance(prediction(my2CV.predict, my2CV.labels), "auc")@y.values)
    } else {
      res <- t(sapply(1:nfold, function(i){
        cr   <- cor(my2CV.predict[[i]], my2CV.labels[[i]])
        rmse <- sqrt(mean((my2CV.predict[[i]] - my2CV.labels[[i]])^2))
        return(c(cr, rmse))
      }))
    }
    metrics <- cbind(lambdas, df, res)
    toreturn <- list()
    toreturn[["metrics"]] <- metrics
    toreturn[["prediction"]] <- my2CV.predict
    toreturn[["labels"]]     <- my2CV.labels
    return(toreturn)
}

#outer loop for nested CV
doubleCV.SGL <- function(y, x, grp, nfold, fam="logit", alpha=0.5, myfold=c()){

    if (length(myfold) == 0){
      avail <- 1:nrow(x)
      nsize <- ceiling(nrow(x)/nfold)
      cv.idx <- rep(nfold, nrow(x))
      for (i in 1:(nfold-1)){
        tmp <- sample(avail, nsize)
        cv.idx[tmp] <- i
        avail <- setdiff(avail, tmp)
      }
    } else {
        cv.idx <- myfold
        nfold  <- max(myfold)
    }

    my2CV.predict <- list()
    my2CV.labels  <- list()
    lambdas <- c()
    df      <- c()
    #run CV
    for(i in 1:nfold){
        message("outer fold: ",i)
        test.y <- y[cv.idx==i]
        test.x <- x[cv.idx==i,]
        train.y <- y[cv.idx!=i]
        train.x <- x[cv.idx!=i,]

        train.data <- list(x=train.x, y=train.y)

        #inner cv.loop
        icv <- cvSGL(train.data, grp, type=fam, nfold=5, alpha=alpha)
        idx.min <- which.min(icv$lldiff)
        nmax <- icv$lldiff[idx.min] + icv$llSD[idx.min]
        idx.1se <- min(which(icv$lldiff < nmax))
        iopt <- idx.1se
        message(iopt)
        plot(icv)
        tmp <- SGL(train.data,grp, type=fam, lambda=icv$fit$lambdas, alpha=alpha)

        my2CV.predict[[i]] <- predictSGL(tmp, as.matrix(test.x), iopt)
        my2CV.labels[[i]] <- test.y
        lambdas <- c(lambdas, icv$fit$lambdas[iopt])
        df      <- c(df, sum(icv$fit$beta[,iopt] != 0))
        print(cbind(my2CV.predict[[i]], my2CV.labels[[i]]))
    }
    res <- unlist(performance(prediction(my2CV.predict, my2CV.labels), "auc")@y.values)
    return(cbind(lambdas, df, res))
}


log.stability <- function(outc, x, nrep=100, alpha=0.5){

    result <- rep(0, ncol(x))
    samp <- 1:nrow(x)
    N <- ceiling(0.9 * nrow(x))
    for(j in 1:nrep){
        message(j)
        sidx <- sample(samp, N, replace=T)
        newX <- x[sidx,]
        newS <- outc[sidx]

        cvmod <- cv.glmnet(as.matrix(newX), newS, family="binomial", nfolds=5, alpha=alpha)
        predmodel <- glmnet(as.matrix(newX), newS, family="binomial", alpha=alpha, lambda=cvmod$lambda)
        idx <- which(cvmod$lambda==cvmod$lambda.1se)
        relevant <- which(predmodel$beta[,idx] != 0)
        print(relevant)
        result[relevant] <- result[relevant] + 1
    }
    result <- result / nrep
    return(result)
}



cox.stability <- function(mysurv, x, nrep=100, alpha=0.5){

    result <- rep(0, ncol(x))
    samp <- 1:nrow(x)
    N <- ceiling(0.9 * nrow(x))
    for(j in 1:nrep){
        message(j)
        sidx <- sample(samp, N, replace=T)
        newX <- x[sidx,]
        newS <- mysurv[sidx]

        cvsurv <- cv.glmnet(as.matrix(newX), newS, family="cox", nfolds=5, alpha=alpha)
        survmodel <- glmnet(as.matrix(newX), newS, family="cox", alpha=alpha, lambda=cvsurv$lambda)
        idx <- which(cvsurv$lambda==cvsurv$lambda.1se)
        relevant <- which(survmodel$beta[,idx] != 0)
        result[relevant] <- result[relevant] + 1
    }
    result <- result / nrep
    return(result)
}




#screen alpha values for a specific set of predictors
screenSet <- function(mypredictors=c()){
  if (length(mypredictors)==0)
    mypredictors <- colnames(predictors)
  tmp <- intersect(mypredictors, colnames(predictors))
  pred.use <- predictors[,tmp]
  screen <- lapply(seq(0,1,0.1), function(al){
    message(al);
    doubleCV(myout, pred.use, nfold=5, alpha=al, myfold=myfold)
  })
  return(screen)
}

plot.screen <- function(results, toplot=c()){

  if (length(toplot)==0)
    toplot <- names(results)

  nplot <- T
  dev.null <- sapply(toplot, function(x){
    aucs <- unlist(lapply(results[[x]], function(y){ mean(y[,3]) }))
    if (nplot){
      #plot(seq(0, 1, 0.1), aucs, ylim=c(0.65,0.9), type="b", col=mycol[x], pch=mypch[x], xlab="alpha", ylab="mean AUC", las=1)
      plot(seq(0, 1, 0.1), aucs, ylim=c(0.45,0.9), type="b", col=mycol[x], pch=mypch[x], xlab="alpha", ylab="mean AUC", las=1)
      nplot <<- F
    } else {
      points(seq(0,1,0.1), aucs, type="b", col=mycol[x], pch=mypch[x])
    }
  })

  legend(0,0.9, toplot, col=mycol[toplot], pch=mypch[toplot], lwd=2)
}

plot.stab <- function(mystab, thresh=0.5){
  plot(mystab, ylim=c(0,1.1), type="h", xlab="feature", ylab="selection frequency", las=1)
  points(mystab, type="p")

  idx <- mystab >= thresh
  text(which(idx), mystab[idx], labels=names(mystab)[idx], pos=3)
  abline(h=thresh, col=2, lwd=2, lty=2)
}