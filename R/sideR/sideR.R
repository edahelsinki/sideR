## sideR
##
## Copyright 2017 Kai Puolamaki <kai.puolamaki@iki.fi>
## Copyright 2017 Finnish Institute of Occupational Health, Helsinki,
## Finland
##
## Permission is hereby granted, free of charge, to any person obtaining
## a copy of this software and associated documentation files (the
## "Software"), to deal in the Software without restriction, including
## without limitation the rights to use, copy, modify, merge, publish,
## distribute, sublicense, and/or sell copies of the Software, and to
## permit persons to whom the Software is furnished to do so, subject to
## the following conditions:
##
## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
## WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

require(fastICA)

#' Pairplot the maximal difference of the current selection and the data
#' 
#' @param data nXd data matrix
#' @param rdata nXd matrix of random data (not used)
#' @param grouped a subset of 1:n that contains the selection
#' @returns As a side effect makes a pair plot of maximum five dimensions in which
#'    the selection differs most from the overall data. Random data is not used.
#'    We use the KL-distance between univariate Gaussians as measure of difference here.
#' 
#' @export
selectionplot <- function(data,rdata=NULL,grouped=NULL) {
    if(length(grouped)>0) {
      data0 <- data[grouped,,drop=FALSE]
      data1 <- data
      m0 <- apply(data0,2,mean)
      m1 <- apply(data1,2,mean)
      if(length(grouped)>1) {
        s0 <- pmax(1e-6,apply(data0,2,var))
        s1 <- pmax(1e-6,apply(data1,2,var))
      } else {
        s0 <- s1 <- rep(1,dim(data)[2])
      }
      v <- (s0/s1+(m1-m0)^2/s1-1+log(s1/s0))/2
      p <- order(v,decreasing=TRUE)[1:min(5,length(v))]
      idx <- 1:dim(data)[1] %in% grouped
      perm <- order(idx)
      pairs(data[perm,p],
            pch=20,
            col=ifelse(idx[perm],"red","black"))
    }
}

#' Makes a scatterplot of the data
#' 
#' @param data nXd data matrix
#' @param rdata nXd matrix of random data (not used)
#' @param currentplot list that contains the current w matrix and name of the axis
#' @param grouped a subset of 1:n that contains the current selection
#' @returns As a side effect makes a pairplot to the directions specified by w.
#' 
#' @export
plotdata <- function(data,rdata,currentplot,grouped=NULL) {
  s <- currentplot$s       # PCA or ICA or axis or whatever
  w <- currentplot$w[,1:2] # Projection vectors
  v <- currentplot$v       # Projection costs
  xy <- data %*% w
  rxy <- rdata %*% w
  plot(rbind(xy,rxy),
       type="n",
       bty="n",
       xlab=sprintf("%s1[%.2g] = %s",
                    s,v[1],largestdim(w[,1],names=colnames(data))),
       ylab=sprintf("%s2[%.2g] = %s",
                    s,v[2],largestdim(w[,2],names=colnames(data))))
  idx <- 1:dim(xy)[1] %in% grouped
  
  for(i in 1:dim(xy)[1]) lines(rbind(rxy[i,],xy[i,]),col=if(i %in% grouped) "pink" else "lightgray")
  points(rxy,pch=1,col=ifelse(idx,"pink3","gray"))
  points(xy,pch=20,col=ifelse(idx,"red","black"))
  if(length(grouped)>0) {
    plotellipse(rxy[grouped,,drop=FALSE],col="blue",lty="dotted",lwd=3)
    plotellipse(xy[grouped,,drop=FALSE],col="blue",lwd=3)
  }
}

#' Makes an ellipse of data
#' 
#' @param data nXd data matrix
#' @param sigma the size of ellipse (default 0.95 confidence area)
#' @param ... graphics parameters
#' @returns As a side effect makes an 95% confidence ellipse based on data matrix.
#' 
#' @export
plotellipse <- function(data,sigma=qnorm(0.975),...) {
  data <- as.matrix(data[,1:2,drop=FALSE])
  s <- svd(mycov(data))
  x <- seq(from=0,to=2*pi,length.out=100)
  lines((rep(1,length(x)) %o% colMeans(data))
        +sigma*sqrt(s$d[1])*(cos(x) %o% s$u[,1])
        +sigma*sqrt(s$d[2])*(sin(x) %o% s$u[,2]),...)
}


#' A version of diag that always returns a matrix with x on diagonal
#' 
#' @param x a vector of length n
#' @returns A nXn matrix with x on diagonal.
#' 
#' @export
diagm <- function(x) { if(length(x)==0) matrix(,0,0) else if(length(x)==1) matrix(x,1,1) else diag(x) }

#' Classic whitening: transform data so that the covariance matrix is a unit matrix
#' 
#' @param X data matrix
#' @return Whitened data.
#' 
#' @export
whitening0 <- function(X) {
  X <- as.matrix(X)
  mu <- colMeans(X)
  s <- svd(cov(X))
  t(s$u %*% diagm(1/sqrt(s$d)) %*% t(s$u) %*% t(X-((rep(1,dim(X)[1]) %o% mu))))
}

mycov <- function(x) {
  x <- x-(rep(1,dim(x)[1]) %o% colMeans(x))
  (t(x) %*% x)/dim(x)[1]
}

#' Makes a cluster constraint with mean and variance constraints.
#'
#' @param data Data matrix
#' @param I items
#' @return cluster constraint
#'
#' @export
clusterconstraint <- function(data,I,maxd=dim(data)[2],eps=.Machine$double.eps^0.25) {
  if(length(I)<1) stop("clusterconstraint: too small set.")
  s <- svd(mycov(data[I,,drop=FALSE]))
  p <- order(s$d,decreasing=TRUE)
  w <- s$u[,p,drop=FALSE]
  ## If there are zero variance dimensions we slow down a lot due to irrelavant constraints.
  ## If variance to some direction is zero then it is usually due to finite number of data points
  ## and not because the distribution would actually be zero to the particular direction.
  # maxd <- min(maxd,sum(s$d>eps))
  unlist(c(lapply(1:maxd,function(i) list(linconstraint(data,I,w[,i]),quadconstraint(data,I,w[,i])))),
         recursive=FALSE)
}

#' Makes a column constraint with mean and variance constraints.
#'
#' @param data Data matrix
#' @param I items
#' @return cluster constraint
#'
#' @export
columnconstraint <- function(data,I=1:dim(data)[1],maxd=dim(data)[2]) {
  w <- diag(maxd)
  unlist(c(lapply(1:maxd,
            function(i) list(linconstraint(data,I,w[,i]),quadconstraint(data,I,w[,i])))),
         recursive=FALSE)
}


#' Makes a character string of a numeric vector
#'
#' @param x numeric vector
#' @return string representation of the vector
#'
#' @export
largestdim <- function(x,names=as.character(1:length(x)),s="") {
  p <- order(abs(x),decreasing=TRUE)[1:min(5,length(x))]
  paste(mapply(function(y,n) sprintf("%+.2f (%s)",y,n),x[p],names[p]),
        collapse=" ")
}

#' Makes a selectionplot parameters, to be fed to plotdata.
#' 
#' @param data nXd data matrix
#' @param grouped a subset of 1:n that contains the current selection
#' @returns A data structure that contains the projection vectors etc.
#' 
#' @export
doselection <- function(data,grouped) {
  if(length(grouped)>0) {
    data0 <- data[grouped,,drop=FALSE]
    data1 <- data
    m0 <- apply(data0,2,mean)
    m1 <- apply(data1,2,mean)
    
    if(length(grouped)>1) {
      s0 <- pmax(1e-6,apply(data0,2,var))
      s1 <- pmax(1e-6,apply(data1,2,var))
    } else {
      s0 <- s1 <- rep(1,dim(data)[2])
    }
  } else {
    m0 <- apply(data,2,mean)
    s0 <- apply(data,2,var)
    m1 <- rep(0,dim(data)[2])
    s1 <- rep(1,dim(data)[2])
  }
  v <- (s0/s1+(m1-m0)^2/s1-1+log(s1/s0))/2
  p <- order(v,decreasing=TRUE)
  list(w=diag(dim(data)[2])[,p],v=v[p])
}

#' Performs PCA on matrix Y in which components with unit variance are
#' "non-interesting"
#' @param Y "whitened" data
#' @return Returns a list of PCA components in decreasing order with
#'     respect to difference to variance 1.
#'
#' @export
dopca <- function(Y) {
  s <- svd(cov(Y))
  w <- s$u
  ## normalize w so that the mean of PCA vectors is positive
  w <- w*(rep(1,dim(w)[1]) %o%
            (apply(w,2,function(x) {
              if(mean(x)>=0) {
                1/sqrt(sum(x^2))
              } else {
                -1/sqrt(sum(x^2))
              }})))
  ## KL divergence between two normal distributions with zero mean
  ## and the other having a unit variance
  v <- (s$d-log(s$d)-1)/2
  p <- order(abs(v),decreasing=TRUE)
  list(w=w[,p],v=v[p])
}

#' Performs fastICA on matrix Y
#' @param Y "whitened" data
#' @param g ICA function
#' @param gv Expected value of g(x) if x is a Gaussian variable
#' @return Returns a list of ICA vectors in decreasing order with
#'     respect to the weight and the weights.
#'
#' @export
dofastica <- function(Y,
                      g=function(u) log(cosh(u)),
                      gv=0.3745,
                      iter=1) {
  ##                gv=mean(g(rnorm(1000000)))) {
  if(iter==1) {
    a <- fastICA(Y,dim(Y)[2])
    w <- a$K %*% a$W
    ## normalize w so that the mean of ICA vectors is positive
    w <- w*(rep(1,dim(w)[1]) %o%
              (apply(w,2,function(x) {
                if(mean(x)>=0) {
                  1/sqrt(sum(x^2))
                } else {
                  -1/sqrt(sum(x^2))
                }})))
    ##v <- apply(g(a$S),2,mean)-gv
    v <- apply(g(whitening0(Y) %*% w),2,mean)-gv
    p <- order(abs(v),decreasing=TRUE)
    list(w=w[,p],v=v[p])
  } else {
    res <- lapply(1:iter,function(x) dofastica(Y,g=g,gv=gv,iter=1))
    res[[which.max(sapply(res,
                          function(x) sort(abs(x$v),
                                           decreasing=TRUE)[2]))]]
  }
}

#' Build a new function with a smaller environment, see
#' http://www.win-vector.com/blog/2015/03/using-closures-as-objects-in-r/
#' @param f input function
#' @param vars names we are allowing to be captured in the closere
#' @return New function with closure restricted to variables in vars.
#'
#' @export
restrictEnvironment <- function(f,vars) {
  oldEnv <- environment(f)
  newEnv <- new.env(parent=parent.env(oldEnv))
  for(v in vars) {
    assign(v,get(v,envir=oldEnv),envir=newEnv)
  }
  environment(f) <- newEnv
  f
}


#' Auxiliary function that indexes a vector.
#'
#' @param rows idx of integers (e.g., row indices related to values),
#'     containing values in [1,m] where m is the number of constraints
#' @return A list whose i'th element is a list containing the indices
#'     of values related to the i'th constraint. So this kind of
#'     inverts a list of indices.
#'
#' @export
findi <- function(idx,maxidx=max(unlist(idx))) {
  a <- unlist(idx)
  b <- unlist(lapply(1:length(idx),
                     function(i) {
                       if(length(idx[[i]])==0)
                         NULL
                       else
                         rep(i,length(idx[[i]]))
                     }))
  p <- order(a)
  a <- a[p]
  b <- b[p]
  start <- c(which(!duplicated(a)),length(a)+1)
  
  res <- vector("list",maxidx)
  res[a[start[-length(start)]]] <- lapply(1:(length(start)-1),
                                          function(i)
                                            sort(unique(b[start[i]:(start[i+1]-1)])))
  res
}

#' Finds the value of a constraint
#' @param data data frame
#' @param fcon constraint function
#' @return Value of the constraint
#'
#' @export
compval <- function(data,fcon) sum(apply(data,1,fcon))


########################################################################
## Definitions related to Gaussian model. Replace them with something
## else for some other exponential family model.

fcon_lin <- function(x,w) sum(x*w)

fcon_quad <- function(x,w,m=0) (sum(x*w)-m)^2

linconstraint <- function(data,I,w) {
  list(
    I=I,
    val=compval(data[I,,drop=FALSE],function(x) fcon_lin(x,w)),
    w=w,
    type="l")
}

quadconstraint <- function(data,I,w,mhat=TRUE) {
  m <- if(mhat) mean(data[I,,drop=FALSE] %*% w) else 0
  list(
    I=I,
    val=compval(data[I,,drop=FALSE],function(x) fcon_quad(x,w,m)),
    w=w,
    m=m,
    type="q")
}

#' Main workhorse function that takes in constraints and produces whitened data and randomized data.
#' 
#' @param data nXd data matrix
#' @param constraints list of constraints
#' @param nn n
#' @param nd d
#' @param maxiter maximum number of iterations
#' @param maxtime time cutoff in seconds
#' @param tol1 tolerance parameter on changes of lambda
#' @param tol2 tolerance parameter on changes of test statistic
#' @return List containing whitened data Y and randomized data rdata.
#' 
#' @export
getwhitenedR <- function(data,
                         constraints,
                         nn=dim(data)[1],
                         nd=dim(data)[2],
                         maxiter=1000,maxtime=10,
                         tol1=0.01,
                         tol2=0.01) {
  cat("getwhitenedR ")
  starttime <- as.numeric(Sys.time())
  X <- sideR(constraints=constraints,data=data,nn=nn,nd=nd,epstol=tol1,tol=tol2*sd(data))
  count <- 1
  aux <- X$updateall(1)
  while(aux<1) {
    cat(".")
    count <- count+1
    aux <- X$updateall(1)
    maxiter <- maxiter-1
    if(maxiter<=0) {
      cat(" maxiter reached")
      break
    }
    if(as.numeric(Sys.time())-starttime>=maxtime) {
      cat(" timeout reached")
      break
    }
  }
  cat(sprintf(" %d iterations, %.1f seconds.\nsigma12 ",count,as.numeric(Sys.time())-starttime))
  X$sigma12() ## Computes matrix needed for whitening and sampling, O(nd^3) operation.
  cat("whitening(Y) ")
  Y <- X$whitening(data)
  cat("sampling(rdata) ")
  rdata <- X$sample()
  cat("done\n")
  
  list(Y=Y,rdata=rdata)
}

#' Main method
#' 
#' @param constraints list of constraints
#' @param data nXd data matrix
#' @param nn n
#' @param nd d
#' @param debug boolean, to have extended debugging output
#' @param eps epsilon parameter
#' @param epstol tolerance parameter on changes of lambda
#' @param tol tolerance parameter on changes of test statistic
#' @param qprior softening prior
#' @return List of methods.
#' 
#' @export
sideR <- function(constraints,data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE,
                  eps=.Machine$double.eps^0.25,epstol=0.01,tol=0.01*sd(data),
                  qprior=0) { 
  rows <- lapply(constraints,function(x) x$I)
  irows <- findi(rows,maxidx=nn)
  urows <- unique(irows)
  idxurows <- sapply(irows,function(x) which(sapply(urows,function(y) identical(x,y)))[1])
  iurows <- findi(idxurows)
  niurows <- sapply(iurows,length)
  
  ww <- sapply(constraints,function(x) x$w,simplify="matrix")
  
  ntheta <- length(urows)
  vals <- sapply(constraints,function(x) x$val)
  linc <- sapply(constraints,function(x) x$type=="l")
  mm <- rep(NA,length(constraints))
  mm[!linc] <- sapply(constraints[!linc],function(x) x$m)
  
  vals[!linc] <- vals[!linc]+qprior*sapply(rows[!linc],length)
  
  ## Initially set everything to zero mean and unit variance
  imu <- mu <- matrix(0,nd,ntheta)
  isigma <- sigma <- sigma12 <- sigma12i <- array(diag(nd),dim=c(nd,nd,ntheta))
  
  updatelin <- function(i) {
    if(constraints[[i]]$type!="l") stop("updatelin: wrong type")
    aux <- tabulate(idxurows[rows[[i]]],nbins=ntheta)
    valshat <- sum(ww[,i]*(mu %*% aux))
    sigmasum <- sum(ww[,i]*(apply(sigma,3,function(x) x %*% ww[,i]) %*% aux))
    lambda <- (vals[i]-valshat)/sigmasum
    imu[,aux>0] <<- imu[,aux>0]+lambda*(ww[,i] %o% rep(1,sum(aux>0)))
    for(j in which(aux>0)) mu[,j] <<- imu[,j] %*% sigma[,,j]
    valshatnew <- sum(ww[,i]*(mu %*% aux))
    if(debug) {
      cat(sprintf("updatelin: valshat = %g -> %g val = %g lambda = %g\n",valshat,valshatnew,vals[i],lambda))
    }
    c(lambda,abs(valshat-valshatnew)/length(rows[[i]]))
  }
  
  updatequad <- function(i) {
    if(constraints[[i]]$type!="q") stop("updatequad: wrong type")
    aux <- tabulate(idxurows[rows[[i]]],nbins=ntheta)
    valshat <- (sum(ww[,i]*(apply(sigma,3,function(x) x %*% ww[,i]) %*% aux)) 
                +sum(aux*((ww[,i] %*% mu)^2+mm[i]^2-2*mm[i]*(ww[,i] %*% mu))))
    cc <- sapply(which(aux>0),function(j) sum(ww[,i]*(sigma[,,j] %*% ww[,i])))
    dd <- sapply(which(aux>0),function(j) sum(imu[,j]*(sigma[,,j] %*% ww[,i])))
    ee <- c(ww[,i] %*% mu[,aux>0])
    # f <- function(x) {
    #   ll <- x/(1+x*cc)
    #   sum(aux[aux>0]*(ll*cc^2+2*ll*cc*dd*ee-ll^2*cc^2*dd^2))+vals[i]-valshat
    # }
    f <- function(x) {
      ll <- x/(1+x*cc)
      ff <- x*mm[i]-ll*dd-x*mm[i]*ll*cc
      sum(aux[aux>0]*(ll*cc^2-ff^2*cc^2+2*ff*cc*(mm[i]-ee)))+vals[i]-valshat
    }
    a <- (eps-1)/max(cc)
    fa <- f(a)
    b <- 1/eps
    fb <- f(b)
    if(sign(fa)*sign(fb)>=0) {
      lambda <- if(abs(fa)<abs(fb)) a else b
    } else {
      fc <- f(0)
      if(sign(fa)*sign(fc)<0) {
        b <- 0
        fb <- fc
      } else {
        a <- 0
        fa <- fc
      }
      lambda <- uniroot(f,lower=a,upper=b,f.lower=fa,f.upper=fb)$root
    }
    for(j in which(aux>0)) {
      isigma[,,j] <<- isigma[,,j]+lambda*(ww[,i] %o% ww[,i])
      imu[,j] <<- imu[,j]+lambda*mm[i]*ww[,i]
      bb <- c(sigma[,,j] %*% ww[,i])
      sigma[,,j] <<- sigma[,,j]-(bb %o% bb)*(lambda/(1+lambda*sum(ww[,i]*bb)))
      mu[,j] <<- imu[,j] %*% sigma[,,j]
    }
    valshatnew <- (sum(ww[,i]*(apply(sigma,3,function(x) x %*% ww[,i]) %*% aux)) 
                   +sum(aux*((ww[,i] %*% mu)^2+mm[i]^2-2*mm[i]*(ww[,i] %*% mu))))
    if(debug) {
      cat(sprintf("updatequad: valshat = %g -> %g val = %g lambda = %g\n",valshat,valshatnew,vals[i],lambda))
    }
    c(lambda,abs(sqrt(valshat)-sqrt(valshatnew))/length(rows[[i]]))
  }
  
  updateall <- function(iter=100,randomorder=FALSE) {
    loss1 <- loss2 <- 0
    while(iter>0) {
      aux <- sapply(if(randomorder) sample.int(length(constraints)) else 1:length(constraints),
                    function(i) if(constraints[[i]]$type=="l") updatelin(i) else updatequad(i),
                    simplify="matrix")
      loss1 <- max(abs(aux[1,]))
      loss2 <- max(aux[2,])
      if(debug) {
        cat("updateall: ")
        for(x in aux[1,]) cat(sprintf("%.2g ",x))
        cat("\n")
        cat(sprintf("updateall: iter = %d loss = %g tol = %g\n",iter,loss1,loss2))
      }
      if(loss1<epstol || loss2<tol) break
      iter <- iter-1
    }
    iter
  }
  
  list(
    sigma12=function() {
      for(i in 1:ntheta) {
        s <- svd(sigma[,,i]) ## O(nd^3)!!
        sigma12[,,i] <<- (sqrt(s$d) %o% rep(1,nd))*t(s$u) 
        sigma12i[,,i] <<- s$u %*% ((sqrt(1/s$d) %o% rep(1,nd))*t(s$v)) ## O(nd^3)!!
      }
    },
    sample=function(random=TRUE) {
      rdata <- t(mu[,idxurows])
      if(random) {
        for(i in 1:ntheta) {
          idx <- iurows[[i]]
          rdata[idx,] <- rdata[idx,]+matrix(rnorm(length(idx)*nd),length(idx),nd) %*% sigma12[,,i]
        }
      }
      rdata
    },
    whitening=function(data,scale=TRUE) {
      Y <- data-t(mu[,idxurows]) ## Data should be zero centered
      if(scale) {
        for(i in 1:ntheta) {
          idx <- iurows[[i]]
          Y[idx,] <- Y[idx,,drop=FALSE] %*% sigma12i[,,i]
        }
      }
      Y
    },
    updateall=updateall,
    updatelin=updatelin,
    updatequad=updatequad,
    list=function() {
      list(idxurows=idxurows,irows=irows,mu=mu,isigma=isigma,sigma=sigma,sigma12=sigma12,sigma12i=sigma12i)
    }
  )
}
