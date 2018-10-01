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


source("../sideR/sideR.R")

set.seed(42)

gensyntdata <- function(n=1000,k=5,d=10,sd=0.5) {
  if(n<2*k) stop("gensyntdata: too small k.\n")
  centroids <- matrix(rnorm(k*d,sd=1),nrow=k,ncol=d)
  class <- sample(c(1:k,1:k,sample.int(k,size=n-2*k,replace=TRUE)))
  list(data=centroids[class,]+matrix(rnorm(n*d,sd=sd),nrow=n,ncol=d),class=class)
}

bencmarkconstraints <- function(data,constraints) {
  t1 <- system.time({ X <- sideR(constraints=constraints,data=data,
                                 nn=dim(data)[1],nd=dim(data)[2],debug=FALSE) 
  },gcFirst=TRUE)["elapsed"]
  t2 <- system.time({ X$updateall(100) },gcFirst=TRUE)["elapsed"]
  t3 <- system.time({ X$sigma12() },gcFirst=TRUE)["elapsed"]
  t4 <- system.time({ Y <- X$whitening(data) },gcFirst=TRUE)["elapsed"]
  t5 <- system.time({ rdata <- X$sample() },gcFirst=TRUE)["elapsed"]
  t6 <- system.time({ pcaplot <- dopca(Y) },gcFirst=TRUE)["elapsed"]
  t7 <- system.time({ pcaplot <- dofastica(Y) },gcFirst=TRUE)["elapsed"]
  data.frame(init=t1,updateall=t2,sigma12=t3,whitening=t4,sample=t5,pca=t6,ica=t7)
}

benchmarkdata <- function(n,k,d) {
  aux <- gensyntdata(n,k,d)
  data <- aux$data
  class <- aux$class
  constraints <- columnconstraint(data)
  if(k>1) {
    for(i in 1:k) {
      constraints <- c(constraints,clusterconstraint(data,which(class==i)))
    }
  }
  bencmarkconstraints(data,constraints)
}

dobm <- function(nn=3,nk=3,nd=4,i=1) {
  ns <- 1024*2^seq.int(nn)
  ks <- c(1,2^seq.int(nk))
  ds <- 8*2^seq.int(nd)
  
  res <- NULL

  idx <- 1
  for(n in ns) {
    for(k in ks) {
      for(d in ds) {
        if(d<n) {
          cat(sprintf("n = %d k = %d d = %d i = %d ",n,k,d,i))
          aux <- data.frame(cbind(data.frame(n=n,k=k,d=d,i=i),benchmarkdata(n,k,d)),row.names=idx)
          junk <- aux[,-(1:4)]
          for(j in colnames(aux)[-(1:4)]) {
            cat(sprintf("%s = %g ",j,aux[,j]))
          }
          cat("\n")
          if(is.null(res)) {
            res <- aux
          } else {
            res <- rbind(res,aux)
          }
          idx <- idx+1
        }
      }
    }
  }
  res
}

for(i in 1:10) {
  set.seed(42+i)
  aux <- dobm(i=i)
  saveRDS(aux,file=sprintf("runtime-alternate-%02d.rds",i))
}
