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

checkonv <- function(X,iter=50) {
  res <- matrix(NA,iter,12)
  j <- X$list()$idxurows
  for(i in 1:iter) {
    for(k in 1:length(j)) {
      for(kk in 1:2) {
        res[i,(j[k]-1)*2+kk] <- X$list()$sigma[kk,kk,j[k]]
        res[i,6+(j[k]-1)*2+kk] <- X$list()$mu[kk,j[k]]
      }
    }
    X$updateall(1)
  }
  res
}

plotellipse2 <- function(xx,a,b,...) {
  x <- seq(from=0,to=2*pi,length.out=100)
  lines((rep(1,length(x)) %o% xx)
        +a*(cos(x) %o% c(1,0))
        +b*(sin(x) %o% c(0,1)),...)
}



data <- matrix(c(1,0,0,0,1,0),3,2)

constraint13 <- list(linconstraint(data,c(1,3),c(1,0)),
                     quadconstraint(data,c(1,3),c(1,0)),
                    linconstraint(data,c(1,3),c(0,1)),
                     quadconstraint(data,c(1,3),c(0,1)))

constraint23 <- list(linconstraint(data,c(1,3),c(0,1)),
                     quadconstraint(data,c(2,3),c(1,0)),
                     linconstraint(data,c(2,3),c(0,1)),
                     quadconstraint(data,c(2,3),c(1,0)))

constraints <- c(constraint13,constraint23)
                    

X13  <- sideR(constraint13,data)
X123 <- sideR(c(constraint13,constraint23),data)

X123 <- sideR(constraints,data)

y13 <- checkonv(X13,iter=1000)[,1]
y123 <- checkonv(X123,iter=1000)[,1]

op <- par(no.readonly=TRUE)


pdf("toyconv.pdf")
par(mar=op$mar+c(0,2,0,0))
x <- 1:length(y13)
plot(c(x,x),c(y13,y123),type="n",bty="n",log="xy",cex.lab=1.5,
     xlab="iterations",ylab=expression((Sigma[1])[11]))
points(x,y13)
lines(x,y13,lty="dotted")
points(x,y123,pch=20,col="red")
lines(x,y123,col="red")
legend("bottomleft",legend=c("Case A","Case B"),lty=c("dotted","solid"),pch=c(1,20),col=c("black","red"))
dev.off()

par(op)

pdf("toy2d.pdf")
par(mar=op$mar+c(0,2,0,0))
eps <- 0.05
plot(c(-0.2,1.2),c(-0.2,1.2),type="n",bty="n",cex.lab=1.5,xlab=expression(x[1]),ylab=expression(x[2]))
text(data,labels=1:3)
plotellipse2(c(0.5,0),a=0.5+eps,b=eps,col="black")
plotellipse2(c(0,0),a=1+eps,b=1+eps,col="black",lty="dashed")
plotellipse2(c(0,1),a=eps,b=eps,col="red")
plotellipse2(c(0,0),a=eps,b=eps,col="red")
plotellipse2(c(1,0),a=eps,b=eps,col="red")
legend("topright",legend=c("Case A, rows 1 and 3","Case A, row 2","Case B"),
       col=c("black","black","red"),lty=c("solid","dashed","solid"))
dev.off()

## junk <- data.frame(i=1:length(y123),y13=y13,y123=y123)

par(op)

