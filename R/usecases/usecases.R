## Examples of using sideR 
##
## Copyright 2017 Emilia Oikarinen <emilia.oikarinen@gmail.com>
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

set.seed(42)
source("../sideR/sideR.R")

## Custom plotting functions

myPalette <- c("#009E73",  "#0072B2", "#CC79A7", "#e79f00","#F0E442", "#9ad0f3")
pch <- 16

plotpairs <- function(data,class=NULL,legend=NULL) {
  n <- dim(data)[1]
  
  panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    
    nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    points((breaks[-1]+breaks[-nB])/2,y)
    lines((breaks[-1]+breaks[-nB])/2,y)
  }
  
  if(!is.null(class)) {
    cols <- myPalette[class] 
  } else {
    cols <- c(rep("black",n))
  }
  
  if(!is.null(legend)) {
    pairs(data, oma=c(8,4,4,4),
          diag.panel=panel.hist,
          col=cols,
          pch=pch)
    par(xpd = TRUE)
    legend("bottom", fill=myPalette, legend=legend, horiz=TRUE, bty = "n")
  }
  else {
    pairs(data, 
          diag.panel=panel.hist,
          col=cols,
          pch=pch)
  }
}

plotdata <- function(data,rdata,currentplot,grouped=NULL,rdata0=NULL,rounding=0) {
  s <- currentplot$s       # PCA or ICA or axis or whatever
  w <- currentplot$w[,1:2] # Projection vectors
  v <- currentplot$v       # Projection costs
  xy <- data %*% w
  rxy <- rdata %*% w
  
  xlim <- NULL
  ylim <- NULL
  if(length(rdata0)>0) {
    rxy0 <- rdata0 %*% w
    xlim <- c(min(rxy0[,1]),max(rxy0[,1]))
    ylim <- c(min(rxy0[,2]),max(rxy0[,2]))
  }
  
  if(rounding!=0) {
    temp <- round(xy / rounding) 
    idxs <- which(!duplicated(temp))
    grouped <- grouped[grouped %in% idxs]
  } 
  else {
    idxs <- 1:dim(xy)[1]
  }
  
  plot(rbind(xy,rxy),xlim=xlim,ylim=ylim,
       type="n",
       bty="n",
       xlab=sprintf("%s1[%.2g] = %s",
                    s,v[1],largestdim(w[,1],names=colnames(data))),
       ylab=sprintf("%s2[%.2g] = %s",
                    s,v[2],largestdim(w[,2],names=colnames(data))))
  
  idx <- idxs %in% grouped
  
  for(i in idxs) lines(rbind(rxy[i,],xy[i,]),col=if(i %in% grouped) "pink" else "lightgray")
  points(rxy[idxs,],pch=1,col=ifelse(idx,"pink3","gray"))
  points(xy[idxs,],pch=20,col=ifelse(idx,"red","black"))
  if(length(grouped)>0) {
    plotellipse(rxy[grouped,,drop=FALSE],col="blue",lty="dotted",lwd=3)
    plotellipse(xy[grouped,,drop=FALSE],col="blue",lwd=3)
  }
  
  dim(xy)[1]/length(idxs)
}

# Toydata X_5
# 
# The toydata X_5 consists of 1000 points. In the first 3 dimensions there are four clusters 
# such that in any 2D-projection of X1, X2, and X3, two of these overlap. In the remaining 
# dimensions, there is cluster stucture of three clusters (randomly chosen centers, and some 
# dependency on the clusters in first three dimensions: with 75% probability points from 
# clusters 2-4(X1-X3) are assigned randomly to clusters 1-2(X4-X5), else the point is assigned 
# to the cluster 3(X4-X5) which has a higher variance).

datasize <- 1000
centroids <- matrix(data=c(0,.75,0,0,0,0,1.5,0,0,0,0,1),nrow=4,ncol=3)
centroids3 <- matrix(data=3*rnorm(6),nrow=3,ncol=2)

# Function for generating points  
genpoint <- function(centroids,n=12) {
  res <- rep(NA,5)
  i <- sample.int(n=dim(centroids)[1],size=1)
  res[1:dim(centroids)[2]] <- centroids[i,]+rnorm(n=dim(centroids)[2],mean=0,sd=0.1)
  
  if (rnorm(1)<0.75 & i>1) {
    j <- sample.int(n=dim(centroids3)[1]-1, size=1)
    res[4:5] <- centroids3[j,]+rnorm(n=dim(centroids3)[2],mean=0,sd=0.1)
  } else {
    j <- 3
    res[4:5] <- centroids3[j,]+rnorm(n=dim(centroids3)[2],mean=0,sd=0.3)
  }
  c(res,i,j)
}

data <- t(replicate(datasize,genpoint(centroids)))
data[,c(-6,-7)] <- scale(data[,c(-6,-7)])
toydata <- data.frame(data)
toydata[,6] <- as.factor(toydata[,6])
toydata[,7] <- as.factor(toydata[,7])
data <- data[,1:5]
colnames(data)<- c('X1', 'X2','X3','X4','X5')

##  Pairplot of X_5
# Colours are used to indicate the cluster identities of the four clusters 
# in the first three dimensions.

idx <- sample.int(dim(data)[1],250)
pdf("toy5.pdf")
plotpairs(data[idx,],toydata[idx,"X6"], legend=c("A","B","C","D"))
dev.off()

## First projection to data, pairplot of whitened data
# Next, we add the column constraints (this is always the first step) 
# and compute then the whitened data as well as the projection to 2D 
# using ICA to find most prominent directions in the whitened data.
#
# In the projection, the user can observe the four clusters in dimensions 1-3. 

# Add the column constraints
constraints <- columnconstraint(data)
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)

# Do whitening and perform ICA on whitened data
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
icaplot0 <- dofastica(Y)
rdata0 <- rdata

png("toy5-ica0.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="ICA"),icaplot0))
dev.off()

pdf("toy5-white0.pdf")
plotpairs(Y[idx,])
dev.off()

# After adding four cluster constraints (one for each of the four clusters visible), 
# the whitened data showns further structure in dimensions 4-5, i.e., the user 
# can now observe the three clusters in dimensions 4-5. 

# Here we use directly cluster identities (for simplicity of the code) to add the 
# four visible cluster, but one should notice that these sets of points could be easily 
# handpicked from the plot without using the cluster identities

class <- toydata$X6
for(i in 1:4) {
  constraints <- c(constraints,clusterconstraint(data,which(class==i)))
}
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)

# Do whitening and perform ICA on whitened data
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
icaplot1 <- dofastica(Y)

png("toy5-ica0_after.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="ICA"),icaplot0),rdata0=rdata0)
dev.off()

png("toy5-ica1.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="ICA"),icaplot1))
dev.off()

pdf("toy5-white1.pdf")
plotpairs(Y[idx,])
dev.off()

# After adding the cluster constraints for the three clusters visible in the 
# view above, the whitened data shownsno further clear structure, and the random 
# sample from the background distribution has the same structure in the most 
# prominent ICA projection as well.

# Add three cluster constraints
class <- toydata$X7
for(i in 1:3) {
  constraints <- c(constraints,clusterconstraint(data,which(class==i)))
}
X <- sideR(constraints=constraints,data=data, nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)

# Do whitening and perform ICA on whitened data
u <- X$updateall(100);
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
icaplot2 <- dofastica(Y)

png("toy5-ica2.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="ICA"),icaplot2))
dev.off()

pdf("toy5-white2.pdf")
plotpairs(Y[idx,])
dev.off()

## 3-dimensional toydata 

# Generation of data
data <- matrix(c(rep(c(0,-1,0),50),rep(c(-1,0,0),50),
                 rep(c(1.01,1, 0.5),25),
                 rep(c(1,1.01,-0.5),25)),nrow=150,ncol=3,byrow=TRUE)
data <- data+0.1*rnorm(450)

toydata <- data.frame(data)
toydata[,4] <- as.factor(c(rep(1,50),rep(2,50),rep(3,25),rep(4,25)))
colnames(data)<- c('X1', 'X2','X3')
colnames(toydata)<- c('X1', 'X2','X3','X4')

# Addition of column constraints
constraints <- columnconstraint(data)
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata0 <- X$sample()
# PCA of the whitened data
pcaplot0 <- dopca(Y)

# The most informative view to data
png("intro_ex1.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata0,c(list(s="PCA"),pcaplot0))
dev.off()

# Cluster constraints are added for the three clusters visible in the current view 
# and background distribution updated
class <- toydata$X4
for(i in 1:2) {
  constraints <- c(constraints,clusterconstraint(data,which(class==i)))
}
constraints <- c(constraints,clusterconstraint(data,c(which(class==3),which(class==4))))

X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
# PCA of the whitened data
pcaplot1 <- dopca(Y)

png("intro_ex2.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot0),rdata0=rdata0)
dev.off()

png("intro_ex3.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot1))
dev.off()

# Further cluster constraints can be added for the two clusters now separated in the view above.
# Background distribution is updated.
for(i in 3:4) {
  constraints <- c(constraints,clusterconstraint(data,which(class==i)))
}
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()

## Use case: British National Corpus data

bnc <- readRDS("../sideR/data/bnc.rds")
data <- as.matrix(bnc[-1])

# We always start by adding the column constraints
constraints <- columnconstraint(data)
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
pcaplot0 <- dopca(Y)

xy1 <- data %*% pcaplot0$w[,1]
xy2 <- data %*% pcaplot0$w[,2]

# The first selection of points 
selection1 <- which(xy1>35.9 & xy2>0 ) 

png("bnc0.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot0), selection1)
dev.off()

# A cluster constraint is added for the selection and background distribution updated
constraints <- c(constraints,clusterconstraint(data,selection1))
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
pcaplot1 <- dopca(Y)

# The second selection
xy1 <- data %*% pcaplot1$w[,1]
xy2 <- data %*% pcaplot1$w[,2]
selection2 <- which(xy1< -1)

png("bnc1.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot1), selection2)
dev.off()

# A cluster constraint is added for the second selection
constraints <- c(constraints,clusterconstraint(data,selection2))
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
pcaplot2 <- dopca(Y)

png("bnc2.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot2))  
dev.off()

## Use case: UCI Segmentation data

segment <- readRDS("../sideR/data/segment.rds")
colnames(segment)<-c("CLASS", "vedge","vedge.sd", "hedge","hedge.sd", "R", "B", "G", "value", "saturation", "hue")
data <- as.matrix(segment[-1])

# We always start by adding the column constraints
constraints <- columnconstraint(data)
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
pcaplot0 <- dopca(Y)

png("segment0.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot0))
dev.off()

# We add 1-cluster constraint
constraints <- c(constraints,clusterconstraint(data,c(1:dim(data)[1])))
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()
rdata0 <- rdata

xy1 <- data %*% pcaplot0$w[,1]
xy2 <- data %*% pcaplot0$w[,2]

# Selections for clustered sets of points
selection1 <- which(xy1>1) 
selection2 <- which((xy1<1 & xy2 >-2.5 & xy1>-0.5) | (xy1> -2.3 & xy2 >-0 & xy1<=-0.5) | (xy1> -1.5 & xy2 >-4 & xy2<=0 & xy1<=-0.5)) 
outlier <-  which(xy2> 0 & xy1 < -3)
selection3 <- setdiff(c(1:dim(data)[1]),c(selection1, selection2, outlier))

png("segment1.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot0), selection1)
dev.off()

png("segment2.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot0), selection2)
dev.off()

png("segment3.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot0), selection3)
dev.off()

# CLuster constraint is added for all three selections and the background distribution updated
constraints <- c(constraints,clusterconstraint(data,selection1), clusterconstraint(data,selection2),         
                 clusterconstraint(data,selection3))
X <- sideR(constraints=constraints,data=data,nn=dim(data)[1],nd=dim(data)[2],debug=FALSE)
u <- X$updateall(100)
s <- X$sigma12()
Y <- X$whitening(data)
rdata <- X$sample()

png("segment4.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot0),rdata0=rbind(data,rdata0))
dev.off()

# The most informative view is updated
pcaplot1 <- dopca(Y)
png("segment5.png", width = 7, height = 7, units = 'in', res = 300)
p <- plotdata(data,rdata,c(list(s="PCA"),pcaplot1))
dev.off()