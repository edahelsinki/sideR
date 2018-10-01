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


require(shiny)
source("sideR.R")

set.seed(42)

## Use the first .rds file found in the data directory
datafile <- list.files("data/",pattern=".*\\.rds$")
kickstart <- function(datafile) {
  origdata <- readRDS(sprintf("data/%s",datafile))
  ## Shuffle rows to avoid any spurious structures in plots due to ordering of rows
  origdata <- origdata[sample.int(dim(origdata)[1]),]

  ## subsets are the initial subsets defined by the factors in the data frame
  subsets <- if(any(sapply(origdata,is.factor))) {
    unlist(
      lapply(
        lapply(
          colnames(origdata)[sapply(origdata,is.factor)],
          function(s) model.matrix(as.formula(sprintf("~0+%s",s)),origdata)),
        function(x) lapply(as.data.frame(x),function(y) which(y==1))),
      recursive=FALSE)
  } else {
    list()
  }
  subsets[["all-column"]] <- 1:dim(origdata)[1]

  ### The below variables are actually used somewhere

  ## data is the real valued part of the data
  data <- as.matrix(origdata[,!sapply(origdata,is.factor)])
  rdata <- data+0.1*matrix(rnorm(length(data)),dim(data)[1],dim(data)[2]) # create some random data...
  constraints <-columnconstraint(data)
  constraintsubsets <- rep("all-column",length(constraints))
  ## currentplot contains plot name s and nX2 projection matrix w, v being the projection gains.
  ## Initially just set it to first two axis.
  currentplot <- c(list(s="SEL"),doselection(data,NULL))
  list(datadescr=sprintf("%s (n=%d d=%d c=%d)",datafile,dim(data)[1],dim(data)[2],sum(sapply(origdata,is.factor))),
       data=data,rdata=rdata,subsets=subsets,constraints=constraints,
       constraintsubsets=constraintsubsets,currentplot=currentplot)
}

aux <- kickstart(datafile[1])
datadescr <- aux$datadescr
data <- aux$data
rdata <- aux$rdata
subsets <- subsets0 <- aux$subsets
constraints <- constraints0 <- aux$constraints
constraintsubsets <- constraintsubsets0 <- aux$constraintsubsets
currentplot <- aux$currentplot

grouped <- NULL
selectioncount <- 0

aux <- getwhitenedR(data,constraints)
rdata <- aux$rdata
Y <- aux$Y

currentplotpca <- c(list(s="PCA"),dopca(Y))
currentplotica <- c(list(s="ICA"),dofastica(Y))
currentplot <- currentplotsel <- c(list(s="PCA"),dopca(data))
currentdatasetname <- datafile[1]

