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

outputfile <- "output-alternate.txt"
sgroup <- function(x,n) {
  s <- sprintf("$\\{%.1f",x[1])
  for(i in 2:n) s <- if(i<=length(x)) sprintf("%s,%.1f",s,x[i]) else sprintf("%s,-",s)
  sprintf("%s\\}$",s)
}
slist <- function(x,n) {
  s <- ""
  for(i in 1:(length(x)-1)) s <- sprintf("%s%s&",s,sgroup(x[[i]],n))
  sprintf("%s%s",s,sgroup(x[[length(x)]],n))
}

aux2 <- Reduce("rbind",lapply(list.files(pattern="^runtime-alternate-..\\.rds$"),function(x) readRDS(x)))
res2 <- aggregate(aux2[,-(1:4)],by=aux2[,c("n","d","k")],FUN=median)
res2 <- res2[order(res2[,"k"]),]
res2 <- res2[order(res2[,"d"]),]
res2 <- res2[order(res2[,"n"]),]
res <- cbind(res2[,1:3],res2[,-(1:3)][,apply(res2[,-(1:3)],2,max)>=2])


cat("$n$&$d$",file=outputfile,append=FALSE)
for(i in colnames(res)[-(1:3)]) cat(sprintf("&%s",i),file=outputfile,append=TRUE)
cat("\\\\\n",file=outputfile,append=TRUE)
n <- d <- -1
nn <- length(unique(res[,"k"]))
for(i in 1:dim(res)[1]) {
  if(res[i,"n"]!=n || res[i,"d"]!=d) {
    if(n>1) cat(sprintf("$%d$&$%d$&%s\\\\\n",n,d,slist(aux,nn)),file=outputfile,append=TRUE)
    aux <- vector("list",dim(res)[2]-3)
    names(aux) <- colnames(res)[-(1:3)]
    n <- res[i,"n"]
    d <- res[i,"d"]
  }
  for(j in names(aux)) {
    aux[[j]] <- c(aux[[j]],res[i,j])
  }
}
cat(sprintf("$%d$&$%d$&%s\\\\\n",n,d,slist(aux,nn)),file=outputfile,append=TRUE)


