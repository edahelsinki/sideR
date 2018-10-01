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
## The global definitions are in global.R!

## The data is in the following variables:
##
## data the real data matrix
##
## subsets0 the original subsets
## 
## subsets the current subsets, list where the name is the name of the subset
##         and the contents are the row indices
## 
## grouped a set containing the current selection
##
## constraints the list containing the constraints

shinyServer(function(input,output,session) {
    updateall <- function() {
      output$datadescr <<- renderText({ 
          s <- ""
          if(length(grouped)>0) {
              aux <- sapply(subsets,function(x) length(intersect(x,grouped))/length(union(x,grouped)))
              p <- order(aux,decreasing=TRUE)
              s <- "matches = "
              for(i in p[1:min(5,length(p))]) {
                  if(aux[i]>0) {
                      s <- sprintf("%s(%s,%.3f) ",
                                   s,names(subsets)[i],aux[i])
                  }
              }
          } 
    
    
          sprintf("%s<br> current selection = %d subsets = %d constraints = %d<br> tol1 = %.2g tol2 = %.2g timeout = %ds<br> %s",
                  datadescr,length(grouped),length(subsets),length(constraints),
                  10^input$tol1,10^input$tol2,input$timeout,s) 
      })
  
      output$plot <<- renderPlot({plotdata(data,rdata,currentplot,grouped)})
  
      output$selectionplot <<- renderPlot({selectionplot(data,grouped=grouped)})
  }
  
    observeEvent(input$dataset,{
        aux <- input$dataset
        if(currentdatasetname!=aux) {
            currentdatasetname <<- aux
            aux <- kickstart(aux)
            datadescr <<- aux$datadescr
            data <<- data <- aux$data
            subsets <<- aux$subsets
            subsets0 <<- aux$subsets
            constraints <<- constraints <- aux$constraints
            constraints0 <<- aux$constraints
            constraintsubsets <<- aux$constraintsubsets
            constraintsubsets0 <<- aux$constraintsubsets
            currentplot <<- aux$currentplot
            grouped <<- NULL
            selectioncount <<- 0
    
            aux <- getwhitenedR(data,constraints,tol1=10^input$tol1,tol2=10^input$tol2,maxtime=input$timeout)
            rdata <<- aux$rdata
            Y <<- Y <- aux$Y
    
            currentplotpca <<- c(list(s="PCA"),dopca(Y))
            currentplotica <<- c(list(s="ICA"),dofastica(Y))
            currentplot <<- currentplotsel <<- c(list(s="PCA"),dopca(data))
    
            updateSelectInput(session,"selectpoints",
                              label="add to current selection:",
                              choices=c("none",names(subsets)))
        }
        updateall()
  })
  
  observeEvent(input$sel1,{
    s <- input$sel1
    xy <- data %*% currentplot$w[,1:2]
    grouped <<- union(grouped,which(s$xmin<=xy[,1] & xy[,1]<=s$xmax &
                                        s$ymin<=xy[,2] & xy[,2]<=s$ymax))
    updateall()
  })
  
  observeEvent(input$selectpoints,{
    if(input$selectpoints!="none") grouped <<- union(grouped,subsets[[input$selectpoints]])
    updateall()
  })
  
  observeEvent(input$clearselection,{
    grouped <<- NULL
    updateall()
  })
  
  observeEvent(input$reverseselection,{
    grouped <<- setdiff(1:dim(data)[1],grouped)
    updateall()
  })
  
  
  observeEvent(input$saveselection,{
    selectioncount <<- selectioncount+1
    subsets[[sprintf("%d",selectioncount)]] <<- grouped
    updateSelectInput(session,"selectpoints",
                      label="add to current selection:",
                      choices=c("none",names(subsets)))
    grouped <<- NULL
    updateall()
  })

  observeEvent(input$removeall,{
    subsets <<- subsets0  
    constraints <<- constraints0
    constraintsubsets <<- constraintsubsets0
    updateSelectInput(session,"selectpoints",
                      label="add to current selection:",
                      choices=c("none",names(subsets)))
    grouped <<- NULL
    updateall()
  })
  
  observeEvent(input$deleteselection,{
    if(input$selectpoints!="none") {
      subsets[[input$selectpoints]] <<- NULL
      idx <- constraintsubsets==input$selectpoints
      constraintsubsets <<- constraintsubsets[!idx]
      constraints <<- constraints[!idx]
      updateSelectInput(session,"selectpoints",
                        label="add to current selection:",
                        choices=c("none",names(subsets)))
      grouped <<- NULL
    }
    updateall()
  })
  
  observeEvent(input$twod,{
    if(length(grouped)>0) {
      selectioncount <<- selectioncount+1
      subsetname <- sprintf("%d-2d",selectioncount)
      subsets[[subsetname]] <<- grouped
      updateSelectInput(session,"selectpoints",
                        label="add to current selection:",
                        choices=c("none",names(subsets)))
      w <- currentplot$w[,1:2]
      # Orthonormal basis for vectors
      w[,1] <- w[,1]/sqrt(sum(w[,1]^2))
      w[,2] <- w[,2]-w[,1]*sum(w[,1]*w[,2])
      w[,2] <- w[,2]/sqrt(sum(w[,2]^2))
      w <- w %*% svd(mycov(data[grouped,] %*% w))$u
    
      newconstraints <- list(linconstraint(data,grouped,w[,1]),
                             linconstraint(data,grouped,w[,2]),
                             quadconstraint(data,grouped,w[,1]),
                             quadconstraint(data,grouped,w[,2]))
      constraints <<- c(constraints,newconstraints)
      constraintsubsets <<- c(constraintsubsets,rep(subsetname,length(newconstraints)))
      grouped <<- NULL
    }
    updateall()
  })
  
  observeEvent(input$cluster,{
    if(length(grouped)>0) {
      selectioncount <<- selectioncount+1
      subsetname <- sprintf("%d-cluster",selectioncount)
      subsets[[subsetname]] <<- grouped
      updateSelectInput(session,"selectpoints",
                        label="add to current selection:",
                        choices=c("none",names(subsets)))
      w <- currentplot$w[,1:2]
      newconstraints <- clusterconstraint(data,grouped)
      constraints <<- c(constraints,newconstraints)
      constraintsubsets <<- c(constraintsubsets,rep(subsetname,length(newconstraints)))
      grouped <<- NULL
    }
    updateall()
  })
  
  observeEvent(input$recompute,{
    input$refresh
    
    if(length(constraints)>0) {
      aux <- getwhitenedR(data,constraints,tol1=10^input$tol1,tol2=10^input$tol2,maxtime=input$timeout)
      rdata <<- aux$rdata
      Y <<- aux$Y
    } else {
      rdata <<- matrix(rnorm(dim(data)[1]*dim(data)[2]),dim(data)[1],dim(data)[2])
      Y <<- data
    }
    updateall()
  })
  
  
  observeEvent(input$computepca,{
    cat("computing pca")
    currentplotpca <<- c(list(s="PCA"),dopca(Y))
    if(input$plot=="pca")
      currentplot <<- currentplotpca
    cat(" done\n")
    updateall()
  })
  
  observeEvent(input$computeica,{
    cat("computing ica")
    currentplotica <<- c(list(s="ICA"),dofastica(Y))
    currentplot <<- currentplotpca
    if(input$plot=="ica")
      currentplot <<- currentplotica
    cat(" done\n")
    updateall()
  })
  
  observeEvent(input$computeselection,{
    cat("computing selection")
    currentplotsel <<- c(list(s="SEL"),doselection(data,grouped))
    if(input$plot=="selection")
      currentplot <<- currentplotica
    cat(" done\n")
    updateall()
  })
  
  observeEvent(input$plot,{
    if(input$plot=="pca")
      currentplot <<- currentplotpca
    else if(input$plot=="ica")
      currentplot <<- currentplotica
    else
      currentplot <<- currentplotsel
    updateall()
  })
})
