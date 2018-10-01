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

shinyUI(fluidPage(
    titlePanel("sideR"),
    sidebarLayout(
        sidebarPanel(
          htmlOutput("datadescr"),
          selectInput("selectpoints",
                      label="add to current selection:",
                      choices=c("none",names(subsets))),
          actionButton("deleteselection",label="delete selection chosen obove"),
          actionButton("removeall",label="delete all saved selections"),
          actionButton("clearselection",label="clear current selection"),
          actionButton("reverseselection",label="reverse current selection"),
          actionButton("saveselection",label="save current selection"),
          actionButton("twod",label="apply 2d constraint to current selection and save"),
          actionButton("cluster",label="apply cluster constraint to current selection and save"),
          actionButton("recompute",label="recompute background"),
          radioButtons("plot",label="plot:",choices=c("pca","ica","selection")),
          actionButton("computepca",label="compute pca projection"),
          actionButton("computeica",label="compute ica projection"),
          actionButton("computeselection",label="compute selection projection"),
          actionButton("refresh",label="refresh all"),
          selectInput("dataset",label="dataset",choices=datafile),
          sliderInput("tol1","log10 of lambda tolerance",min=-5,max=1,value=-2,step=0.1),
          sliderInput("tol2","log10 of sigma tolerance",min=-5,max=1,value=-2,step=0.1),
          sliderInput("timeout","timeout (s)",min=1,max=300,value=10,step=1),
          HTML("<p>Kai Puolam&auml;ki &lt;kai.puolamaki@iki.fi&gt. Published under the MIT license. See the README.txt file and <a href=\"https://arxiv.org/abs/1710.08167\">arXiv:1710.08167</a> [stat.ML] for details.</p>")
        ),
        mainPanel(
            plotOutput("plot",brush=brushOpts(id="sel1"),click="plot_click"),
            plotOutput("selectionplot")
        )
    )
))

