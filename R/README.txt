sideR

This is an attachment to the following manuscript:

Kai Puolamäki, Emilia Oikarinen, Bo Kang, Jefrey Lijffijt, Tijl De
Bie. Interactive Visual Data Exploration with Subjective Feedback: An
Information-Theoretic Approach. arXiv:1710.08167
[stat.ML]. https://arxiv.org/abs/1710.08167

Reference to sideR:

Kai Puolamäki. sideR - a tool for subjective and interactive visual
data exploration in R. Downloaded from
http://www.iki.fi/kaip/sider.html on 24 October 2017.

* *

You need to have the following software installed before running
sideR:
* R 3.4.0 or higher with the following packages which can be installed
  from CRAN (see https://www.r-bloggers.com/installing-r-packages/ for
  instructions): shiny, fastICA.
* A modern web browser such as Google Chrome.

We have tested sideR with a unix system (Mac OS X) but it should work
in any system that has R with the above mentioned packages installed.

NB: R 3.4.0 includes somes substantial speedups to matrix computations
which may affect the performance of the sideR. If you have an older
version of R installed we therefore strongly recommended you to update
to R 3.4.0 or newer!

* * *

Run sideR by executing:
Rscript --vanilla run.R
Then point your web browser to the link shown in the terminal. The
link should be of form http://127.0.0.1:xxxx where xxxx are random
four digits.

Run the runtime experiment by executing in directory runtime:
Rscript --vanilla runtime.R
Create the LaTeX table by executing in the same directory:
Rscript --vanilla report.R
Produce convergence figure by executing in directory toyexample:
Rscript --vanilla toyexample.R

* * *

How to get started by sideR:

sideR loads the BNC dataset by default. The BNC dataset contains
n=1335 documents which are represented by d=100 dimensional word
vectors. You can, e.g., mark a cluster you see on top right hand side
with your mouse, painting it red. If you make a mistake you can push
"clear current selection". Mark the selection as cluster constraint by
pressing "apply cluster constraint to current selection and
save". Wait a few seconds and hit "recompute background", after which
the background distribution is updated. When the computation is
finished you should see how the background distribution changed. You
can compute new PCA and ICA projections by pressing "compute pca
projection" and "compute ica projection", respectively. You can choose
the view (pca/ica/selection) by radio buttons.

You can view the classes and previously saved groupings by the
drop-down menu "add to current selection". The status line at the top
left corner shows concisely the Jaccard distance between the current
selection and clusters. E.g., from the top right cluster should
correspond quite well with the texts of the class
"Cconversation". Notice that from the axis of the scatterplot you can
see the projection vector in terms of the words (since we are here
dealing with text documents).

Use cases are described in the manuscript in more detail.

PS. Notice that due in Shiny the ordering of the update of global
variables is hard to enforce. This means that sometimes to update a
plot you should press "refresh all" button or swap between projections
(with radio buttons pca/ica/selection) a couple of time to force the
update.

* * *

Adding a data set:

You can easily add a data set. You should make a dataset a R data
frame with no missing values. The columns that are factors
(is.factor(x) true) are automatically converted to class variables
visible in the dropdown menu. The other columns are assumed to be real
valued attributes. You should save the data file in directory
sideR/data and the data file should end with .rds after which it will
automatically be read into the dataset dropdown menu when sideR
starts.

* * *

Data set descriptions and references, see the paper for more
discussion:

bnc.rds
Extract from the British National Corpus.
BNC 2007. The British National Corpus, version 3 (BNC XML
Edition). Distributed by Oxford University Computing Services on
behalf of the BNC Consortium. (2007). http://www.natcorp.ox.ac.uk/.

segment.rds
UCI Image Segmentation Data Set,
https://archive.ics.uci.edu/ml/datasets/Image+Segmentation
Lichman, M. (2013). UCI Machine Learning Repository
[http://archive.ics.uci.edu/ml]. Irvine, CA: University of California,
School of Information and Computer Science.

toy3.rds
toydata.rds
Artificial data sets created by us and described in the paper.

* * *

sideR

Copyright 2017 Kai Puolamaki <kai.puolamaki@iki.fi>
Copyright 2017 Finnish Institute of Occupational Health, Helsinki,
Finland

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
