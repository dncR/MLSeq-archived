## Introduction
[![Build Status](http://bioconductor.org/shields/build/release/bioc/MLSeq.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/MLSeq/)
[![InBioc](http://bioconductor.org/shields/years-in-bioc/MLSeq.svg)](http://bioconductor.org/packages/devel/bioc/html/MLSeq.html#since)

<br>

MLSeq is an R/BIOCONDUCTOR package, which provides over 90 algorithms including support vector machines (SVM),random forest (RF), classification and regression trees (CART), Poisson and Negative Binomial Linear Discriminant Analysis (PLDA, NBLDA) and voom-based classifiers (voomDLDA, voomNSC, etc.) for the classification of sequencing data. MLSeq requires a count table as an input which contains the number of reads mapped to each transcript for each sample. This kind of count data can be obtained from RNA-Seq experiments, also from other sequencing experiments such as DNA or ChIP-sequencing. MLSeq includes both normalization (e.g deseq median ratio, trimmed mean of M values) and transformation (variance stabiliation transformation, regularized logarithmic transformation, etc.) techniques which can be performed through classification process. Although the main purpose of MLSeq is to classify samples using a count matrix from RNA-Sequencing data, some of the classifiers which are called sparse classifiers such as PLDA and voomNSC can be used to detect significant features. 

To install the MLSeq package in R:

```{r, eval = FALSE, message=FALSE, warning=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MLSeq")
```

If you use MLSeq package in your research, please cite it as below:

> Goksuluk D, Zararsiz G, Korkmaz S, Eldem V, Zararsiz GE, Ozcetin E, Ozturk A, Karaagaoglu AE. MLSeq: Machine
  learning interface for RNA-sequencing data. Computer Methods and Programs in Biomedicine. 2019, 175:223-231.


To get BibTeX entry for LaTeX users, type the following:

```{r, eval = FALSE}
citation("MLSeq")
```

<br>

Please contact us, if you have any questions or suggestions:

  gokmenzararsiz@hotmail.com <br>
  dincer.goksuluk@gmail.com <br>
  selcukorkmaz@gmail.com

## News:

#### Major changes in version 2.x.y

* Functions are reconstructed using S4 systems and new classes such as `MLSeq`, `MLSeqMetaData` and `MLSeqModelInfo`.
* New classifiers from [caret](https://cran.r-project.org/web/packages/caret/index.html) package are now available for MLSeq. These functions can be used for transformed continuous data using one of transformation techniques which are provided by MLSeq's classification algorithms.
* A complete list of available classifiers can be viewed using `availableMethods()` and `printAvailableMethods()`.
* New setter and getter functions are included.
* Predictions are now evaluated usin generic function `predict(...)`. The older function `predictClassify(...)` can also be used for predictions.
* For more details see package manuals.
