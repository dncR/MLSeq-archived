\name{classify}
\alias{classify}

\title{
Fitting Classification Models to Sequencing Data}
\description{
This function fits classification algorithms to sequencing data and measures model performances using various statistics
}
\usage{
classify(data, method = c("svm", "bagsvm", "randomforest", "cart"), normalize = c("deseq", "none", "tmm"), 
deseqTransform = c("vst", "voom"), cv = 5, rpt = 10, B = 100, ref=NULL, ...)
}

\arguments{
  \item{data}{
\code{DESeqDataSet} instance
}
  \item{method}{
There are four methods available to perform classification: \code{svm}: support vector machines using radial-based kernel function, \code{bagsvm}: support vector machines with bagging ensemble, \code{randomForest}: random forest algorithm, \code{cart}: classification and regression trees algorithm.}
  \item{normalize}{
Normalization of count data for classification. \code{none}: Normalization is not applied. Count data is used for classification. \code{deseq}: deseq normalization. \code{tmm}: Trimmed mean of \code{M} values. 
}
  \item{deseqTransform}{
Transformation method applied after normalization.\code{vst}: variance stabilizing transformation. \code{voom}: voom transformation.

}
  \item{cv}{
Number of cross-validation folds.}
  \item{rpt}{
Number of complete sets of folds for computation.
}
  \item{B}{
Number of bootstrap samples for \code{bagsvm} method.
}
  \item{ref}{
User defined reference class. Default is \code{NULL}.
}
  \item{\dots}{
Optional arguments for \code{train()} function from \code{caret} package.
}
}
\details{
In RNA-Seq studies, normalization is used to adjust between-sample differences for further analysis. In this package, "deseq" and "tmm" normalization methods are available. "deseq" estimates the size factors by dividing each sample by the geometric means of the transcript counts. "tmm" trims the lower and upper side of the data by log fold changes to minimize the log-fold changes between the samples and by absolute intensity. After normalization, it is useful to transform the data for classification. \code{MLSeq} package has "voom" and "vst" transformation methods. "voom" transformation applies a logarithmic transformation to normalized count data and computes gene weights using the mean-dispersion relationship. Second transformation method is the "vst" transformation and this approach uses an error modeling and the concept of variance stabilizing transformations to estimate the mean-dispersion relationship of data.

For model validation, k-fold cross-validation ("cv" option in \code{MLSeq} package) is a widely used technique. Using this technique, training data is randomly splitted into k non-overlapping and equally sized subsets. A classification model is trained on (k-1) subsets and tested in the remaining subsets. \code{MLSeq} package also has the repeat option as "rpt" to obtain more generalizable models. Giving a number of m repeats, cross validation concept is applied m times. 

For more details, see the vignette.
}
\value{
\item{model }{fitted classification model}
\item{method }{used classification method}
\item{normalization }{used normalization method}
\item{deseqTransform }{deseq transformation if \code{deseq} normalization is used}
\item{confusionMat }{cross-tabulation of observed and predicted classes and corresponding statistics}
\item{ref}{reference class}
}
\references{
Kuhn M. (2008). Building predictive models in R using the caret package. Journal of Statistical Software, (http://www.jstatsoft.org/v28/i05/).

Anders S. Huber W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11:R106

Witten DM. (2011). Classification and clustering of sequencing data using a poisson model. The Annals of Applied Statistics, 5(4), 2493:2518.

Charity WL. et al. (2014) Voom: precision weights unlock linear model analysis tools for RNA-Seq read counts, Genome Biology, 15:R29, doi:10.1186/gb-2014-15-2-r29

Witten D. et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical biomarker analysis in a collection of cervical tumours and matched controls. BMC Biology, 8:58

Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-Seq data. Genome Biology, 11:R25, doi:10.1186/gb-2010-11-3-r25
}
\author{
Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
}

\seealso{
\code{\link[MLSeq]{predictClassify}}
}
\examples{
data(cervical)

data = cervical[c(1:150),]  # a subset of cervical data with first 150 features.

class = data.frame(condition=factor(rep(c("N","T"),c(29,29))))# defining sample classes.

n = ncol(data)  # number of samples
p = nrow(data)  # number of features

nTest = ceiling(n*0.2)  # number of samples for test set (20% test, 80% train).
ind = sample(n,nTest,FALSE)

# train set
data.train = data[,-ind]
data.train = as.matrix(data.train + 1)
classtr = data.frame(condition=class[-ind,])

# train set in S4 class
data.trainS4 = DESeqDataSetFromMatrix(countData = data.train,
colData = classtr, formula(~ condition))
data.trainS4 = DESeq(data.trainS4, fitType="local")

# Support Vector Machines (SVM) Classification
svm = classify(data = data.trainS4, method = "svm", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref="T")
svm

# Random Forest (RF) Classification
rf = classify(data = data.trainS4, method = "randomforest", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref="T")
rf
}

\keyword{ RNA-seq classification }
