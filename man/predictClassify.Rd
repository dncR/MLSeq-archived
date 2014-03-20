\name{predictClassify}
\alias{predictClassify}
\title{
Extract Predictions From \code{classify()} objects
}
\description{
This function predicts the class labels of test data for a given model.
}
\usage{
predictClassify(model, test.data)
}
\arguments{
  \item{model}{
a model of \code{MLSeq} class 
}
  \item{test.data}{
a \code{DESeqDataSet} instance of new observations.
}
}
\value{
\item{predicted }{a vector of predicted classes of test data. See details.}
}
\details{
\code{predictClassify} function gives a vector of predicted classes of data set. This vector is in \code{factor} class. 
}
\references{
Kuhn M. (2008). Building predictive models in R using the caret package. Journal of Statistical Software, (http://www.jstatsoft.org/v28/i05/).

Anders S. Huber W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11:R106

Witten DM. (2011). Classification and clustering of sequencing data using a poisson model. The Annals of Applied Statistics, 5(4), 2493:2518.

Charity WL. et al. (2014) Voom: precision weights unlock linear model analysis tools for RNA-seq read counts, Genome Biology, 15:R29, doi:10.1186/gb-2014-15-2-r29

Witten D. et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical biomarker analysis in a collection of cervical tumours and matched controls. BMC Biology, 8:58

Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-Seq data. Genome Biology, 11:R25, doi:10.1186/gb-2010-11-3-r25
}
\author{
Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
}

\seealso{
\code{\link[MLSeq]{classify}}
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
data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
colData = classtr, formula(~ condition))
data.trainS4 <- DESeq(data.trainS4, fitType="local")

# test set
data.test = data[,ind]
data.test = as.matrix(data.test + 1)
classts = data.frame(condition=class[ind,])

# test set in S4 
data.testS4 = DESeqDataSetFromMatrix(countData = data.test,
colData = classts, formula(~ condition))
data.testS4 = DESeq(data.testS4, fitType="local")

## Number of repeats (rpt) might change model accuracies ##

# Support Vector Machines (SVM) Classification
svm = classify(data = data.trainS4, method = "svm", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref="T")
svm

# Random Forest (RF) Classification
rf = classify(data = data.trainS4, method = "randomforest", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref="T")
rf

# predicted classes of test samples for SVM method
pred.svm = predictClassify(svm, data.testS4)
pred.svm

# predicted classes of test samples for RF method
pred.rf = predictClassify(rf, data.testS4)
pred.rf
}

