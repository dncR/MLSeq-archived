\name{ref-methods}
\docType{methods}
\alias{ref}
\alias{ref,MLSeq-method}
\title{ Accessors for the 'ref' slot of an MLSeq object  }
\usage{
  \S4method{ref}{MLSeq}(object)

}

\arguments{
  \item{object}{an \code{MLSeq} object}

}
\description{
The reference class category for the trained model using \code{classify} function.
}
\details{
Reference class category is important while calculating the statistical measures for the confusion matrix obtained from classification models.
}

\author{
Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
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

ref(svm)
}