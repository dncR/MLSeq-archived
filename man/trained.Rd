\name{trained-methods}
\docType{methods}
\alias{trained}
\alias{trained,MLSeq-method}
\title{Accessors for the 'trained' slot of an MLSeq object }
\usage{
  \S4method{trained}{MLSeq}(object)

}
\arguments{
  \item{object}{an \code{MLSeq} object}

}

\description{
Details about the training model information which is obtained \code{classify} function.

}
\details{
\code{trained} slot stores information about the training process such as optimum model parameters and resampling properties on the fitted classification model.
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

trained(svm)

}
