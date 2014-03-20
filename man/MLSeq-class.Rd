\name{MLSeq-class}
\Rdversion{1.1}
\docType{class}
\alias{MLSeq-class}

\title{MLSeq object}
\description{
For classification, this is the main class for the \code{MLSeq} package.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ClassifySeq", ...)}.

This type of objects is created as a result of \code{classify} function of \code{MLSeq} package. It is then used in \code{predictClassify} function for predicting the class labels of new samples.
}
\section{Slots}{
  \describe{
    \item{\code{method}:}{stores the name of used classification method in the classification model}
    \item{\code{deseqTransform}:}{stores the name of used transformation method in the classification model}
    \item{\code{normalization}:}{stores the name of used normalization method in the classification model}
    \item{\code{confusionMat}:}{stores the information of classification performance results}
    \item{\code{trained}:}{stores the information about training process and model parameters that used in the corresponding model}
  	\item{ref}{stores user defined reference class}
  }
}

\author{
Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk}
\note{
An \code{MLSeq} class stores the results of \code{classify} function and offers further slots that are populated during the analysis. The slot \code{confusionMat} stores the information of classification performance results. These results contain the classification table and several statistical measures including accuracy rate, sensitivity, specifity, positive and negative predictive rates, etc. \code{method}, \code{normalization} and \code{deseqTransform} slots store the name of used classification method, normalization method and transformation method in the classification model respectively. Lastly, the slot \code{trained} stores the information about training process and model parameters that used in the corresponding model.
}

\examples{

# See the vignette

}