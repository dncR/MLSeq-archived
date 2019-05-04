#' @include all_classes.R


######## selectedGenes #######
#' @title Accessors for the 'selectedGenes'.
#'
#' @description This slot stores the name of selected genes which are used in the classifier.
#' The trained model is stored in slot \code{trainedModel}. See \code{\link{trained}} for details.
#'
#' @docType methods
#' @name selectedGenes
#' @rdname selectedGenes
#'
#' @param object an \code{MLSeq} object.
#'
#' @seealso \code{\link{trained}}
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' selectedGenes(cart)
#'}
#'
#' @export
setGeneric("selectedGenes", function(object) standardGeneric("selectedGenes"))


######## input #######
#' @title Accessors for the 'inputObject' slot of an \code{MLSeq} object
#'
#' @description \code{MLSeq} package benefits from \code{DESeqDataSet} structure from bioconductor package \code{DESeq2} for storing gene
#' expression data in a comprehensive structure. This object is used as an input for classification task through \code{\link{classify}}.
#' The input is stored in \code{inputObject} slot of \code{MLSeq} object.
#'
#' @docType methods
#' @name input
#' @rdname input
#'
#' @param object an \code{MLSeq} object.
#'
#' @seealso \code{\link{classify}}, \code{\link[DESeq2]{DESeqDataSet}}
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' input(cart)
#'}
#'
#' @export
setGeneric("input", function(object) standardGeneric("input"))


######## preProcessing #######
#' @title Accessors for the 'preProcessing' slot of an \code{MLSeq} object
#'
#' @description \code{MLSeq} package benefits from \code{DESeqDataSet} structure from bioconductor package \code{DESeq2} for storing gene
#' expression data in a comprehensive structure. This object is used as an input for classification task through \code{\link{classify}}.
#' The input is stored in \code{inputObject} slot of \code{MLSeq} object.
#'
#' @docType methods
#' @name preProcessing
#' @rdname preProcessing
#'
#' @param object an \code{MLSeq} object.
#'
#' @seealso \code{\link{classify}}, \code{\link[DESeq2]{DESeqDataSet}}
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' preProcessing(cart)
#'}
#'
#' @export
setGeneric("preProcessing", function(object) standardGeneric("preProcessing"))

#' @param value a character string. Which preProcessing should be replaced with current one?
#' @rdname preProcessing
#' @export
setGeneric("preProcessing<-", function(object, value) standardGeneric("preProcessing<-"))



######## method #######
#' Accessors for the 'method'.
#'
#' This slot stores the name of selected model which is used in \code{classify} function.
#' The trained model is stored in slot \code{trainedModel}.
#' See \code{\link{trained}} for details.
#'
#' \code{method} slot stores the name of the classification method such as "svmRadial" for Radial-based Support Vector Machines, "rf" for Random Forests, "voomNSC" for
#' voom-based Nearest Shrunken Centroids, etc. For the complete list of available methods, see \code{\link{printAvailableMethods}} and \code{\link{availableMethods}}.
#'
#' @docType methods
#' @name method
#' @rdname method
#'
#' @param object an \code{MLSeq} object.
#'
#' @seealso \code{\link{trained}}
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' method(cart)
#'}
#'
#'@export
setGeneric("method", function(object) standardGeneric("method"))

#' @param value a character string. One of the available classification methods to replace with current method stored in MLSeq object.
#' @rdname method
#' @export
setGeneric("method<-", function(object, value) standardGeneric("method<-"))


######## transformation #######
#' Accessors for the 'transformation' slot.
#'
#' This slot stores the name of transformation method which is used while transforming the count data (e.g "vst", "rlog", etc.)
#'
#' @docType methods
#' @name transformation
#' @rdname transformation
#'
#' @param object an \code{MLSeq} or \code{MLSeqModelInfo} object.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' transformation(cart)
#'}
#'
#'@export
setGeneric("transformation", function(object) standardGeneric("transformation"))


######## normalization #######
#' Accessors for the 'normalization' slot.
#'
#' This slot stores the name of normalization method which is used while normalizing the count data such
#' as "deseq", "tmm" or "none"
#'
#' @docType methods
#' @name normalization
#' @rdname normalization
#'
#' @param object an \code{MLSeq} or \code{MLSeqModelInfo} object.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' normalization(cart)
#'}
#'
#'@export
setGeneric("normalization", function(object) standardGeneric("normalization"))

#' @param value a character string. One of the available normalization methods for voom-based classifiers.
#' @rdname normalization
#' @export
setGeneric("normalization<-", function(object, value) standardGeneric("normalization<-"))


######## confusionMat #######
#' Accessors for the 'confusionMat' slot.
#'
#' This slot stores the confusion matrix for the trained model using \code{classify} function.
#'
#' \code{confusionMat} slot includes information about cross-tabulation of observed and predicted classes
#' and corresponding statistics such as accuracy rate, sensitivity, specifity, etc. The returned object
#' is in \code{confusionMatrix} class of caret package. See \code{\link[caret]{confusionMatrix}} for details.
#'
#' @docType methods
#' @name confusionMat
#' @rdname confusionMat
#'
#' @param object an \code{MLSeq} or \code{MLSeqModelInfo} object.
#'
#' @seealso \code{\link[caret]{confusionMatrix}}
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' confusionMat(cart)
#'}
#'
#'@export
setGeneric("confusionMat", function(object) standardGeneric("confusionMat"))



######## trainedModel #######
#' Accessors for the 'trainedModel' slot.
#'
#' This slot stores the trained model. This object is returned from \code{train} function in caret package.
#' Any further request using caret functions is available for \code{trainedModel} since this object is in the
#' same class as the returned object from \code{train}. See \code{\link[caret:train]{train}} for details.
#'
#' @docType methods
#' @name trained
#' @rdname trained
#'
#' @param object an \code{MLSeq} or \code{MLSeqModelInfo} object.
#'
#' @seealso \code{\link[caret:train]{train.default}}, \code{\link{voom.train-class}}, \code{\link{discrete.train-class}}
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' trained(cart)
#'}
#'
#'@export
setGeneric("trained", function(object) standardGeneric("trained"))



######## ref #######
#' Accessors for the 'ref' slot.
#'
#' This slot stores the information about reference category. Confusion matrix and related statistics are calculated using
#' the user-defined reference category.
#'
#' @docType methods
#' @name ref
#' @rdname ref
#'
#' @param object an \code{MLSeq} or \code{MLSeqModelInfo} object.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' ref(cart)
#'}
#'
#'@export
setGeneric("ref", function(object) standardGeneric("ref"))

#' @param value a character string. Select reference category for class labels.
#' @rdname ref
#' @export
setGeneric("ref<-", function(object, value) standardGeneric("ref<-"))



######## control #######
#' Accessors for the 'control' slot.
#'
#' This slot stores the information about control parameters of selected classification model.
#'
#' @docType methods
#' @name control
#' @rdname control
#'
#' @param object an \code{MLSeq} or \code{MLSeqModelInfo} object.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' control(cart)
#'}
#'
#'@export
setGeneric("control", function(object) standardGeneric("control"))

#' @param value a list with elements for controlling trained model. It should be a list returned from one of 
#' \code{\link{discreteControl}}, \code{\link{voomControl}}, \code{\link[caret]{trainControl}} functions.
#' @rdname control
#' 
#' @seealso \code{\link{discreteControl}}, \code{\link{voomControl}}, \code{\link[caret]{trainControl}}
#' 
#' @export
setGeneric("control<-", function(object, value) standardGeneric("control<-"))



######## trainParameters #######
#' Accessors for the 'trainParameters' slot.
#'
#' This slot stores the transformation and normalization parameters from train set. These parameters are used to normalize and
#' transform test set using train set parameters.
#'
#' @docType methods
#' @name trainParameters
#' @rdname trainParameters
#'
#' @param object an \code{MLSeq} or \code{MLSeqModelInfo} object.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' trainParameters(cart)
#'}
#'
#'@export
setGeneric("trainParameters", function(object) standardGeneric("trainParameters"))



######## modelInfo #######
#' Accessors for the 'modelInfo' slot of an \code{MLSeq} object
#'
#' This slot stores all the information about classification model.
#'
#' @docType methods
#' @name modelInfo
#' @rdname modelInfo
#'
#' @param object an \code{MLSeq} object.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' modelInfo(cart)
#'}
#'
#'@export
setGeneric("modelInfo", function(object) standardGeneric("modelInfo"))



######## metaData #######
#' Accessors for the 'metaData' slot of an \code{MLSeq} object
#'
#' This slot stores metada information of \code{MLSeq} object.
#'
#' @docType methods
#' @name metaData
#' @rdname metaData
#'
#' @param object an \code{MLSeq} object.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' metaData(cart)
#'}
#'
#'@export
setGeneric("metaData", function(object) standardGeneric("metaData"))



######## isUpdated / isModified #######
#' Checks if MLSeq object is updated/modified or not.
#'
#' These functions are used to check whether the \code{MLSeq} object is modified and/or updated. It is possible to update
#' classification parameters of \code{MLSeq} object which is returned by \code{classify()} function.
#'
#' @return a logical.
#'
#' @docType methods
#' @name isUpdated
#' @rdname isUpdated
#'
#' @param object an \code{MLSeq} object.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (30% test, 70% train).
#' nTest <- ceiling(n*0.3)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition = class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                    colData = classtr, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' isUpdated(cart)
#' isModified(cart)
#'}
#'
#' @export
setGeneric("isUpdated", function(object) standardGeneric("isUpdated"))

#' @param value a logical. Change the state of update info.
#' @rdname isUpdated
#' @export
setGeneric("isUpdated<-", function(object, value) standardGeneric("isUpdated<-"))

#' @rdname isUpdated
#' @export
setGeneric("isModified", function(object) standardGeneric("isModified"))

#' @rdname isUpdated
#' @export
setGeneric("isModified<-", function(object, value) standardGeneric("isModified<-"))








