###### 1. Class: discrete.train #####
#' \code{discrete.train} object
#'
#' @description This object is the subclass for the \code{MLSeq.train} class. It contains trained model information for discrete
#' classifiers such as Poisson Linear Discriminant Analysis (PLDA) and Negative Binomial Linear Discriminant Analysis (NBLDA).
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{inputs}:}{a list with elements used as input for classification task.}
#'   \item{\code{control}:}{a list with control parameters for discrete classifiers, e.g. PLDA, PLDA2 and NBLDA.}
#'   \item{\code{crossValidatedModel}:}{a list. It stores the results for cross validation.}
#'   \item{\code{finalModel}:}{a list. This is the trained model with optimum parameters.}
#'   \item{\code{tuningResults}:}{a list. It stores the results for tuning parameter if selected classifier has
#'   one or more parameters to be optimized.}
#'   \item{\code{callInfo}:}{a list. call info for selected method.}
#' }
#'
#' @docType class
#' @name discrete.train-class
#' @rdname discrete.train-class
#'
#' @exportClass discrete.train
setClass("discrete.train",
         slots = c(inputs = "list",
                   control = "list",
                   crossValidatedModel = "list",
                   finalModel = "list",
                   tuningResults = "list",
                   callInfo = "list"))

setValidity("discrete.train", function(object){
  TRUE    ##### SLOTLARIN VALIDITY'SI YAPILACAK.
})

###### 2. Class: voom.train #####
#' \code{voom.train} object
#'
#' @description This object is the subclass for the \code{MLSeq.train} class. It contains trained model information for voom based
#' classifiers, i.e. "voomDLDA", "voomDQDA" and "voomNSC".
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{weigtedStats}:}{a list with elements of weighted statistics which are used for training the model. Weights are calculated from voom
#'   transformation. }
#'   \item{\code{foldInfo}:}{a list containing information on cross-validated folds.}
#'   \item{\code{control}:}{a list with control parameters for voom based classifiers.}
#'   \item{\code{tuningResults}:}{a list. It stores the cross-validation results for tuning parameter(s).}
#'   \item{\code{finalModel}:}{a list. It stores results for trained model with optimum parameters.}
#'   \item{\code{callInfo}:}{a list. call info for related function.}
#' }
#'
#' @docType class
#' @name voom.train-class
#' @rdname voom.train-class
#'
#' @exportClass voom.train
setClass("voom.train",
         slots = c(weigtedStats = "list",
                   foldInfo = "list",
                   control = "list",
                   tuningResults = "list",
                   finalModel = "list",
                   callInfo = "list"))

# setValidity("voom.train", function(object){
#   TRUE    ##### SLOTLARIN VALIDITY'SI YAPILACAK.
# })


# setOldClass(c("confusionMatrix","train", "DESeqDataSet"))
setOldClass(c("confusionMatrix", "train"))
setClassUnion("MLSeq.train", c("train", "voom.train", "discrete.train"))    ## Details on this class is given under "MLSeqModelInfo" documentation.
setClassUnion("confMat", c("confusionMatrix"))

###### 3. Class: MLSeqModelInfo #####
#' \code{MLSeqModelInfo} object
#'
#' @description For classification, this is the subclass for the \code{MLSeq} class. This object contains all the information about classification model.
#'
#' @details Objects can be created by calls of the form \code{MLSeqModelInfo(...)}. This type
#' of objects is created as a result of \code{classify} function of \code{MLSeq} package.
#' It is then used in \code{predictClassify} function for predicting the class labels of new samples.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{method, transformation, normalization}:}{these slots store the classification method, transformation technique and
#'   normalization method respectively. See notes for details. }
#'   \item{\code{preProcessing}:}{See \code{\link{classify}} for details.}
#'   \item{\code{ref}:}{a character string indicating the reference category for cases (diseased subject, tumor sample, etc.)}
#'   \item{\code{control}:}{a list with controlling parameters for classification task.}
#'   \item{\code{confusionMat}:}{confusion table and accuracy measures for the predictions.}
#'   \item{\code{trainedModel}:}{an object of \code{MLSeq.train} class. It contains the trained model. See notes for details.}
#'   \item{\code{trainParameters}:}{a list with training parameters from final model. These parameters are used for test set before predicting class labels.}
#'   \item{\code{call}:}{a call object for classification task.}
#' }
#'
#' @note \code{method, transformation, normalization} slots give the information on classifier, transformation and normalization techniques.
#' Since all possible pairs of transformation and normalization are not available in practice, we specify appropriate transformations and
#' normalization techniques with \code{preProcessing} argument in \code{\link{classify}} function. Finally, the information on normalization and transformation
#' is extracted from preProcessing argument.
#'
#' \code{MLSeq.train} is a union class of \code{train} from caret package, \code{voom.train} and \code{discrete.train} from MLSeq package. See related class
#' manuals for details.
#'
#' @docType class
#' @name MLSeqModelInfo-class
#' @rdname MLSeqModelInfo-class
#'
#' @seealso \code{\link[caret]{train}}, \code{\link{voom.train-class}}, \code{\link{discrete.train-class}}
#'
#' @exportClass MLSeqModelInfo
setClass("MLSeqModelInfo",
         slots = c(method = "character",
                   transformation = "character",
                   normalization = "character",
                   preProcessing = "character",
                   ref = "character",
                   control = "list",
                   confusionMat = "confMat",
                   trainedModel = "MLSeq.train",
                   trainParameters = "list",
                   call = "list"),
         prototype = prototype(confusionMat = structure(list(), class = "confMat"),
                               trainedModel = structure(list(), class = "MLSeq.train"),
                               trainParameters = structure(list(), class = "list"),
                               call = structure(list(), class = "list")))


# setValidity("MLSeqModelInfo", function(object){
#
#   if (!(method(object)  %in% availableMethods(method = NULL))){
#     return("Error: 'method' slot must be in one of the available methods in MLSeq. See \"availableMethods()\" for details.")
#   }
#
#   # preProcessing <- c("deseq-vst", "deseq-rlog", "deseq-logcpm", "tmm-logcpm", "logcpm")
#   if (!(normalization(object)  %in% c("deseq", "none", "tmm"))){
#     return("Error: 'normalization' slot must be in one of the following: \"deseq\", \"none\", \"tmm\" ")
#   }
#
#   if (!(transformation(object)  %in% c("vst", "rlog", "logcpm", "none"))){
#     return("Error: 'transformation' slot must be in one of the following: \"vst\", \"logcpm\", \"rlog\", \"none\" ")
#   }
#
#   if (!is.character(ref(object))){
#     return("Error: 'ref' slot must be a character ")
#   }
#
#   TRUE
# })



###### 4. Class: MLSeqMetaData #####
#' \code{MLSeqMetaData} object
#'
#' @description This object is a subclass for the \code{MLSeq} class. It contains metadata information, i.e. information on modified and/or
#' updated elements, raw data etc..
#'
#' @details Objects can be created by calls of the form \code{new("MLSeqMetaData", ...)}. This type
#' of objects is created as a result of \code{classify} function of \code{MLSeq} package.
#' It is then used in \code{update} function for updating the object in given object.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{updated, modified}:}{a logical. See notes for details.}
#'   \item{\code{modified.elements}:}{a list containing the modified elements in \code{MLSeq} obejct.}
#'   \item{\code{rawData.DESeqDataSet}:}{raw data which is used for classification.}
#'   \item{\code{classLabel}:}{a character string indicating the name of class variable.}
#' }
#'
#' @note The function \code{\link{update}} is used to re-run classification task with modified elements in \code{MLSeq} object. This function is
#' useful when one wish to perform classification task with modified options without running \code{classify} function from the beginning.
#' \code{MLSeqMetaData} object is used to store information on updated and/or modified elements in MLSeq object.
#'
#' If an \code{MLSeq} object is modified, i.e. one or more elements in MLSeq object is replaced using related setter functions such as
#' \code{\link{method}}, \code{\link{ref}} etc., the slot \code{modified} becomes TRUE. Similarly, the slot \code{updated} stores the
#' information that the MLSeq object is updated (or classification task is re-runned) or not. If updated slot is FALSE and modified slot is TRUE, one
#' should run \code{\link{update}} to obtain the classification results by considering the modified elements.
#'
#' @docType class
#' @name MLSeqMetaData-class
#' @rdname MLSeqMetaData-class
#'
#' @seealso \code{\link{update}}, \code{\link{isUpdated}}, \code{\link{isModified}}
#'
#' @exportClass MLSeqMetaData
setClass("MLSeqMetaData",
         slots = c(updated = "logical",
                   modified = "logical",
                   modified.elements = "list",
                   rawData.DESeqDataSet = "DESeqDataSet",
                   classLabel = "character"),
         prototype = prototype(updated = FALSE,
                               modified = FALSE,
                               modified.elements = structure(list(modifiedElements = NULL), class = "list"),
                               classLabel = NA_character_,
                               rawData.DESeqDataSet = structure(list(), class = "DESeqDataSet")))
#
# setValidity("MLSeqMetaData", function(object){
#   TRUE    ##### SLOTLARIN VALIDITY'SI YAPILACAK.
# })

#
#

###### 5. Class: MLSeq #####
#' \code{MLSeq} object
#'
#' @description For classification, this is the main class for the \code{MLSeq} package. It contains all the information including trained model,
#' selected genes, cross-validation results, etc.
#'
#' @details Objects can be created by calls of the form \code{new("MLSeq", ...)}. This type
#' of objects is created as a result of \code{classify} function of \code{MLSeq} package.
#' It is then used in \code{\link{predict}} or \code{\link{predictClassify}} function for predicting the class labels of new samples.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{inputObject}:}{stores the data in \code{\link[DESeq2]{DESeqDataSet}} object.}
#'   \item{\code{modelInfo}:}{stores all the information about classification model. The object is from subclass \code{MLSeqModelInfo}. See \code{\link{MLSeqModelInfo-class}} for details.}
#'   \item{\code{metaData}:}{metadata for MLSeq object. The object is from subclass \code{MLSeqMetaData}. See \code{\link{MLSeqMetaData-class}} for details.}
#' }
#'
#' @note An \code{MLSeq} class stores the results of \code{classify} function and offers further slots that are populated
#' during the analysis. The slot \code{inputObject} stores the raw and transformed data throughout the classification. The slot
#' \code{modelInfo} stores all the information about classification model. These results may contain the classification table
#' and performance measures such as accuracy rate, sensitivity, specifity, positive and negative predictive values, etc. It also
#' contains information on classification method, normalization and transformation used in the classification model.
#' Lastly, the slot \code{metaData} stores the information about modified or updated slots in MLSeq object.
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Bernd Klaus, Ahmet Ozturk and Ahmet Ergun Karaagaoglu
#'
#' @docType class
#' @name MLSeq-class
#' @rdname MLSeq-class
#'
#' @seealso \code{\link{MLSeqModelInfo-class}}, \code{\link{MLSeqMetaData-class}}
#'
#' @exportClass MLSeq
setClass("MLSeq",    ## MLSEQ CLASS'i yeni yapiya gore tekrar getter fonksiyonlari ve slotlari duzenlenecek.
         slots = c(inputObject = "list",
                   modelInfo = "MLSeqModelInfo",
                   metaData = "MLSeqMetaData"),
         prototype = prototype(inputObject = structure(list(), class = "list"),
                               metaData = structure(list(), class = "MLSeqMetaData")))
#
# setValidity("MLSeq", function(object){
#   if ((class(input(object)) != "list")){
#     return("Error: 'inputObject' is not a \"list\" object.")
#   }
#
#   if ((class(modelInfo(object)) != "MLSeqModelInfo")){
#     return("Error: 'modelInfo' is not a \"MLSeqModelInfo\" object.")
#   }
#
#   # if ((class(metaData(object)) != "environment")){    ## metada icin bir getter yazilinca buradaki degerde guncellenecek. modelInfo(...)
#   #   return("Error: 'metaData' is not an \"environment\" object.")
#   # }
#
#   TRUE
# })
#





