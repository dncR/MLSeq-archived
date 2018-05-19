#' @include predict.R
#' @include classify.R


######  1. Setter/Getter Functions for S4 Object slots  ##############

#### MODIFIED ELEMENTS için bir getter fonksiyon yazılacak.


#' @rdname selectedGenes
#' @export
setMethod("selectedGenes", signature(object = "MLSeq"), function(object){
  if (method(object) %in% c("PLDA", "PLDA2")){
    selectedGenesInfo <- object@modelInfo@trainedModel@finalModel$SelectedGenes
    if (is.null(selectedGenesInfo$selectedGenesIdx)){
      tmp <- metaData(object)@rawData.DESeqDataSet
      rownames(counts(tmp))
    } else {
      selectedGenesInfo$selectedGenesNames
    }
  } else if (method(object) == "voomNSC"){
    tmp <- unlist(object@modelInfo@trainedModel@finalModel$model$SelectedGenes)
    names(tmp) <- NULL
    tmp
  }
})


#' @rdname input
#' @export
setMethod("input", signature(object = "MLSeq"), function(object) object@inputObject)

#### inputObject ICIN BIR SETTER YAPILACAK MI? BU KONUYU TARTISALIM. INPUT OBJECT DEGISTIGI ANDA BUTUN
#### ANALIZ SONUCLARININ UPDATE EDILMESI GEREKIYOR. BOYLE BIR SETTER VERMEK YERINE KISININ INPUT VERILERINI DEGISTIRIP
#### BUTUN ANALIZI YENIDEN YAPMASINI BEKLEMEK DAHA MI MANTIKLI OLUR



#' @rdname preProcessing
#' @export
setMethod("preProcessing", signature(object = "MLSeq"), function(object){object@modelInfo@preProcessing})

#' @rdname preProcessing
#' @exportMethod "preProcessing<-"
setReplaceMethod("preProcessing", signature(object = "MLSeq", value = "character"),
                function(object, value){
                  # Change current method with a given value.
                  # object <- updateObject(object)
                  # if (class(trained(object)) == "train"){
                  if ((method(object) %in% c("voomNSC", "voomDLDA", "voomDQDA", "PLDA", "PLDA2", "NBLDA"))){
                    cat(" WARNING: 'preProcessing' slot can not be changed due to one of followings:", "\n")
                    cat("  (1) It is an option to be used in \"caret-based\" classifiers.", "\n")
                    cat("  (2) 'method' should be set before 'preProcessing'.", "\n")
                    cat("  (3) 'normalization' should be used for \"voom-based\" and \"discrete\" classifiers.", "\n\n")
                  } else {
                    object@modelInfo@preProcessing <- value
                    if (value == "deseq-vst"){
                      object@modelInfo@normalization <- "deseq"
                      object@modelInfo@transformation <- "vst"
                    } else  if (value == "deseq-rlog"){
                      object@modelInfo@normalization <- "deseq"
                      object@modelInfo@transformation <- "rlog"
                    } else  if (value == "deseq-logcpm"){
                      object@modelInfo@normalization <- "deseq"
                      object@modelInfo@transformation <- "logcpm"
                    } else  if (value == "tmm-logcpm"){
                      object@modelInfo@normalization <- "tmm"
                      object@modelInfo@transformation <- "logcpm"
                    } else  if (value == "logcpm"){
                      object@modelInfo@normalization <- "none"
                      object@modelInfo@transformation <- "logcpm"
                    } else {
                      stop(warning("Unknown method for 'preProcessing'."))
                    }

                    object@metaData@modified <- TRUE
                    object@metaData@updated <- FALSE

                    ## add modified element to metadata
                    temp <- object@metaData@modified.elements$modifiedElements
                    if (!("preProcessing" %in% temp)){
                      temp <- c(temp, "preProcessing")
                    }
                    object@metaData@modified.elements$modifiedElements <- temp
                  }
                  # else {
                  #   stop(warning("'preProcessing' is used only for 'caret-based' classifiers."))
                  # }

                  # validObject(object)
                  return(object)
                })



#' @rdname method
#' @export
setMethod("method", signature(object = "MLSeq"), function(object) object@modelInfo@method)

#' @rdname method
#' @aliases method
#' @export
setMethod("method", signature(object = "MLSeqModelInfo"), function(object) object@method)

#' @rdname method
#' @exportMethod "method<-"
setReplaceMethod("method", signature(object = "MLSeq", value = "character"),
  function(object, value){
    # Change current method with a given value.
    # object <- updateObject(object)
    if (value %in% availableMethods()){
      object@modelInfo@method <- value

      object@metaData@modified <- TRUE
      object@metaData@updated <- FALSE

      ## add modified element to metadata
      temp <- object@metaData@modified.elements$modifiedElements
      if (!("method" %in% temp)){
        temp <- c(temp, "method")
      }

      object@metaData@modified.elements$modifiedElements <- temp
      validObject(object)
    } else {
      stop(warning("Selected method is not available in 'MLSeq'. Run 'availableMethods()' for a list of available methods in MLSeq."))
    }
    return(object)
})





#' @rdname transformation
#' @export
setMethod("transformation", signature(object = "MLSeq"), function(object) object@modelInfo@transformation)

#' @rdname transformation
#' @aliases transformation
#' @export
setMethod("transformation", signature(object = "MLSeqModelInfo"), function(object) object@transformation)






#' @rdname normalization
#' @export
setMethod("normalization", signature(object = "MLSeq"), function(object) object@modelInfo@normalization)

#' @rdname normalization
#' @aliases normalization
#' @export
setMethod("normalization", signature(object = "MLSeqModelInfo"), function(object) object@normalization)

#' @rdname normalization
#' @exportMethod "normalization<-"
setReplaceMethod("normalization", signature(object = "MLSeq", value = "character"),
                function(object, value){
                  # Change current method with a given value.
                  # object <- updateObject(object)

                  # if (class(trained(object)) == "train"){
                  if (!(method(object) %in% c("voomNSC", "voomDLDA", "voomDQDA", "PLDA", "PLDA2", "NBLDA"))){
                    cat(" WARNING: 'normalization' slot can not be changed due to one of followings:", "\n")
                    cat("  (1) It is an option to be used in \"voom-based\" and \"discrete\" classifiers.", "\n")
                    cat("  (2) 'method' should be set before 'normalization'.", "\n")
                    cat("  (3) 'preProcessing' should be used for \"caret-based\" classifiers.", "\n\n")
                  } else {
                    if (value %in% c("deseq", "TMM", "none")){
                      object@modelInfo@normalization <- value

                      object@metaData@modified <- TRUE
                      object@metaData@updated <- FALSE

                      ## add modified element to metadata
                      temp <- object@metaData@modified.elements$modifiedElements
                      if (!("normalization" %in% temp)){
                        temp <- c(temp, "normalization")
                      }

                      object@metaData@modified.elements$modifiedElements <- temp
                      validObject(object)
                    } else {
                      stop(warning("Selected normalization method is not available in 'MLSeq'. See classify(...) manual for available options."))
                    }
                  }

                  return(object)
                })




#' @rdname confusionMat
#' @export
setMethod("confusionMat", signature(object = "MLSeq"), function(object){
  object@modelInfo@confusionMat
})

#' @rdname confusionMat
#' @export
setMethod("confusionMat", signature(object = "MLSeqModelInfo"), function(object){
  object@confusionMat
})





#' @rdname trained
#' @export
setMethod("trained", signature(object = "MLSeq"), function(object) object@modelInfo@trainedModel)

#' @rdname trained
#' @aliases trained
#' @export
setMethod("trained", signature(object = "MLSeqModelInfo"), function(object) object@trainedModel)





#' @rdname ref
#' @export
setMethod("ref", signature(object = "MLSeq"), function(object) object@modelInfo@ref)

#' @rdname ref
#' @aliases ref
#' @export
setMethod("ref", signature(object = "MLSeqModelInfo"), function(object) object@ref)

#' @rdname ref
#' @exportMethod "ref<-"
setReplaceMethod("ref", signature(object = "MLSeq", value = "character"),
                function(object, value){
                  # Change current method with a given value.
                  # object <- updateObject(object)
                  object@modelInfo@ref <- value

                  object@metaData@modified <- TRUE
                  object@metaData@updated <- FALSE

                  ## add modified element to metadata
                  temp <- object@metaData@modified.elements$modifiedElements
                  if (!("ref" %in% temp)){
                    temp <- c(temp, "ref")
                  }

                  object@metaData@modified.elements$modifiedElements <- temp
                  validObject(object)
                  return(object)
                })




#' @rdname control
#' @export
setMethod("control", signature(object = "MLSeq"), function(object){
  if (class(trained(object)) == "voom.train"){
    object@modelInfo@trainedModel@control
  } else if (class(trained(object)) == "train"){
    object@modelInfo@trainedModel$control
  } else if (class(trained(object)) == "discrete.train") {
    object@modelInfo@control
  }
})

#' @rdname control
#' @exportMethod "control<-"
setReplaceMethod("control", signature(object = "MLSeq", value = "list"),
                 function(object, value){
                   # Change current method with a given value.
                   # object <- updateObject(object)

                   if (class(trained(object)) == "voom.train"){
                     object@modelInfo@trainedModel@control <- value
                   } else if (class(trained(object)) == "train"){
                     object@modelInfo@trainedModel$control <- value
                   } else if (class(trained(object)) == "discrete.train"){
                     object@modelInfo@control <- value
                   }


                   object@metaData@modified <- TRUE
                   object@metaData@updated <- FALSE

                   ## add modified element to metadata
                   temp <- object@metaData@modified.elements$modifiedElements
                   if (!("control" %in% temp)){
                     temp <- c(temp, "control")
                   }

                   object@metaData@modified.elements$modifiedElements <- temp
                   validObject(object)
                   return(object)
                 })




#' @rdname trainParameters
#' @export
setMethod("trainParameters", signature(object = "MLSeq"), function(object) object@modelInfo@trainParameters)

#' @rdname trainParameters
#' @aliases trainParameters
#' @export
setMethod("trainParameters", signature(object = "MLSeqModelInfo"), function(object) object@trainParameters)

##### trainParameters icin bir setter hazirlanacak.
##### BURADAKI SETTERLAR EN SON DUZENLEMELERDEN SONRA YAPILACAK.




#' @rdname modelInfo
#' @export
setMethod("modelInfo", signature(object = "MLSeq"), function(object) object@modelInfo)




#' @rdname metaData
#' @export
setMethod("metaData", signature(object = "MLSeq"), function(object) object@metaData)   ## METADATA environment oldugu icin bu kisim incelenecek.




#' @rdname isUpdated
#' @export
setMethod("isUpdated", signature(object = "MLSeq"), function(object){
  object@metaData@updated
})

#' @rdname isUpdated
#' @export
setReplaceMethod("isUpdated", signature(object = "MLSeq", value = "logical"),
                 function(object, value){
                   # Change current method with a given value.
                   # object <- updateObject(object)
                   object@metaData@updated <- value
                   validObject(object)

                   return(object)
                 })




#' @rdname isUpdated
#' @export
setMethod("isModified", signature(object = "MLSeq"), function(object){
  object@metaData@modified
})

#' @rdname isUpdated
#' @export
setReplaceMethod("isModified", signature(object = "MLSeq", value = "logical"),
                 function(object, value){
                   # Change current method with a given value.
                   # object <- updateObject(object)
                   object@metaData@modified <- value
                   validObject(object)

                   return(object)
                 })



######  2. Show/Print method for S4 objects   ##########
####### show ########
#' @title Show method for MLSeq objects
#'
#' @description Prints out the information from the trained model using \code{classify} function.
#'
#' @docType methods
#' @name show
#' @rdname show
#'
#' @param object an \code{MLSeq} object returned from \code{classify} function.
#'
#' @seealso \code{\link{classify}}
show.MLSeq <- function(object){
  # object <- x ## use arguement object in the function if "show" method is specified.
  if (dim(confusionMat(object)$table)[1] == 2){
    cat("\n", sep = " ")
    cat("  An object of class ", "\"", class(object), "\"","\n", sep = "")
    cat("  Model Description: ", modelDescription(method(object)), " (", method(object), ")", "\n\n", sep = "")
    if (isModified(object) & !isUpdated(object)){
      cat("  NOTE: MLSeq object is modified but not updated.", "\n")
      cat("        Update 'MLSeq' object to get true classification accuracies.", "\n\n")
    }
    cat("            Method  : ", method(object), "\n\n")
    cat("       Accuracy(%)  : ", round(confusionMat(object)$overall[1],4)*100, "\n")
    cat("    Sensitivity(%)  : ", round(confusionMat(object)$byClass[1],4)*100, "\n")
    cat("    Specificity(%)  : ", round(confusionMat(object)$byClass[2],4)*100, "\n\n")
    cat("  Reference Class   : ", ref(object), "\n\n")
  } else {
    cat("\n", sep = " ")
    cat("  An object of class ", class(object), "\n\n", sep = " ")
    cat("            Method  : ", method(object), "\n\n")
    cat("       Accuracy(%)  : ", round(confusionMat(object)$overall[1],4)*100, "\n")
    cat("    Sensitivity(%)  : ", round(confusionMat(object)$byClass[1,1],4)*100, "\n")
    cat("    Specificity(%)  : ", round(confusionMat(object)$byClass[1,2],4)*100, "\n\n")
    cat("  Reference Class   : ", ref(object), "\n\n")
    invisible(NULL)
  }
}

#' @rdname show
setMethod(f = "show", signature = signature(object = "MLSeq"), definition = show.MLSeq)


## Belirli bir class'a gore yazilan fonksiyonlar .<classname> formunda yazildiktan sonra bu fonksiyonun
##  setMethod ile nasil davranacagi belirleniyor. export edilme zorunlulugu yok.
####### TODO ######
## MLSeqModelInfo için düzenleme yapılacak. ###
# show.MLSeqModelInfo <- function(object){
#   # object <- x
#   cat("  An object of class ", "\"", class(object), "\"","\n", sep = "")
#   cat("  Model Description: ", modelDescription(method(object)),"\n\n")
#   cat("            Method  : ", method(object), "\n\n")
#   cat("       Accuracy(%)  : ", round(confusionMat(object)$overall[1],4)*100, "\n")
# }

show.MLSeqModelInfo <- function(object){
  return(NULL)
}

#' @rdname show
setMethod(f = "show", signature = signature(object = "MLSeqModelInfo"), definition = show.MLSeqModelInfo)



## Belirli bir class'a gore yazilan fonksiyonlar .<classname> formunda yazildiktan sonra bu fonksiyonun
##  setMethod ile nasil davranacagi belirleniyor. export edilme zorunlulugu yok.
show.MLSeqMetaData <- function(object){
  # object <- x
  cat("class:  ", class(object), ifelse(isS4(object), ", in S4 class", ""), sep = "", "\n")
  cat("Updated:  ", ifelse(object@updated, "YES", "NO"), sep = "", "\n")
  cat("Modified:  ", ifelse(object@modified, "YES", "NO"), sep = "", "\n")

  elements <- if (length(object@modified.elements$modifiedElements) > 3){
    paste(paste(object@modified.elements$modifiedElements[1:3], sep = "", collapse = ", "), ", ...", sep = "")
  } else {
    paste(object@modified.elements$modifiedElements, sep = "", collapse = ", ")
  }
  nModifiedElements <- length(object@modified.elements$modifiedElements)
  modified.text <- paste(" (", nModifiedElements, ")", sep = "")

  cat("Modified Elements", ifelse(nModifiedElements != 0, paste(modified.text, ":  ", sep = ""), ":  "), ifelse(is.null(object@modified.elements$modifiedElements), "(NULL)", elements), sep = "", "\n")
  cat("Initial Data:  A ", class(object@rawData.DESeqDataSet), " object", sep = "", "\n\n")
}

#' @rdname show
setMethod(f = "show", signature = signature(object = "MLSeqMetaData"), definition = show.MLSeqMetaData)



## voom ile egitilen objenin "voom.train" ekranda goruntulenme yontemi
show.voom.train <- function(object){
  # object <- x
  n <- ceiling(object@weigtedStats$n)
  p <- ceiling(object@weigtedStats$p)
  nc <- ceiling(object@weigtedStats$nclass)
  ndigits <- length(unlist(strsplit(as.character(max(n, p, nc)), "")))
  class.names <- object@callInfo$class.names
  cNames <- paste(sapply(class.names, function(x)paste("\'", x, "\'", sep = "")), collapse = ", ")

  # PART 1. Model description and summary of sample/feature sizes. Class labels info.
  cat("\n", modelDescription(object@callInfo$method), " (", object@callInfo$method, ")", sep = "", "\n\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " samples", sep = ""), n), "\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " predictors", sep = ""), p), "\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " classes: ", cNames, "  (Reference category: '", object@callInfo$ref, "')", sep = ""), nc), "\n\n")

  normalizationInfo <- if (object@callInfo$normalize == "none"){
    "Normalization is NOT performed."
  } else if (object@callInfo$normalize == "deseq"){
    "DESeq median ratio."
  } else {
    "Trimmed-mean of M values."
  }

  # PART 2. Normalization and Resampling (Cross validation) info
  cat("Normalization: ", normalizationInfo, "\n", sep = "")
  cat("Resampling: Cross-Validated (", object@control$number, " fold, repeated ", object@control$repeats, " times)", sep = "", "\n")

  foldSampleSize <- unlist(lapply(object@foldInfo$foldIndex$indexIn, length))
  foldSampleSizeText <- if (length(foldSampleSize) > 5){
    paste(c(foldSampleSize[1:5], "..."), collapse = ", ", sep = "")
  } else {
    paste(c(foldSampleSize), collapse = ", ", sep = "")
  }

  cat("Summary of sample sizes: ", foldSampleSizeText, "\n")

  selectedGenesInfo <- list(selectedGenesNames = NULL)
  if (!is.null(object@finalModel$model$opt.threshold) && object@finalModel$model$opt.threshold != 0){
    selectedGenesInfo$selectedGenesNames <- object@finalModel$model$SelectedGenes[[1]]
    selectedGenesInfo$selectedGenesIdx <- object@finalModel$model$SelectedGenesIndex[[1]]
  }
  featureSelectionText <- if (is.null(selectedGenesInfo$selectedGenesNames)){
    "All features are selected."
  } else {
    paste(length(selectedGenesInfo$selectedGenesIdx), " out of ", p, " features are selected.", sep = "")
  }

  cat("Summary of selected features: ", featureSelectionText, "\n")

  # PART 3. Tuning parameter info.
  # voomDLDA and voomDQDA methods have no tuning parameter.
  tuningResults <- object@tuningResults$results
  if (object@callInfo$method != "voomNSC"){
    cat("\n")
    cat(sprintf("%10s", "Model"), " ", sprintf("%10-s", "Accuracy"), "\n")
    cat(sprintf("%10s", object@callInfo$method), " ", sprintf("%10.7-f", object@finalModel$accuracy), "\n\n")

  } else {
    cat("\n")
    cat(sprintf("%13s", "threshold"), " ", sprintf("%11s", "NonZeroFeat."), " ", sprintf("%10s", "Accuracy"),"\n")
    for (i in 1:nrow(tuningResults)){
      cat(sprintf("%13.5f", tuningResults[i, 1]), " ", sprintf("%12.2f", tuningResults[i, 2]), " ", sprintf("%10.4f", tuningResults[i, 3]), " ", sprintf("%9.4f", tuningResults[i, 4]), "\n")
    }
  }

  # PART 4. Final notes on optimum model parameters.
  if (object@callInfo$method != "voomNSC"){
    cat("There is no tuning parameter for selected method.", "\n")
    cat("Cross-validated model accuracy is given.", "\n\n")

  } else {
    cat("\n")
    cat("The optimum model is obtained when threshold = ",  sprintf("%6.5f", object@tuningResults$opt.threshold), " with an overall accuracy","\n", sep = "")
    cat("of Accuracy = ",  sprintf("%6.4f", object@tuningResults$overallAccuracy), " over folds. On the average ",
        object@tuningResults$results[object@tuningResults$opt.threshold.idx, "nonzero"], " out of ", p, " features was used", "\n", sep = "")
    cat("in the classifier." ,"\n\n", sep = "")

    cat("NOTES: The optimum model is selected using tuning parameter 'threshold' which \n")
    cat(" achieves the highest classification accuracy (Accuracy). The 'NonZeroFeat.'", "\n")
    cat(" is the average number of non-zero features (the selected variables in the", "\n")
    cat(" classification task) over cross-validated folds. As the number of non-zero", "\n")
    cat(" features decreases, the model becomes more sparse.", "\n")
  }

}

#' @rdname show
setMethod(f = "show", signature = signature(object = "voom.train"), definition = show.voom.train)


## PLDA ile egitilen objenin "discrete.train" ekranda goruntulenme yontemi
show.discrete.train <- function(object){
  # x <- object
  n <- ceiling(length(object@inputs$y))
  p <- ceiling(ncol(object@inputs$x))
  nc <- ceiling(length(unique(object@inputs$y)))
  ndigits <- length(unlist(strsplit(as.character(max(n, p, nc)), "")))
  class.names <- object@callInfo$class.names
  cNames <- paste(sapply(class.names, function(x)paste("\'", x, "\'", sep = "")), collapse = ", ")

  # PART 1. Model description and summary of sample/feature sizes. Class labels info.
  cat("\n", modelDescription(object@callInfo$method), " (", object@callInfo$method, ")", sep = "", "\n\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " samples", sep = ""), n), "\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " predictors", sep = ""), p), "\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " classes: ", cNames, "  (Reference category: '", object@callInfo$ref, "')", sep = ""), nc), "\n\n")

  normalizationInfo <- if (object@callInfo$normalize == "none"){
    "Normalization is NOT performed."
  } else if (object@callInfo$normalize == "deseq"){
    "DESeq median ratio."
  } else {
    "Trimmed-mean of M values."
  }

  # PART 2. Normalization and Resampling (Cross validation) info
  cat("Normalization: ", normalizationInfo, "\n", sep = "")
  if (object@callInfo$method != "NBLDA"){
    cat("Power transformation is ", ifelse(object@callInfo$method == "PLDA2", "performed.", " NOT performed."), "\n", sep = "")
  }
  cat("Resampling: Cross-Validated (", object@control$number, " fold, repeated ", object@control$repeats, " times)", sep = "", "\n")

  foldIdx <- if (object@callInfo$method != "NBLDA"){
    object@control$foldIdx
  } else {
    object@control$foldIdx$indexIn
  }

  foldSampleSize <- if (object@callInfo$method != "NBLDA"){
    n - unlist(lapply(unlist(foldIdx, recursive = FALSE), length))
  } else {
    unlist(lapply(foldIdx, length))
  }


  foldSampleSizeText <- if (length(foldSampleSize) > 5){
    paste(c(foldSampleSize[1:5], "..."), collapse = ", ", sep = "")
  } else {
    paste(c(foldSampleSize), collapse = ", ", sep = "")
  }

  cat("Summary of sample sizes: ", foldSampleSizeText, "\n", sep = "")

  if (object@callInfo$method != "NBLDA"){
    selectedGenesInfo <- object@finalModel$SelectedGenes
    featureSelectionText <- if (is.null(selectedGenesInfo$selectedGenesNames)){
      "All features are selected."
    } else {
      paste(length(object@finalModel$SelectedGenes$selectedGenesIdx), " out of ", p, " features are selected.", sep = "")
    }

    cat("Summary of selected features: ", featureSelectionText, "\n")
  }


  # PART 3. Tuning parameter info.
  # PLDA and PLDA2 methods have tuning parameter.
  tuningResults <- object@tuningResults$results
  if (object@callInfo$method != "NBLDA"){
    cat("\n")
    cat(sprintf("%10s", "rho"), " ", sprintf("%8s", "Avg.Error"), " ", sprintf("%14s", "Avg.NonZeroFeat."), " ", sprintf("%8s", "Accuracy"),"\n")
    for (i in 1:nrow(tuningResults)){
      cat(sprintf("%10.5f", tuningResults[i, 1]), " ", sprintf("%9.2f", tuningResults[i, 2]), " ", sprintf("%15.2f", tuningResults[i, 3]), " ", sprintf("%9.4f", tuningResults[i, 4]), "\n")
    }
  } else {
    cat("\n")
    cat(sprintf("%10s", "Model"), " ", sprintf("%10-s", "Accuracy"), "\n", sep = "")
    cat(sprintf("%10s", object@callInfo$method), " ", sprintf("%10.5f", object@finalModel$overallAccuracy), "\n\n", sep = "")

  }

  # PART 4. Final notes on optimum model parameters.
  if (object@callInfo$method != "NBLDA"){
    cat("\n")
    cat("The optimum model is obtained when rho = ",  sprintf("%6.5f", object@tuningResults$bestrho), " with an overall accuracy of","\n", sep = "")
    cat("Accuracy = ",  sprintf("%6.4f", object@tuningResults$bestAccuracy), " over folds. On the average ",
        object@tuningResults$bestNonZeroFeat, " out of ", p, " features was used", "\n", sep = "")
    cat("in the classifier." ,"\n\n", sep = "")

    cat("NOTES: The optimum model is selected using tuning parameter 'rho' which achieves \n")
    cat(" the lowest classification error (Avg.Error) or the highest Accuracy. The ", "\n")
    cat(" classification error is given as the average number of misclassified samples", "\n")
    cat(" over cross-validated folds. Similarly, the 'Avg.NonZeroFeat.' is the average ", "\n")
    cat(" number of non-zero features (the selected variables in the classification task)  ", "\n")
    cat(" over cross-validated folds. As the number of non-zero features decreases, the ", "\n")
    cat(" model becomes more sparse.", "\n")

  } else {
    ### Show method for NBLDA fitted model.
    cat("There is no tuning parameter for selected method.", "\n")
    cat("Cross-validated model accuracy is given.", "\n\n")
  }
}

#' @rdname show
setMethod(f = "show", signature = signature(object = "discrete.train"), definition = show.discrete.train)


#' @title Print method for confusion matrix
#'
#' @description This function prints the confusion matrix of the model.
#'
#' @param x an object of class \code{confMat}
#' @param mode see \code{\link[caret]{print.confusionMatrix}}
#' @param digits see \code{\link[caret]{print.confusionMatrix}}
#' @param \dots further arguments to be passed to \code{print.table}
#'
#' @method print confMat
#' @export
#'
#' @rdname print
print.confMat <- function(x, ..., mode = x$mode, digits = max(3, getOption("digits") - 3)){
  cat("NOTE: The classification table given below is a table of predictions which
      are aggregated over folds by combining fold-out samples in each fold. The cell
      counts are calculated by taking the average of aggregated predictions over
      repeats. We provide 2 tables below: (i) average cell counts, Table 1,
      and (ii) integer valued confusion matrix to be used in hypothesis tests and
      confidence interval estimations. Hence, the performance measures might be
      slightly different as a result of rounding cell counts to nearest integer value.", "\n\n", sep = "")

  # below code chunk is copied from package "caret":  caret:::print.confusionMatrix(x, ...)
  if (is.null(mode)){
    mode <- "sens_spec"
  }
  if (!(mode %in% c("sens_spec", "prec_recall", "everything"))){
    stop("`mode` should be either 'sens_spec', 'prec_recall', or 'everything'")
  }
  cat("Confusion Matrix and Statistics (adapted from package \"caret\")\n\n")
  cat("Table 1: Confusion matrix (Not rounded)\n")
  print(x$tableRounded, ...)
  cat("\n")
  cat("Table 2: Confusion matrix (Rounded for calculations)\n")
  print(x$table, ...)
  cat("\n")

  tmp <- round(x$overall, digits = digits)
  pIndex <- grep("PValue", names(x$overall))
  tmp[pIndex] <- format.pval(x$overall[pIndex], digits = digits)
  overall <- tmp
  accCI <- paste("(", paste(overall[c("AccuracyLower",
                                      "AccuracyUpper")], collapse = ", "), ")", sep = "")
  overallText <- c(paste(overall["Accuracy"]), accCI, paste(overall[c("AccuracyNull",
                                                                      "AccuracyPValue")]), "", paste(overall["Kappa"]),
                   paste(overall["McnemarPValue"]))
  overallNames <- c("Accuracy", "95% CI", "No Information Rate",
                    "P-Value [Acc > NIR]", "", "Kappa", "Mcnemar's Test P-Value")
  if (dim(x$table)[1] > 2) {
    cat("\nOverall Statistics\n")
    overallNames <- ifelse(overallNames == "", "", paste(overallNames,
                                                         ":"))
    out <- cbind(format(overallNames, justify = "right"),
                 overallText)
    colnames(out) <- rep("", ncol(out))
    rownames(out) <- rep("", nrow(out))
    print(out, quote = FALSE)
    cat("\nStatistics by Class:\n\n")
    if (mode == "prec_recall")
      x$byClass <- x$byClass[, !grepl("(Sensitivity)|(Specificity)|(Pos Pred Value)|(Neg Pred Value)",
                                      colnames(x$byClass))]
    if (mode == "sens_spec")
      x$byClass <- x$byClass[, !grepl("(Precision)|(Recall)|(F1)",
                                      colnames(x$byClass))]
    print(t(x$byClass), digits = digits)
  } else {
    if (mode == "prec_recall")
      x$byClass <- x$byClass[!grepl("(Sensitivity)|(Specificity)|(Pos Pred Value)|(Neg Pred Value)",
                                    names(x$byClass))]
    if (mode == "sens_spec")
      x$byClass <- x$byClass[!grepl("(Precision)|(Recall)|(F1)",
                                    names(x$byClass))]
    overallText <- c(overallText, "", format(x$byClass,
                                             digits = digits))
    overallNames <- c(overallNames, "", names(x$byClass))
    overallNames <- ifelse(overallNames == "", "", paste(overallNames,
                                                         ":"))
    overallNames <- c(overallNames, "", "'Positive' Class :")
    overallText <- c(overallText, "", x$positive)
    out <- cbind(format(overallNames, justify = "right"),
                 overallText)
    colnames(out) <- rep("", ncol(out))
    rownames(out) <- rep("", nrow(out))
    out <- rbind(out, rep("", 2))
    print(out, quote = FALSE)
  }

  invisible(x)
}


#' @rdname print
#' @export
setMethod(f = "print", signature = signature(x = "confMat"), definition = print.confMat)


######  3. Method for Predictions: predict(...)   ########

#' @rdname predict
#' @aliases predict,MLSeq-method
#'
#' @export
setMethod(f = "predict", signature = signature(object = "MLSeq"), definition = predict.MLSeq)

######  4. Summary functions  #########
## summary icin bir generic fonksiyon yazilacak.
## MLSeq ve MLSeqModelInfo icin ayri ayri yazilabilir.


######  5. Update method for MLSeq object: update(...)  ######

#' @title Update \code{MLSeq} objects returnd from \code{classify()}
#'
#' @description This function updates the MLSeq object. If one of the options is changed inside MLSeq object, it should be updated
#' to pass its effecs into classification results.
#'
#' @param object a model of \code{MLSeq} class returned by \code{\link{classify}}
#' @param env an environment. Define the environment where the trained model is stored.
#' @param ... optional arguements passed to \code{\link{classify}} function.
#'
#' @return same object as an MLSeq object returned from \code{classify}.
#'
#' @note When an \code{MLSeq} object is updated, new results are updated on the given object. The results before update process are
#' lost when update is done. To keep the results before update, one should copy the MLSeq object to a new object in global environment.
#'
#' @seealso \code{\link{classify}}
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
#' # test set
#' data.test <- data[ ,ind]
#' data.test <- as.matrix(data.test + 1)
#' classts <- data.frame(condition=class[ind, ])
#'
# test set in S4
#' data.testS4 <- DESeqDataSetFromMatrix(countData = data.test,
#'                                       colData = classts, formula(~ 1))
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "rpart",
#'           ref = "T", preProcessing = "deseq-vst",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#' cart
#'
#' # Change classification model into "Random Forests" (rf)
#' method(cart) <- "rf"
#' rf <- update(cart)
#'
#' rf
#'}
#'
#' @name update
#' @rdname update
#'
#' @method update MLSeq
update.MLSeq <- function(object, ..., env = .GlobalEnv){

  ### Use this object name to update model without using left assignment.
  object.name <- deparse(substitute(object))

  if (!isUpdated(object) & isModified(object)){
    # MLSeq object is modified but not updated.
    # Update MLSeq object using parameters stored in metaData

    modifiedElements <- metaData(object)@modified.elements$modifiedElements

    # Update classify(...) object
    ## If classification method or preProcessing is changed.
    if (any(c("method", "preProcessing", "control", "normalization") %in% modifiedElements)){
      inputData <- metaData(object)@rawData.DESeqDataSet

      #### CONTROL objesi kullanilan sınıflamaya göre farklı olacak. Burada trained(model) sınıfına göre
      control.object <- control(object)
      if (!("controlClass" %in% names(control.object))){
        ## trControl ile gelen obje içinde class yer almadığı için bunu kontrol etmek gerekiyor.
        control.object$controlClass <- "train.control"
      }

      if (method(object) %in% c("voomDLDA", "voomDQDA", "voomNSC") & control.object$controlClass != "voom.control"){
        stop(warning("Incorrect elements in 'control' argument. It should be defined using 'voomControl(...)' function."))
      } else if (method(object) %in% c("PLDA", "PLDA2", "NBLDA") & control.object$controlClass != "discrete.control"){
        stop(warning("Incorrect elements in 'control' argument. It should be defined using 'discreteControl(...)' function."))
      } else if (!(method(object) %in% c("PLDA", "PLDA2", "NBLDA", "voomDLDA", "voomDQDA", "voomNSC")) & control.object$controlClass != "train.control"){
        stop(warning("Incorrect elements in 'control' argument. It should be defined using 'trainControl(...)' function."))
      }

      #### if bloğu gerekiyor.
      if (!(method(object) %in% c("PLDA", "PLDA2", "NBLDA", "voomDLDA", "voomDQDA", "voomNSC"))){   ### train
        if ("control" %in% modifiedElements){
          control.object <- control(object)
          folds <- foldIndex(n = dim(inputData)[2], nFolds = control.object$number, repeats = control.object$repeats)
          control.object$index <- folds$indexIn
          control.object$indexOut <- folds$indexOut

          ## Remove seeds to prevent seed error.
          control.object$seeds <- NULL

        } else {
          control.object <- control(object)

          ## Remove seeds to prevent seed error.
          control.object$seeds <- NULL
        }

        ### Get call from trained model.
        call.args <- object@modelInfo@call$args
        call.args$data <- inputData
        call.args$method <- method(object)
        call.args$ref <- ref(object)
        call.args$preProcessing <- preProcessing(object)
        call.args$control <- control.object

      } else if (method(object) %in% c("voomDLDA", "voomDQDA", "voomNSC")){ # voom.train
        control.object <- control(object)

        ### Get call from trained model. (voom.train)
        call.args <- object@modelInfo@call$args
        call.args$data <- inputData
        call.args$normalize <- normalization(object)
        call.args$method <- method(object)
        call.args$ref <- ref(object)
        call.args$control <- control.object

      } else if (method(object) %in% c("PLDA", "PLDA2", "NBLDA")){ # discrete.train
        if ("control" %in% modifiedElements){
          if (method(object) != "NBLDA"){
            control.object <- control(object)

            y <- colData(object@metaData@rawData.DESeqDataSet)[ ,object@modelInfo@call$args$class.labels]

            ###### FOLD KODLARI STANDART BIR FORMDA OLACAK. CARET ve PLDA'da kullanilan fold fonksiyonlari farkli.
            foldIdx <- lapply(1:control.object$repeats, function(u){
              tmp.fold <- balanced.folds(y = y, nfolds = control.object$number)
              names(tmp.fold) <- paste("Fold.", 1:length(tmp.fold), sep = "")
              tmp.fold
            })

            names(foldIdx) <- paste("Repeat.", 1:control.object$repeats, sep = "")
            control.object$foldIdx <- foldIdx

          } else {
            control.object <- control(object)

            foldIdx <- foldIndex(n = dim(inputData)[2], nFolds = control.object$number, repeats = control.object$repeats)
            control.object$foldIdx <- foldIdx
          }

        } else {
          control.object <- control(object)
        }

        ### Get call from trained model.
        call.args <- object@modelInfo@call$args
        call.args$data <- inputData
        call.args$normalize <- normalization(object)
        call.args$method <- method(object)
        call.args$ref <- ref(object)
        call.args$control <- control.object

      }

      call.args <- c(call.args, ...)
      new.object <- do.call(what = "classify", args = call.args)

      isModified(new.object) <- TRUE
      isUpdated(new.object) <- TRUE
      new.object@metaData@modified.elements$modifiedElements <- NULL

    } else if ("ref" %in% modifiedElements & length(modifiedElements) == 1){
      new.object <- object
      confM <- confusionMat(object)
      tbl <- confM$table

      confM.updated <- confusionMatrix(tbl, positive = ref(object))
      confM.updated$tableRounded <- confM$tableRounded
      attr(confM.updated, "class") <- "confMat"
      new.object@modelInfo@confusionMat <- confM.updated

      isModified(new.object) <- TRUE
      isUpdated(new.object) <- TRUE

      new.object@metaData@modified.elements$modifiedElements <- NULL
    }

    cat("\n\n Update is successfull... \n\n")
  } else {
    warning("Either 'MLSeq' object is not modified or already updated. No update is performed on given object.")
  }

  # assign(object.name, new.object, envir = env)   ## Bu kısım iptal edilecek.
  return(new.object)
}


#' @rdname update
#' @export
setMethod(f = "update", signature = signature(object = "MLSeq"), definition = update.MLSeq)



######### 6. Plots ########
# PLOT icin graphics altındaki plot(...) generic fonksiyonunu kullanalim. caret'in plot fonksiyonunu
# maskeliyoruz kendi generic fonksiyonumuzu yazdigimiz durumda.

#' @title Plot accuracy results from 'MLSeq' object
#'
#' @description This generic function is used to plot accuracy results from 'MLSeq' object returned by
#' \code{classify} function.
#'
#' @docType methods
#' @name plot
#' @rdname plot
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Bernd Klaus and Ahmet Ozturk
#'
#' @param x an \code{MLSeq} object returned from \code{classify} function.
#' @param y this parameter is not used. Deprecated.
#' @param ... further arguements. Deprecated.
#'
#' @method plot MLSeq
#' @export
#'
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom caret 'plot.train'
#' @importFrom grDevices rgb
plot.MLSeq <- function(x, y, ...){
  object <- x
  object.trained <- trained(object)
  if (method(object) %in% c("PLDA", "PLDA2")){
    bestRho <- object@modelInfo@trainedModel@tuningResults$bestrho
    bestAccuracy <- object@modelInfo@trainedModel@tuningResults$bestAccuracy
    plotData <- object.trained@tuningResults$results
    # plot(plotData$rho, plotData$Accuracy, xlab = "rho", ylab = "Accuracy", type = "l", pch = 16)
    ggplot(data = plotData, aes_string(x = "rho", y = "Accuracy")) +
      geom_vline(xintercept = bestRho, lty = 2) +
      geom_hline(yintercept = bestAccuracy, lty = 2) +
      geom_point(colour = rgb(0.1, 0.3, 0.8, alpha = 0.8)) +
      geom_line(colour = rgb(0.1, 0.3, 0.8, alpha = 0.5)) +
      theme_bw(base_size = 12) +
      xlab("Tuning parameter (rho)") +
      theme(axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
            axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
            axis.text.x = element_text(margin = margin(5, 0, 0, 0)),
            axis.text.y = element_text(margin = margin(0, 5, 0, 0)))
  } else if (method(object) %in% c("voomDLDA", "voomDQDA", "voomNSC")){
    bestThreshold <- object@modelInfo@trainedModel@tuningResults$opt.threshold
    bestAccuracy <- object@modelInfo@trainedModel@tuningResults$overallAccuracy
    plotData <- object@modelInfo@trainedModel@tuningResults$results
    # plot(plotData$rho, plotData$Accuracy, xlab = "rho", ylab = "Accuracy", type = "l", pch = 16)
    ggplot(data = plotData, aes_string(x = "threshold", y = "accuracy")) +
      geom_vline(xintercept = bestThreshold, lty = 2) +
      geom_hline(yintercept = bestAccuracy, lty = 2) +
      geom_point(colour = rgb(0.1, 0.3, 0.8, alpha = 0.8)) +
      geom_line(colour = rgb(0.1, 0.3, 0.8, alpha = 0.5)) +
      theme_bw(base_size = 12) +
      xlab("Tuning parameter (threshold)") +
      ylab("Accuracy") +
      theme(axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
            axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
            axis.text.x = element_text(margin = margin(5, 0, 0, 0)),
            axis.text.y = element_text(margin = margin(0, 5, 0, 0)))
  } else if (class(trained(object)) == "train"){
    plot.train(trained(object))
  }
}

#' @rdname plot
#' @export
setMethod(f = "plot", signature = signature(x = "MLSeq"), definition = plot.MLSeq)

