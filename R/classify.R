globalVariables("ii")  ## to remove warnings from 'foreach' package.

#' @include all_generics.R
#' @include voomFunctions.R

#' @title Fitting classification models to sequencing data
#'
#' @description This function fits classification algorithms to sequencing data and measures model performances using various statistics.
#'
#' @details MLSeq consists both microarray-based and discrete-based classifiers along with the preprocessing approaches. These approaches
#' include both normalization techniques, i.e. deseq median ratio (Anders et al., 2010) and trimmed mean of M values (Robinson et al., 2010)
#' normalization methods, and the transformation techniques, i.e. variance- stabilizing transformation (vst)(Anders and Huber, 2010),
#' regularized logarithmic transformation (rlog)(Love et al., 2014), logarithm of counts per million reads (log-cpm)(Robinson et al., 2010)
#' and variance modeling at observational level (voom)(Law et al., 2014). Users can directly upload their raw RNA-Seq count data, preprocess
#' their data, build one of the numerous classification models, optimize the model parameters and evaluate the model performances.
#'
#' MLSeq package consists of a variety of classification algorithms for the classification of RNA-Seq data. These classifiers are
#' categorized into two class: i) microarray-based classifiers after proper transformation, ii) discrete-based classifiers. First option
#' is to transform the RNA-Seq data to bring it hierarchically closer to microarrays and apply microarray-based algorithms. These methods
#' are implemented from the caret package. Run availableMethods() for a list of available methods. Note that voom transformation both
#' exports transformed gene-expression matrix as well as the precision weight matrices in same dimension. Hence, the classifier should
#' consider these two matrices. Zararsiz (2015) presented voom-based diagonal discriminant classifiers and the sparse voom-based nearest
#' shrunken centroids classifier. Second option is to build new discrete-based classifiers to classify RNA-Seq data. Two methods are
#' currently available in the literature. Witten (2011) considered modeling these counts with Poisson distribution and proposed sparse
#' Poisson linear discriminant analysis (PLDA) classifier. The authors suggested a power transformation to deal with the overdispersion
#' problem. Dong et al. (2016) extended this approach into a negative binomial linear discriminant analysis (NBLDA) classifier.
#' More detailed information can be found in referenced papers.
#'
#' @param data a \code{DESeqDataSet} object, see the constructor functions \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}}, \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSetFromMatrix}},
#' \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSetFromHTSeqCount}} in DESeq2 package.
#'
#' @param method a character string indicating the name of classification method. Methods are implemented from the \code{caret} package.
#' Run \code{availableMethods()} for a list of available methods.
#'
#' @param normalize a character string indicating the type of normalization. Should be one of 'deseq', 'tmm' and 'none'. Default is 'deseq'. This option
#' should be used with discrete and voom-based classifiers since no transformation is applied on raw counts. For caret-based classifiers,
#' the argument 'preProcessing' should be used.
#'
#' @param preProcessing a character string indicating the name of the preprocessing method. This option consists both the normalization and
#' transformation of the raw sequencing data. Available options are:
#' \itemize{
#' \item \code{deseq-vst}: Normalization is applied with deseq median ratio method. Variance stabiling transformation is applied to the normalized data.
#' \item \code{deseq-rlog}: Normalization is applied with deseq median ratio method. Regularized logarithmic transformation is applied to the normalized data.
#' \item \code{deseq-logcpm}: Normalization is applied with deseq median ratio method. Log of counts-per-million transformation is applied to the normalized data.
#' \item \code{tmm-logcpm}: Normalization is applied with trimmed mean of M values (TMM) method. Log of counts-per-million transformation is applied to the
#' normalized data.
#' \item \code{logcpm}: Normalization is not applied. Log of counts-per-million transformation is used for the raw counts.
#' }
#' \bold{IMPORTANT}: See Details for further information.
#'
#' @param class.labels a character string indicating the column name of colData(...). Should be given as "character". The column from colData()
#' which matches with given column name is used as class labels of samples. If NULL, first column is used as class labels. Default is NULL.
#'
#' @param control a list including all the control parameters passed to model training process. This arguement should be defined using wrapper functions
#' \code{\link[caret]{trainControl}} for caret-based classifiers, \code{\link{discreteControl}} for discrete classifiers (PLDA, PLDA2 and NBLDA) and
#' \code{\link{voomControl}} for voom-based classifiers (voomDLDA, voomDQDA and voomNSC). See related functions for further details.
#'
#' @param B an integer. It is the number of bootstrap samples for bagging classifiers, for example "bagFDA" and "treebag". Default is 25.
#'
#' @param ref a character string indicating the user defined reference class. Default is \code{NULL}. If NULL is selected,
#' first category of class labels is used as reference.
#'
#' @param \dots optional arguments passed to selected classifiers.
#'
#' @return an \code{MLSeq} object for trained model.
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Bernd Klaus, Ahmet Ozturk and Ahmet Ergun Karaagaoglu
#'
#' @references
#'
#' Kuhn M. (2008). Building predictive models in R using the caret package. Journal of Statistical Software, (http://www.jstatsoft.org/v28/i05/)
#'
#' Anders S. Huber W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11:R106
#'
#' Witten DM. (2011). Classification and clustering of sequencing data using a poisson model. The Annals of Applied Statistics, 5(4), 2493:2518
#'
#' Law et al. (2014) Voom: precision weights unlock linear model analysis tools for RNA-Seq read counts, Genome Biology, 15:R29, doi:10.1186/gb-2014-15-2-r29
#'
#' Witten D. et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical biomarker analysis in a collection of
#' cervical tumours and matched controls. BMC Biology, 8:58
#'
#' Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-Seq data. Genome Biology,
#' 11:R25, doi:10.1186/gb-2010-11-3-r25
#'
#' M. I. Love, W. Huber, and S. Anders (2014). Moderated estimation of fold change and dispersion for rna-seq data with deseq2. Genome Biol,
#' 15(12):550,. doi: 10.1186/s13059-014-0550-8.
#'
#' Dong et al. (2016). NBLDA: negative binomial linear discriminant analysis for rna-seq data. BMC Bioinformatics, 17(1):369, Sep 2016.
#' doi: 10.1186/s12859-016-1208-1.
#'
#' Zararsiz G (2015). Development and Application of Novel Machine Learning Approaches for RNA-Seq Data Classification. PhD thesis,
#' Hacettepe University, Institute of Health Sciences, June 2015.
#'
#' @keywords RNA-seq classification
#'
#' @seealso \code{\link{predictClassify}}, \code{\link[caret]{train}}, \code{\link[caret]{trainControl}},
#' \code{\link{voomControl}}, \code{\link{discreteControl}}
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
#' ## Number of repeats (repeats) might change model accuracies
#' ## 1. caret-based classifiers:
#' # Random Forest (RF) Classification
#'  rf <- classify(data = data.trainS4, method = "rf",
#'          preProcessing = "deseq-vst", ref = "T",
#'          control = trainControl(method = "repeatedcv", number = 5,
#'                                 repeats = 2, classProbs = TRUE))
#' rf
#'
#' # 2. Discrete classifiers:
#' # Poisson Linear Discriminant Analysis
#' pmodel <- classify(data = data.trainS4, method = "PLDA", ref = "T",
#'                    class.labels = "condition",normalize = "deseq",
#'                    control = discreteControl(number = 5, repeats = 2,
#'                                              tuneLength = 10, parallel = TRUE))
#' pmodel
#'
#' # 3. voom-based classifiers:
#' # voom-based Nearest Shrunken Centroids
#' vmodel <- classify(data = data.trainS4, normalize = "deseq", method = "voomNSC",
#'                    class.labels = "condition", ref = "T",
#'                    control = voomControl(number = 5, repeats = 2, tuneLength = 10))
#' vmodel
#' }
#'
#' @name classify
#' @rdname classify
#'
#' @import methods DESeq2 Biobase
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom stats model.matrix predict relevel xtabs quantile as.dist rexp rnbinom rnorm runif sd
#' @importFrom Biobase ExpressionSet exprs
#' @importFrom limma voom
#' @importFrom caret bagControl confusionMatrix svmBag trainControl bag train
#' @importFrom SummarizedExperiment colData 'colData<-'
#' @importFrom graphics segments
#'
#' @export
classify <- function(data, method = "rpart", B = 25, ref = NULL, class.labels = NULL,
                     preProcessing = c("deseq-vst", "deseq-rlog", "deseq-logcpm", "tmm-logcpm", "logcpm"),
                     normalize = c("deseq", "TMM", "none"), control = NULL, ...){

  # Define the name of class.labels from DESeqDataSet
  if (is.null(class.labels)){
    class.labels <- colnames(colData(data))[1]
  }

  if (!is.null(ref)){
    if (!is.character(ref)){
      stop("Reference class should be \"character\"")
    }
  }

  # If reference category is not specified, first category is selected as reference category
  if (is.null(ref)){
    ref = levels(colData(data)[ ,class.labels])[1]
  }

  # Data class olarak DESeqDataSet'e bagli kalmayalim.
  if (class(data)[1] != "DESeqDataSet") {
    stop("Data should be a \"DESeqDataSet Object\" of S4 class.")
  }

  if (is.null(method)){
    stop("Classification method is not specified.")
  }

  if (!(method %in% availableMethods())){
    stop("Requested method is not available in \"MLSeq\".")
  }

  if (method %in% c("voomDLDA", "voomDQDA", "voomNSC")){
    ## voom-based classifications:
    result <- classify.voom(data = data, normalize = normalize, method = method,
                            class.labels = class.labels, ref = ref, control = control)
  } else if (method %in% c("PLDA", "PLDA2", "NBLDA")){
    result <- classify.discrete(data = data, method = method, normalize = normalize, ref = ref,
                                class.labels = class.labels, control = control)
  } else {
    result <- classify.continous(data, method = method, B = B, ref = ref, class.labels = class.labels,
                                 preProcessing = preProcessing, control = control, ...)
  }

  return(result)
}

# data = data.trainS4
# method = "svmRadial"
# class.labels = "condition"
# ref = "T"
# preProcessing = "deseq-vst"
# control = trainControl(method = "repeatedcv", number = 5,
#                        repeats = 2, classProbs = TRUE)
# tuneLength = 15


#' @importFrom caret 'confusionMatrix.train'
classify.continous <- function(data, method = "rpart", B = 25, ref = NULL, class.labels = NULL,
                               preProcessing = c("deseq-vst", "deseq-rlog", "deseq-logcpm", "tmm-logcpm", "logcpm"),
                               control = trainControl(method = "repeatedcv", number = 5, repeats = 10), ...){
  ## Classification function for microarray type classifiers (caret's built-in library is imported here.)

  control$controlClass <- "train.control"
  call <- match.call()
  arg.list <- list(data = data, method = method, B = B, ref = ref, class.labels = class.labels,
                   preProcessing = preProcessing,
                   control = control, ...)

  rawData <- data
  transformedData <- NULL

  preProcessing <- match.arg(preProcessing)

  if (preProcessing == "deseq-vst"){
    normalization <- "deseq"
    transformation <- "vst"

    dataSF <- estimateSizeFactors(data)    # Estimate Size factors:
    dataDisp <- estimateDispersions(dataSF, fitType = "local")  # Estimate dispersions:
    transformedData <- varianceStabilizingTransformation(dataDisp, fitType = "local")

    dataVST <- t(as.matrix(assay(transformedData)))
    input <- dataVST   ## Input data from transformed expression data.
    output <- colData(data)[ ,class.labels]   ## Output: class labels of samples.

    trainParameters <- list(disperFunc = dispersionFunction(dataDisp),
                            sizeFactors = sizeFactors(dataDisp))

  } else if (preProcessing == "deseq-rlog"){
    normalization <- "deseq"
    transformation <- "rlog"

    dataSF <- estimateSizeFactors(data)    # Estimate Size factors:
    dataDisp <- estimateDispersions(dataSF, fitType = "local")  # Estimate dispersions:

    transformedData = rlog(dataDisp, blind = FALSE)  ### RLOG donusumu uygulanmis olan bu nesnenin de MLSeq objesi icerisinde predictClassify'a gonderilmesi gerekiyor.
    dataRLOG = t(assay(transformedData))
    input <- dataRLOG   ## Input data from transformed expression data.
    output <- colData(data)[ ,class.labels]   ## Output: class labels of samples.

    # dataexpTrainRLOG = t(assay(trS4expRLOG))
    # betaPriorVar = attr(trS4expRLOG,"betaPriorVar")
    # intercept = mcols(trS4expRLOG)$rlogIntercept

    trainParameters <- list(betaPriorVar = attr(dataRLOG,"betaPriorVar"),
                            intercept = attr(dataRLOG, "intercept"),
                            dispFit = mcols(transformedData)$dispFit)

  } else if (preProcessing == "deseq-logcpm"){
    normalization <- "deseq"
    transformation <- "logcpm"

    rawCounts = counts(data, normalized = FALSE)
    countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
    countsDGE.normalized <- calcNormFactors(countsDGE, method = "RLE")   ## RLE: DESeq mantigi ile normalize ediyor.
    countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek.

    rawData <- countsDGE
    classes <- colData(data)[ ,class.labels]
    rawData[[class.labels]] <- classes

    transformedData <- list(normalizedData = countsDGE.normalized, transformedData = countsDGE.transformed)
    transformedData[[class.labels]] <- classes

    input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
    output <- rawData[[class.labels]]   ## Output: class labels of samples.

    trainParameters <- list(NULL)

  } else if (preProcessing == "tmm-logcpm"){
    normalization <- "tmm"
    transformation <- "logcpm"

    rawCounts = counts(data, normalized = FALSE)
    countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
    countsDGE.normalized <- calcNormFactors(countsDGE, method = "TMM")   ## RLE: DESeq mantigi ile normalize ediyor.
    countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek.

    rawData <- countsDGE
    classes <- colData(data)[ ,class.labels]
    rawData[[class.labels]] <- classes

    transformedData <- list(normalizedData = countsDGE.normalized, transformedData = countsDGE.transformed)
    transformedData[[class.labels]] <- classes

    input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
    output <- rawData[[class.labels]]   ## Output: class labels of samples.

    # This function is used to find reference sample.
    # Codes are copied from edgeR and modified here.
    findRefSample <- function (rawCounts, lib.size, p = 0.75){
      y <- t(t(rawCounts) / lib.size)
      f75 <- apply(y, 2, function(x){
        quantile(x, p = p)
      })
      refSample <- which.min(abs(f75-mean(f75)))
      return(refSample)
    }

    refSampleID <- findRefSample(rawCounts, lib.size = colSums(rawCounts), p = 0.75)
    trainParameters <- list(refSample = rawCounts[ ,refSampleID])

  } else if (preProcessing == "logcpm"){
    normalization <- "none"
    transformation <- "logcpm"

    rawCounts = counts(data, normalized = FALSE)
    countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
    countsDGE.normalized <- calcNormFactors(countsDGE, method = "none")   ## RLE: DESeq mantigi ile normalize ediyor.
    countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek.

    rawData <- countsDGE
    classes <- colData(data)[ ,class.labels]
    rawData[[class.labels]] <- classes

    transformedData <- list(normalizedData = countsDGE.normalized, transformedData = countsDGE.transformed)
    transformedData[[class.labels]] <- classes

    input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
    output <- rawData[[class.labels]]  ## Output: class labels of samples.

    trainParameters <- list(NULL)
  }

  trainedModel <- train(input, output, method = method, trControl = control, ...)

  tmp <- confusionMatrix.train(trainedModel, norm = "none")
  tbl <- tmp$table / control$repeats

  confM = confusionMatrix(round(tbl, 0), positive = ref)
  confM$tableRounded <- round(tbl, 2)
  confM$table <- round(tbl, 0)

  attr(confM, "class") <- "confMat"

  ## All the information about classification model.
  modelInfo <- new("MLSeqModelInfo",
                   method = method,
                   transformation = transformation,
                   normalization = normalization,
                   preProcessing = preProcessing,
                   ref = ref,
                   control = trainedModel$control,
                   confusionMat = confM,
                   trainedModel = trainedModel,
                   trainParameters = trainParameters,
                   call = list(call = call, args = arg.list))

  result = new("MLSeq",
               modelInfo = modelInfo,
               inputObject = list(rawData = rawData, transformedData = transformedData))

  ## Metadata for MLSeq object.
  meta <- new("MLSeqMetaData", classLabel = class.labels, rawData.DESeqDataSet = data)
  result@metaData <- meta

  result
}

### voom classifiers:
### TODO's
# 1. method argumaninda yer alan runAvailableMethods() fonksiyonunda bir "type" argumani olacak. Bu argumana kisi "voom", "microarray" ...
# girisler yapinca uc farkli classify fonksiyonuna ait methodlar listelenebilecek. setter ve getterlar'da buna gore belirlenecek.

### Classification for voom-based methods.

# data <- data.trainS4
# normalize = "deseq"   # c("deseq", "tmm", "none")
# method = "voomNSC"   # c("voomNSC", "voomDLDA", "voomDQDA")
# class.labels = "condition"
# ref = "T"
# control = voomControl(number = 3, repeats = 2)

#' @importFrom plyr ldply
classify.voom <- function(data, normalize = c("deseq", "TMM", "none"), method = c("voomNSC", "voomDLDA", "voomDQDA"),
                         class.labels = NULL, ref = NULL, control = voomControl(), ...){

  normalize <- match.arg(normalize)
  method <- match.arg(method)
  call <- match.call()
  arg.list <- list(data = data, normalize = normalize, method = method,
                   class.labels = class.labels, ref = ref, control = control, ...)

  ### Calculate train set performances:
  folds <- foldIndex(n = nrow(colData(data)), nFolds = control$number, repeats = control$repeats)

  ## k-fold r-repeat:
  # train ve predict işlemleri tüm foldlar ve repeatler için yapılacak.
  if (method %in% c("voomDLDA", "voomDQDA")){
    # Samples for train in this fold.
    counts <- counts(data)
    output <- colData(data)[ ,class.labels]

    trainResults.AllFolds <- lapply(folds$indexIn, function(idx){
      input.foldIn <- counts[ ,idx]
      output.foldIn <- output[idx]

      # Samples for validation in this fold.
      input.foldOut <- counts[ ,-idx]
      output.foldOut <- output[-idx]

      trainResults.fold <- voomDDA.train(input = input.foldIn, output = output.foldIn,
                                         pooled.var = ifelse(method == "voomDLDA", TRUE, FALSE),
                                         normalize = normalize, class.labels = class.labels)

      ### Predictions for train set:
      pred.fold <- factor(predict.voomDDA(object = trainResults.fold, newdata = input.foldOut),
                          levels = levels(output))
      tbl <- table(Prediction = pred.fold, Reference = output.foldOut)
      acc.fold <- sum(diag(tbl))/sum(tbl)

      return(list(accuracies = acc.fold, predictions = pred.fold, confTables = tbl))
    })

    acc <- ldply(lapply(trainResults.AllFolds, function(x)x$accuracies), rbind, .id = "Folds")
    colnames(acc)[2] <- "Accuracy"

    overallAccuracy <- mean(acc[ ,"Accuracy"])

    ## Final Model
    trainedOverall <- voomDDA.train(input = counts(data), output = colData(data)[ ,class.labels],
                                    pooled.var = ifelse(method == "voomDLDA", TRUE, FALSE), normalize = normalize,
                                    class.labels = class.labels)

    # predOverall <- factor(predict.voomDDA(object = trainedOverall, newdata = counts(data)))
    # tblOverall <- table(Predicted = predOverall, Actual = colData(data)[ ,class.labels])

    # Confusion matrix over folds and repeats
    aggregatedTable <- lapply(trainResults.AllFolds, function(x){
      x$confTables
    })

    tmp <- as.matrix(aggregatedTable[[1]])
    for (i in 2:length(aggregatedTable)){
      tmp <- tmp + as.matrix(aggregatedTable[[i]])
    }

    tmp <- tmp / control$repeats

    confM <- confusionMatrix(round(tmp, 0), positive = ref)
    confM$tableRounded <- round(tmp, 2)
    confM$table <- round(tmp, 0)
    attr(confM, "class") <- "confMat"

    results.tune <- list()  ## No tuning results for "voomDLDA" and "voomDQDA"
  } else {  ### voom NSC calculations:
    fit.initial <- voomNSC.train(counts = counts(data), conditions = colData(data)[ ,class.labels],
                                 normalization = normalize, n.threshold = control$tuneLength, getInitialThresholds = FALSE)

    initial.thresholds <- fit.initial$threshold
    counts <- counts(data)
    output <- colData(data)[ ,class.labels]

    trainResults.AllFolds <- lapply(folds$indexIn, function(idx){
      # # Samples for train in this fold.
      testdata <- counts[ ,-idx]
      test.conditions <- output[-idx]

      trainResults.fold <- voomNSC.train(counts = counts, conditions = output,
                                         normalization = normalize, thresholds = initial.thresholds, idxIn = idx)


      # ### Predictions for train set:
      # pred.fold <- factor(predict.voomNSC(fit = trainResults.fold, newdata = testdata))
      # tbl <- table(test.conditions, pred.fold)
      # acc.fold <- sum(diag(tbl))/sum(tbl)
      #
      # return(list(accuracies = acc.fold, predictions = pred.fold))

      trainResults.fold[["test.conditions"]] <- test.conditions
      return(trainResults.fold)
    })

    ### Find optimum threshold value:
    for (i in 1:length(trainResults.AllFolds)){
      if (i == 1){
        modelResAll <- trainResults.AllFolds[[1]]$modelRes
      } else {
        modelResAll <- modelResAll + trainResults.AllFolds[[i]]$modelRes
      }
    }

    modelResAll <- modelResAll / length(folds$indexIn)
    # opt.threshold.idx <- max(which(modelResAll$avg.errors == min(modelResAll$avg.errors)))
    opt.threshold.idx <- max(which(modelResAll$accuracy == max(modelResAll$accuracy)))
    opt.threshold <- max(modelResAll$threshold[opt.threshold.idx])

    # Return the fold accuracies using optimum thresholds:
    acc.Folds <- lapply(trainResults.AllFolds, function(x){
      Actual <- x$test.conditions
      Preds <- x$predictions[ ,opt.threshold.idx]
      tbl <- table(Prediction = Preds, Reference = Actual)
      list(errors = x$errors[opt.threshold.idx], predictions = x$predictions[ ,opt.threshold.idx],
           accuracy = x$accuracy[opt.threshold.idx], confTables = tbl)
    })


    acc <- ldply(lapply(acc.Folds, function(x)x$accuracy), rbind, .id = "Folds")
    colnames(acc)[2] <- "Accuracy"

    overallAccuracy <- mean(acc[ ,"Accuracy"])


    ### BU KISIM PREDICT KISMINA EKLENECEK.
    ## Confusion Matrix mantığı caret ile benzer şekilde düzenlenecek.
    trainedOverall <- voomNSC.train(counts = counts(data), conditions = colData(data)[ ,class.labels],
                                    normalization = normalize, idxIn = NULL, thresholds = opt.threshold)

    # trainedOverall <- fit.initial

    # Confusion matrix over folds and repeats
    aggregatedTable <- lapply(acc.Folds, function(x){
      x$confTables
    })

    tmp <- as.matrix(aggregatedTable[[1]])
    for (i in 2:length(aggregatedTable)){
      tmp <- tmp + as.matrix(aggregatedTable[[i]])
    }

    tmp <- tmp / control$repeats

    confM <- confusionMatrix(round(tmp, 0), positive = ref)
    confM$tableRounded <- round(tmp, 2)
    confM$table <- round(tmp, 0)
    attr(confM, "class") <- "confMat"

    results.tune <- list(results = modelResAll, opt.threshold = opt.threshold,
                         opt.threshold.idx = opt.threshold.idx, overallAccuracy = overallAccuracy)
  }

  # ## voomTrained adı ile bir S4 class yazılacak.
  trainedModel <- new("voom.train",
                      weigtedStats = trainedOverall$weightedStats,
                      foldInfo = list(foldAccuracies = acc,
                                      foldPredictions = lapply(trainResults.AllFolds, function(x)x$predictions),
                                      foldIndex = folds),
                      control = control,
                      tuningResults = results.tune,
                      finalModel = list(model = trainedOverall, accuracy = overallAccuracy, finalPredictions = list()),
                      callInfo = list(method = method, normalize = normalize,
                                      class.labels = class.labels, class.names = as.character(unique(colData(data)[ ,class.labels])),
                                      ref = ref))


  ## All the information about classification model.
  modelInfo <- new("MLSeqModelInfo",
                   method = method,
                   transformation = "voom",
                   normalization = normalize,
                   preProcessing = NA_character_,
                   ref = ref,
                   control = trainedModel@control,
                   confusionMat = confM,
                   trainedModel = trainedModel,
                   call = list(call = call, args = arg.list))

  result = new("MLSeq",
               modelInfo = modelInfo,
               inputObject = list(rawData = trainedOverall$rawData, transformedData = trainedOverall$transformedData))

  ## Metadata for MLSeq object.
  meta <- new("MLSeqMetaData", classLabel = class.labels, rawData.DESeqDataSet = data)
  result@metaData <- meta

  return(result)
}



####### PLDA, NBLDA ##########

#' Define controlling parameters for discrete classifiers (NBLDA and PLDA)
#'
#' This function sets the control parameters for discrete classifiers (PLDA and NBLDA) while training the model.
#'
#' @param method validation method. Support repeated cross validation only ("repeatedcv").
#' @param number a positive integer. Number of folds.
#' @param repeats a positive integer. Number of repeats.
#' @param tuneLength a positive integer. If there is a tuning parameter in the classifier, this value
#' is used to define total number of tuning parameter to be searched.
#' @param rho a single numeric value. This parameter is used as tuning parameter in PLDA classifier.
#' It does not effect NBLDA classifier.
#' @param rhos a numeric vector. If optimum parameter is searched among given values, this option shpould be used.
#' @param beta parameter of Gamma distribution. See PLDA for details.
#' @param prior prior probabilities of each class
#' @param alpha a numeric value in the interval 0 and 1. It is used to apply power transformation through PLDA method.
#' @param truephi a numeric value. If true value of genewise dispersion is known and constant for all genes, this
#' parameter should be used.
#' @param foldIdx a list including the fold indexes.
#' @param parallel if TRUE, parallel computing is performed.
#' @param ... further arguments. Deprecated.
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Bernd Klaus, Ahmet Ozturk and Ahmet Ergun Karaagaoglu
#'
#' @keywords RNA-seq classification
#'
#' @seealso \code{\link{classify}}, \code{\link[caret]{trainControl}}, \code{\link{discreteControl}}
#'
#' @examples
#' 1L
#'
#' @rdname discreteControl
#'
#' @export
discreteControl <- function(method = "repeatedcv", number = 5, repeats = 10, rho = NULL, rhos = NULL, beta = 1,
                            prior = NULL, alpha = NULL, truephi = NULL, foldIdx = NULL, tuneLength = 30,
                            parallel = FALSE, ...){
  list(method = method, number = number, repeats = repeats, rho = rho, rhos = rhos, beta = beta,
       prior = prior, alpha = alpha, truephi = truephi, foldIdx = foldIdx, tuneLength = tuneLength,
       parallel = parallel, controlClass = "discrete.control", ...)
}


#' @importFrom foreach foreach '%do%' '%dopar%'
#' @importFrom utils globalVariables
trainPLDA <- function(x, y, rhos = NULL, beta = 1, type = c("mle", "deseq", "quantile", "none", "TMM"),
                      folds = NULL, transform = TRUE, alpha = NULL, prior = NULL, nfolds = 5, repeats = 10,
                      tuneLength = 10, parallel = FALSE, ...){
  # This function performs cross-validation to find the optimum value of tuning parameter "rho"
  # Args:
  #   x: A n-by-p training data matrix; n observations and p features. Used to train the classifier.
  #   y: A numeric vector of class labels of length n: 1, 2, ...., K if there are K classes.
  #      Each element of y corresponds to a row of x; i.e. these are the class labels for the observations in x.
  #   rhos: A vector of tuning parameters that control the amount of soft thresholding performed. If "rhos" is
  #         provided then a number of models will be fit (one for each element of "rhos"), and a number of
  #         predicted class labels will be output (one for each element of "rhos").
  #   beta: A smoothing term. A Gamma(beta,beta) prior is used to fit the Poisson model.
  #         Recommendation is to just leave it at 1, the default value.
  #   nfolds: The number of folds in the cross-validation; default is 5-fold cross-validation.
  #   type: How should the observations be normalized within the Poisson model, i.e. how should the size factors be estimated?
  #         Options are "quantile" or "deseq" (more robust) or "mle" (less robust).
  #         In greater detail: "quantile" is quantile normalization approach of Bullard et al 2010 BMC Bioinformatics
  #         "deseq" is median of the ratio of an observation to a pseudoreference obtained by taking the geometric mean,
  #         described in Anders and Huber 2010 Genome Biology and implemented in Bioconductor package "DESeq", and "mle" is the
  #         sum of counts for each sample; this is the maximum likelihood estimate under a simple Poisson model.
  #   folds: Instead of specifying the number of folds in cross-validation, one can explicitly specify the folds. To do this,
  #          input a list of length r (to perform r-fold cross-validation). The rth element of the list should be a vector
  #          containing the indices of the test observations in the rth fold.
  #   transform: Should data matrices x and xte first be power transformed so that it more closely fits the Poisson model?
  #              TRUE or FALSE. Power transformation is especially useful if the data are overdispersed relative to the Poisson model.
  #   alpha: If transform=TRUE, this determines the power to which the data matrices x and xte are transformed. If alpha = NULL then
  #          the transformation that makes the Poisson model best fit the data matrix x is computed. (Note that alpha is computed
  #          based on x, not based on xte). Or a value of alpha, 0<alpha<=1, can be entered by the user.
  #   prior: Vector of length equal to the number of classes, representing prior probabilities for each class. If NULL then
  #          uniform priors are used (i.e. each class is equally likely).

  # ii <- NULL
  if (!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
  if (transform && is.null(alpha)) alpha <- FindBestTransform(x)

  ## Apply power transform here.
  if (transform){
    if (alpha <= 0 || alpha > 1){
      stop("alpha must be between 0 and 1")
    }
    x <- x^alpha
  }

  if (is.null(rhos)){
    ns <- NullModel(x, type = type)$n  ## Normalized counts.
    uniq <- sort(unique(y))
    maxrho <- rep(NA, length(uniq))

    for (k in 1:length(uniq)){
      a <- colSums(x[y == uniq[k], ]) + beta
      b <- colSums(ns[y == uniq[k], ]) + beta
      maxrho[k] <- max(abs(a/b - 1)*sqrt(b), na.rm = TRUE)
    }

    rhos <- seq(0, max(maxrho, na.rm = TRUE)*(2/3), len = tuneLength)   ## len will be included in the control arguement as "tuneLength"
  }

  ### f-fold r-repeats fold indices.
  if (is.null(folds)){
    foldIdx <- lapply(1:repeats, function(u){
      tmp.fold <- balanced.folds(y, nfolds = nfolds)
      names(tmp.fold) <- paste("Fold.", 1:length(tmp.fold), sep = "")
      tmp.fold
    })

    names(foldIdx) <- paste("Repeat.", 1:repeats, sep = "")
  } else {
    foldIdx <- folds
  }

  ### FOREACH ILE PARALEL KODLAMA:
  # c("rhos", "x", "y", "beta", "prior", "Classify",
  #   "NullModel", "NullModelTest", "GetD", "Soft", "foldIdx")
  if (parallel){
    nn <- length(foldIdx)
    foldResults <- foreach(ii = 1:nn, .inorder = TRUE, .multicombine = TRUE,
                           .export = c("Classify", "NullModel", "NullModelTest", "GetD", "Soft")) %dopar% {
     f <- foldIdx[[ii]]  # ii.th repeat
     nfolds <- length(f)  # number of folds in ii.th repeat
     errs <- nnonzero <- accuracy <- matrix(NA, nrow = nfolds, ncol = length(rhos))

     confusionTables.fold <- confusionTables.fold.rho <- list()
     for (i in 1:nfolds){  ### loop for each fold within ii'th repeat
       tr <- -f[[i]]
       te <- f[[i]]

       # Type olarak neden sadece "quantile" üzerinden classify edilmiş??
       out <- Classify(x[tr, ], y[tr], x[te, ], rhos = rhos, beta = beta,
                       type = "quantile", prior = prior, transform = FALSE) # Have already power-transformed x, so don't need to do it again!!!

       for (j in 1:length(rhos)){
         errs[i, j] <- sum(out[[j]]$ytehat != y[te])
         nnonzero[i, j] <- sum(colSums(out[[j]]$ds != 1) != 0)
         accuracy[i, j] <- (length(y[te]) - sum(out[[j]]$ytehat != y[te])) / length(y[te])

         confusionTables.fold.rho[[j]] <- table(Prediction = out[[j]]$ytehat, Reference = y[te])
       }
       names(confusionTables.fold.rho) <- paste0("rho.", rhos)
       confusionTables.fold[[i]] <- confusionTables.fold.rho
     }
     names(confusionTables.fold) <- paste0("Fold.", 1:nfolds)

     return(list(errs = errs, nnonzero = nnonzero, accuracy = accuracy, rhos = rhos, confTables = confusionTables.fold))
   }
    names(foldResults) <- paste("Repeat.", 1:length(foldResults), sep = "")
  } else {
    nn <- length(foldIdx)
    foldResults <- foreach(ii = 1:nn, .inorder = TRUE, .multicombine = TRUE,
                           .export = c("Classify", "NullModel", "NullModelTest", "GetD", "Soft")) %do% {
       f <- foldIdx[[ii]]
       nfolds <- length(f)
       errs <- nnonzero <- accuracy <- matrix(NA, nrow = nfolds, ncol = length(rhos))

       confusionTables.fold <- confusionTables.fold.rho <- list()
       for (i in 1:nfolds){
         tr <- -f[[i]]
         te <- f[[i]]

         # Type olarak neden sadece "quantile" üzerinden classify edilmiş??
         out <- Classify(x[tr, ], y[tr], x[te, ], rhos = rhos, beta = beta,
                         type = "quantile", prior = prior, transform = FALSE) # Have already power-transformed x, so don't need to do it again!!!
         for (j in 1:length(rhos)){
           errs[i, j] <- sum(out[[j]]$ytehat != y[te])
           nnonzero[i, j] <- sum(colSums(out[[j]]$ds != 1) != 0)
           accuracy[i, j] <- (length(y[te]) - sum(out[[j]]$ytehat != y[te])) / length(y[te])

           confusionTables.fold.rho[[j]] <- table(Prediction = out[[j]]$ytehat, Reference = y[te])
         }
         confusionTables.fold[[i]] <- confusionTables.fold.rho
       }

       names(confusionTables.fold) <- paste0("Fold.", 1:nfolds)

       return(list(errs = errs, nnonzero = nnonzero, accuracy = accuracy, rhos = rhos, confTables = confusionTables.fold))
     }
    names(foldResults) <- paste("Repeat.", 1:length(foldResults), sep = "")
  }

  ### Accuracy and non-zero features avaraged over repeated folds.
  errs.average <- colMeans(plyr::ldply(lapply(foldResults, function(u){
    colMeans(u$errs)
  }), "rbind", .id = NULL))

  nnonzero.average <- colMeans(plyr::ldply(lapply(foldResults, function(u){
    colMeans(u$nnonzero)
  }), "rbind", .id = NULL))

  acc.average <- colMeans(plyr::ldply(lapply(foldResults, function(u){
    colMeans(u$accuracy)
  }), "rbind", .id = NULL))

  tuneResults <- data.frame(rho = rhos, Avg.Error = errs.average, Avg.NonZeroFeat. = nnonzero.average, Accuracy = acc.average)

  idx <- max(which(acc.average == max(acc.average)))
  bestrho <- rhos[idx]
  bestAccuracy <- tuneResults[idx, "Accuracy"]
  bestNonZeroFeat <- tuneResults[idx, "Avg.NonZeroFeat."]

  # tuneResults <- as.matrix(tuneResults)
  ## BURADA PLDA için yazılmış olan sınıfa göre bir S4 sonuç dönecek.
  ## discrete.train

  res <- list(inputs = list(x = x, y = y, x.test = NULL, y.test = NULL), control = list(foldIdx = foldIdx),
              tuningResults = list(rhos = rhos, bestrho = bestrho, bestrho.idx = idx, bestAccuracy = bestAccuracy,
                                   bestNonZeroFeat = bestNonZeroFeat, results = tuneResults, alpha.transform = alpha),
              selectedGenes = NULL,
              foldResults = foldResults)
  return(res)
}

#' @importFrom sSeq getT getAdjustDisp rowVars
trainNBLDA <- function(x, y, xte = NULL, beta = 1, type = c("mle", "deseq", "quantile", "none", "TMM"),
                       prior = NULL, truephi = NULL, ...){
  # Args:
  #   x: a data.frame or matrix with counts. Samples are in the rows and genes are in the columns.
  #   y: a vector (numeric or factor) including the class labels of each subject.
  #   xte: a data frame or matrix of test samples. If NULL, train set is used as test set.
  #   beta:
  #   type:
  #   prior:
  #   truephi:

  if (is.null(prior)){
    prior <- rep(1/length(unique(y)), length(unique(y)))
    ## prior vectorunun sinif sayisi ile esit uzunlukta olmasi lazım.
    ## Bu durumu kontrol edip gerekirse uyari verilecek.
  }

  null.out <- NullModel(x, type = type)
  ns <- null.out$n
  nste <- NullModelTest(null.out, x, xte, type = type)$nste

  uniq <- sort(unique(y))
  ds <- GetDnb(ns, x, y, beta)

  ## Dispersion estimates.
  if (!is.null(truephi)){
    if (length(truephi) >= 2 & (length(truephi) != ncol(ns))){
      truephi <- truephi[1]
      disperhat <- rep(truephi, ncol(ns))
      warning("The length of \"truephi\" should be the same as number of features. Only the first element is used and replicated for all features.")

    } else if (length(truephi) == 1){
      disperhat <- rep(truephi, ncol(ns))

    } else if (length(truephi) >= 2 & (length(truephi) == ncol(ns))){
      disperhat <- truephi
    }

  } else {
    ####### TO DO #######
    # tt değerinin 0'a eşit olması durumu ile analiz yapılıyor ama normalde üstteki değere shrink ediliyor. Yu et al çalışmasındaki gibi.
    # Bunun nedeni araştırılacak. 0 Değerini almak problem yaratır mı diye incelenecek.
    tt <- getT(t(x), sizeFactors = rep(1, nrow(x)), verbose = FALSE, propForSigma = c(0, 1))$target

    ### Moment estimation of gene-wise dispersions.
    rM = rowMeans(t(x))
    rV = rowVars(t(x))

    disp = (rV - rM)/rM^2   ## Dispersion estimates using method-of-moments.
    # disp0 = numeric()

    ## Negative dispersions are set to 0.
    disp0 <- sapply(disp, function(x){
      max(0, x)
    })

    disperhat <- getAdjustDisp(disp0, shrinkTarget = tt, verbose = FALSE)$adj
    disperhat <- pmax(rep(0, length(disperhat)), disperhat) ## Negative dispersions are set to 0.
  }

  # if (!is.null(truephi)){
  #   if (length(truephi) >= 2){
  #     truephi <- truephi[1]
  #     warning("\"truephi\" should be given as a single value. Only first element is used.")
  #   }
  #   disperhat <- rep(truephi, ncol(nste))
  #
  # } else {
  #   tt = getT(t(x), sizeFactors = rep(1, nrow(x)), verbose = FALSE)$target
  #
  #   ### Moment estimation of gene-wise dispersions.
  #   rM = rowMeans(t(x))
  #   rV = rowVars(t(x))
  #
  #   disp = (rV - rM)/rM^2
  #   disp0 = numeric()
  #
  #   for (i in 1:length(disp)){
  #     a1 = 0
  #     a2 = disp[i]
  #     disp0[i] = max(a1, a2)  ## Negative dispersions are set to 0.
  #   }
  #   names(disp0) <- names(disp)
  #
  #   disperhat = getAdjustDisp(disp0, shrinkTarget = tt, verbose = FALSE)$adj
  # }

  phihat <- as.numeric(disperhat)
  discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

  # Replace Infinity with zero dispersions.
  inv.phihat <- (1/phihat)
  if (any(inv.phihat == Inf | inv.phihat == -Inf)){
    id.inf <- which(inv.phihat == Inf | inv.phihat == -Inf)
    inv.phihat[id.inf] <- 0
  } else {
    id.inf <- NULL
  }

  for (k in 1:length(uniq)){
    for (l in 1:nrow(xte)){
      dstar = ds[k, ]
      part2 = 1 + nste[l, ] * dstar * phihat
      part1 = dstar/part2

      part3 <- (1/phihat) * log(part2)

      if (!is.null(id.inf)){
        part3.limits <- nste[l, ] * dstar
        part3[id.inf] <- part3.limits[id.inf]
      }

      discriminant[l, k]<- sum(xte[l, ] * log(part1)) - sum(part3) + log(prior[k])
    }
  }

  preds <- uniq[apply(discriminant, 1, which.max)]  ## Predicted class labels.

  result <- structure(list(ns = ns, nste = nste, ds = ds, discriminant = discriminant,
                           ytehat = preds, x = x, y = y, xte = xte, type = type,
                           disp = disperhat, null.out = null.out),
                      class = "list")

  return(result)
}

# #
# data <- data.trainS4
# method = "NBLDA"       # c("PLDA", "PLDA2", "NBLDA")
# normalize = "deseq"        # c("deseq", "TMM", "none")
# ref = "T"
# class.labels = "condition"
# control = discreteControl(number = 5, repeats = 1, tuneLength = 10, parallel = TRUE)
# type <- normalize

## Discrete classifiers: PLDA, PLDA2, NBLDA
## data: a DESeqDataSet object.
classify.discrete <- function(data, method = c("PLDA", "PLDA2", "NBLDA"), normalize = c("deseq", "TMM", "none"),
                              ref = NULL, class.labels = NULL, control = discreteControl(), ...){

  call <- match.call()
  arg.list <- list(data = data, method = method, normalize = normalize,
                   ref = ref, class.labels = class.labels, control = control, ...)

  # extract arguements from control list.
  rho <- control$rho
  rhos <- control$rhos
  beta <- control$beta
  prior <- control$prior
  alpha <- control$alpha
  folds <- control$foldIdx
  tuneLength <- control$tuneLength

  # size factor estimation method. should be one of "deseq", "TMM" and "none"
  type <- match.arg(normalize)
  method <- match.arg(method)

  transform <- FALSE
  if (method == "PLDA2") transform <- TRUE   ### PLDA2 is used fr power transformed PLDA.

  # Raw counts from DESeqDataSet object. Samples in the rows and genes in the columns.
  x <- t(counts(data))

  # Class labels.
  y <- colData(data)[ ,class.labels]

  if (method %in% c("PLDA", "PLDA2")){
    model.fit <- trainPLDA(x, y, rhos = rhos, beta = beta, nfolds = control$number, repeats = control$repeats, type = type,
                           folds = folds, transform = transform, alpha = alpha, prior = prior, tuneLength = tuneLength, parallel = control$parallel)

    # # CONFUSION MATRIX MANTIGI CARET ILE BENZER OLACAK
    model.fit.final <- if (transform){
      Classify(x, y, xte = x, rho = model.fit$tuningResults$bestrho, beta = beta, rhos = NULL,
               type = type, prior = prior, transform = TRUE,
               alpha = model.fit$tuningResults$alpha.transform)
    } else {
      Classify(x, y, xte = x, rho = model.fit$tuningResults$bestrho, beta = beta, rhos = NULL,
               type = type, prior = prior, transform = FALSE)
    }

    # export "control" list into plda train object.
    if (is.null(control$foldIdx)){
      fold.idx <- model.fit$control$foldIdx
      model.fit$control <- control
      model.fit$control$foldIdx <- fold.idx
    } else {
      model.fit$control <- control
    }

    # bestrho kullanılarak train set yeniden train edilerek modele seçilen genlerin isimleri belirlenecek.
    # Bu genler selectedFeatures olarak train objesi içerisinde verilecek.
    ds <- model.fit.final$ds
    selectedGenesIdx <- which(!apply(ds == 1, 2, all))
    selectedGenesNames <- colnames(x)[selectedGenesIdx]
    if (model.fit$tuningResults$bestrho != 0){
      model.fit.final$SelectedGenes <- list(selectedGenesIdx = selectedGenesIdx, selectedGenesNames = selectedGenesNames)
    } else {
      model.fit.final$SelectedGenes <- list(selectedGenesIdx = NULL, selectedGenesNames = NULL)
    }

    ### Confusion matrix daha sonra overall predicted values üzerinden düzenlenecek. Aggregated over folds.
    ### Confusion matrix "caret" mantığı ile veriliyor. Bütün repeat ve fold'lara ait classification table'lar
    ### kullanılarak bir aggregated sonuç verilecek.

    bestrho.index <- model.fit$tuningResults$bestrho.idx

    aggregatedTable <- unlist(lapply(model.fit$foldResults, function(x){
      tmp <- x$confTables
      lapply(tmp, function(y){
        y[[bestrho.index]]
      })
    }), recursive = FALSE)

    tmp <- as.matrix(aggregatedTable[[1]])
    for (i in 2:length(aggregatedTable)){
      tmp <- tmp + as.matrix(aggregatedTable[[i]])
    }

    ## BURADA ELDE EDILEN DEGERLER TAMSAYI OLMADIGI DURUMDA HATA VERIYOR. BU SEBEPLE OVERALL DEGERLERIN
    ## KULLANILMASI DAHA UYGUNDUR. ANCAK ELDE EDILEN SONUCLARDAKI P-DEGERLERI VE GUVEN ARALIKLARININ
    ## DIKKATE ALINMAMASI GEREKIYOR. BU BILGININ BIR NOT OLARAK VERILMESI UYGUN OLABILIR.
    tmp <- tmp / model.fit$control$repeats

    tmp.attr <- attributes(tmp)
    names(tmp.attr$dimnames) <- c("Prediction", "Reference")
    attributes(tmp) <- tmp.attr

    confM <- confusionMatrix(round(tmp, 0), positive = ref)
    confM$tableRounded <- round(tmp, 2)
    confM$table <- round(tmp, 0)

    attr(confM, "class") <- "confMat"

  } else if (method == "NBLDA"){

    ### Calculate train set performances using repeated cross-validation:
    folds <- foldIndex(n = nrow(colData(data)), nFolds = control$number, repeats = control$repeats)
    parallel <- control$parallel

    if (parallel){
      nn <- length(folds$indexIn)
      fold.tmp <- folds$indexIn
      trainResults.AllFolds <- foreach(i = 1:nn, .inorder = FALSE, .multicombine = TRUE,
                                       .export = c("NullModel", "NullModelTest", "GetDnb", "Soft",
                                                  "trainNBLDA", "x", "y", "fold.tmp", "type")) %dopar% {
        idx <- fold.tmp[[i]]
        input.foldIn <- x[idx, ]
        output.foldIn <- y[idx]

        # Samples for validation in this fold.
        input.foldOut <- x[-idx, ]
        output.foldOut <- y[-idx]

        trainResults.fold <- trainNBLDA(x = input.foldIn, y = output.foldIn, xte = input.foldOut, type = type)

        ### Predictions for train set:
        pred.fold <- factor(trainResults.fold$ytehat, levels = levels(y))
        tbl <- table(Prediction = pred.fold, Reference = output.foldOut)
        acc.fold <- sum(diag(tbl))/sum(tbl)

        return(list(accuracies = acc.fold, predictions = pred.fold, confTables = tbl))
      }
      names(trainResults.AllFolds) <- names(folds$indexIn)
    } else {
      trainResults.AllFolds <- lapply(folds$indexIn, function(idx){
        input.foldIn <- x[idx, ]
        output.foldIn <- y[idx]

        # Samples for validation in this fold.
        input.foldOut <- x[-idx, ]
        output.foldOut <- y[-idx]

        trainResults.fold <- trainNBLDA(x = input.foldIn, y = output.foldIn, xte = input.foldOut, type = type)

        ### Predictions for train set:
        pred.fold <- factor(trainResults.fold$ytehat, levels = levels(y))
        tbl <- table(Prediction = pred.fold, Reference = output.foldOut)
        acc.fold <- sum(diag(tbl))/sum(tbl)

        return(list(accuracies = acc.fold, predictions = pred.fold, confTables = tbl))
      })
    }

    acc <- ldply(lapply(trainResults.AllFolds, function(x)x$accuracies), rbind, .id = "Folds")
    colnames(acc)[2] <- "Accuracy"

    overallAccuracy <- mean(acc[ ,"Accuracy"])

    ## Final Model
    trainedOverall <- trainNBLDA(x = t(counts(data)), y = colData(data)[ ,class.labels], xte = x, type = normalize)
    model.fit <- trainedOverall

    model.fit[["tuningResults"]] <- list()
    model.fit.final <- model.fit
    model.fit.final[["overallAccuracy"]] <- overallAccuracy
    model.fit.final[["trainParameters"]] <- model.fit$null.out

    control$foldIdx <- folds
    model.fit[["control"]] <- control
    model.fit[["inputs"]] <- list(x = model.fit$x, y = model.fit$y)

    model.fit[["foldInfo"]] <- list(foldAccuracies = acc,
                                    foldPredictions = lapply(trainResults.AllFolds, function(x)x$predictions),
                                    foldIndex = folds)

    # Confusion matrix over folds and repeats
    aggregatedTable <- lapply(trainResults.AllFolds, function(x){
      x$confTables
    })

    tmp <- as.matrix(aggregatedTable[[1]])
    for (i in 2:length(aggregatedTable)){
      tmp <- tmp + as.matrix(aggregatedTable[[i]])
    }

    tmp <- tmp / control$repeats

    confM <- confusionMatrix(round(tmp, 0), positive = ref)
    confM$tableRounded <- round(tmp, 2)
    confM$table <- round(tmp, 0)
    attr(confM, "class") <- "confMat"
  }


  # ## discrete.train adı ile bir S4 class yazılacak.
  trainedModel <- new("discrete.train",
                      inputs = model.fit$inputs,
                      control = model.fit$control,
                      crossValidatedModel = model.fit,
                      tuningResults = model.fit$tuningResults,
                      finalModel = model.fit.final,
                      callInfo = list(method = method, normalize = normalize,
                                      class.labels = class.labels, class.names = as.character(unique(colData(data)[ ,class.labels])),
                                      ref = ref))

  ## All the information about classification model.
  modelInfo <- new("MLSeqModelInfo",
                   method = method,
                   transformation = NA_character_,
                   normalization = normalize,
                   preProcessing = NA_character_,
                   ref = ref,
                   control = trainedModel@control,   ## PLDA'da kullanılan discreteControl(...) objesi
                   confusionMat = confM,   ## Cross-validated confusion matrix
                   trainedModel = trainedModel,
                   trainParameters = trainedModel@finalModel$trainParameters,
                   call = list(call = call, args = arg.list))

  if (transform){
    transformedData <- model.fit.final$x
  } else {
    transformedData <- NULL
  }

  result = new("MLSeq",
               modelInfo = modelInfo,
               inputObject = list(rawData = x, transformedData = transformedData, conditions = y))

  ## Metadata for MLSeq object.
  meta <- new("MLSeqMetaData", classLabel = class.labels, rawData.DESeqDataSet = data)
  result@metaData <- meta

  return(result)
}
