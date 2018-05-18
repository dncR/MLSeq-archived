## voom ve NSC turu sonuclarda predict edebilmesi icin MLSeq sinifi metaData icerisine bir type argumani eklenecek
## Bu argumana göre predict farkli sekillerde davranacak.

#' @importFrom caret predict.train
#' @importFrom stats median
#' @importFrom SummarizedExperiment assay mcols 'mcols<-'
predict.MLSeq <- function (object, test.data, ...){
  # Args:
  #   object: a model of "MLSeq" class returned by classify(...)
  #   test.data: a DESeqDataSet instance of new observations.

  test.pred <- NULL
  model <- object
  if (class(model) != "MLSeq"){
    stop("'model' arguement should be an \"MLSeq\" object.")
  }

  if (!("DESeqDataSet" %in% class(test.data))){
    stop("'test.data' arguement should be a \"DESeqDataSet\" object.")
  }

  ## Prediction steps for "caret" based classifications.
  if (class(trained(model)) == "train"){
    if (normalization(model) == "deseq" & transformation(model) == "vst"){  # "deseq-vst"
      counts.test <- counts(test.data)
      counts.train <- counts(model@inputObject$rawData)    ### Raw Counts icin getter yazilinca burada kullanilacak.  rawCounts(...) adinda bir getter kullanilabilir.

      #Calculation of test set size factors using geometric means from training data
      #Genes in the row, samples in the column
      geomts = counts.test / exp(rowMeans(log(counts.train)))   ## Geometric mean of test data using estimators from train data
      sizeF.ts = apply(geomts, 2, function(x) median(x))

      test.dataSF <- estimateSizeFactors(test.data)    # Estimate Size factors:
      sizeFactors(test.dataSF) <- sizeF.ts             # Replace size factors with size factors which are estimates using train set parameters.

      ## Change dispersion function of test data with dispersion function of train data
      dispersionFunction(test.dataSF) <- trainParameters(model)$disperFunc

      transformedData <- varianceStabilizingTransformation(test.dataSF, fitType = "local", blind = FALSE)

      dataVST <- t(as.matrix(assay(transformedData)))
      input <- dataVST   ## Input data from transformed expression data.

    } else if (normalization(model) == "deseq" & transformation(model) == "rlog"){   # "deseq-rlog"
      counts.test <- counts(test.data)
      counts.train <- counts(model@inputObject$rawData)    ### Raw Counts icin setter yazilinca burada kullanilacak.  rawCounts(...)

      #Calculation of test set size factors using geometric means from training data
      #Genes in the row, samples in the column
      geomts = counts.test / exp(rowMeans(log(counts.train)))   ## Geometric mean of test data using estimators from train data
      sizeF.ts = apply(geomts, 2, function(x) median(x))

      test.dataSF <- estimateSizeFactors(test.data)    # Estimate Size factors:
      sizeFactors(test.dataSF) <- sizeF.ts

      test.dataDisp <- estimateDispersions(test.dataSF, fitType = "local")

      ## Required train parameters for test data
      mcols(test.dataDisp)$dispFit <- trainParameters(model)$dispFit
      betaPrior <- trainParameters(model)$betaPrior
      intercept <- trainParameters(model)$intercept

      test.dataRLOG <- rlog(test.dataDisp, blind = FALSE,
                            intercept = intercept, betaPriorVar = betaPrior)

      dataRLOG = t(assay(test.dataRLOG))
      input <- dataRLOG   ## Input data from transformed expression data.

    } else if (normalization(model) == "deseq" & transformation(model) == "logcpm"){  ## deseq-logcpm

      rawCounts = counts(test.data, normalized = FALSE)
      countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
      counts.train <- counts(model@metaData$metaData@rawData.DESeqDataSet)    ### Raw Counts icin setter yazilinca burada kullanilacak.  rawCounts(...)

      #Calculation of test set size factors using geometric means from training data
      #Genes in the row, samples in the column
      geomts = rawCounts / exp(rowMeans(log(counts.train)))   ## Geometric mean of test data using estimators from train data
      sizeF.ts = apply(geomts, 2, function(x) median(x))

      RLE <- sizeF.ts / colSums(rawCounts)
      normFactors <- RLE / exp(mean(log(RLE)))

      countsDGE.normalized <- calcNormFactors(countsDGE, method = "RLE")   ## RLE: DESeq mantigi ile normalize ediyor.
      countsDGE.normalized$samples$norm.factors <- normFactors   ### Replace normFactors using normFactors obtained from train set parameters.
      countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek

      input <- t(countsDGE.transformed)   ## Input data from transformed expression data.

    } else if (normalization(model) == "none" & transformation(model) == "logcpm"){
      rawCounts = counts(test.data, normalized = FALSE)
      countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
      countsDGE.transformed <- cpm(countsDGE, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek

      input <- t(countsDGE.transformed)   ## Input data from transformed expression data.

    } else if (normalization(model) == "tmm" & transformation(model) == "logcpm"){

      rawCounts = counts(test.data, normalized = FALSE)
      referenceSample <- trainParameters(model)$refSample

      ## Chech if the feature names of reference sample are in the same order with rawCounts.
      if (identical(rownames(rawCounts), names(referenceSample))){
        rawCounts <- cbind(rawCounts, referenceSample)
      } else {
        stop(warning("Reference sample either does not have same features or the features are not in the same order as test set. Calculation stops.."))
      }

      countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))

      ## Reference column is selected from train set.
      countsDGE.normalized <- calcNormFactors(countsDGE, method = "TMM", refColumn = ncol(countsDGE))   ## RLE: DESeq mantigi ile normalize ediyor.
      countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek.

      input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
      input <- input[-nrow(input), ]   ## Remove reference sample from test data.
    }

    test.pred = predict(trained(model), input, ...)    ## Predicted classes.

  } else if (class(trained(model)) == "voom.train"){
    if (method(model) %in% c("voomDLDA", "voomDQDA")){
      test.pred <- predict.voomDDA(object = trained(model)@finalModel$model, newdata = counts(test.data), ...)
      test.pred <- as.factor(test.pred)
    } else {
      test.pred <- predict.voomNSC(fit = trained(model)@finalModel$model, newdata = counts(test.data), ...)
      test.pred <- as.factor(test.pred)
    }
  } else if (class(trained(model)) == "discrete.train"){
    ## Predictions for PLDA and PLDA2
    if (method(model) %in% c("PLDA", "PLDA2")){
      args <- list(x = input(model)$rawData,
                   y = input(model)$conditions,
                   xte = t(counts(test.data)),
                   rho = model@modelInfo@trainedModel@finalModel$rho,
                   beta = control(model)$beta,
                   type <- normalization(model),
                   prior <- control(model)$prior,
                   transform <- ifelse(method(model) == "PLDA", FALSE, TRUE),
                   alpha <- model@modelInfo@trainedModel@finalModel$alpha,
                   null.out <- model@modelInfo@trainParameters,
                   ds = model@modelInfo@trainedModel@finalModel$ds)

      test.pred <- do.call("predictPLDA", args)

    } else if (method(model) == "NBLDA"){
      args <- list(x = input(model)$rawData,
                   y = input(model)$conditions,
                   xte = t(counts(test.data)),
                   beta = control(model)$beta,
                   type = normalization(model),
                   prior = control(model)$prior,
                   truephi = control(model)$truephi,
                   null.out = model@modelInfo@trainedModel@finalModel$trainParameters,
                   ds = model@modelInfo@trainedModel@finalModel$ds,
                   disperhat = model@modelInfo@trainedModel@finalModel$disp)

      test.pred <- do.call("predictNBLDA", args)
    }
  }

  # res <- list(MLSeqObject = model, Predictions = test.pred)
  return(test.pred)
}



## Previous function predictClassify
#' @rdname predict
#' @export
predictClassify <- predict.MLSeq
