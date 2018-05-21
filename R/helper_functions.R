#' @title Available classification/regression methods in \code{MLSeq}
#'
#' @description This function returns a character vector of available classification/regression methods in \code{MLSeq}. These methods
#' are imported from \code{caret} package. See details below.
#'
#' @details There are 200+ methods available in \code{caret}. We import approximately 85 methods which are available for "classification" task.
#' Some of these methods are available for both classification and regression tasks. \code{availableMethods()} returns a character vector
#' of available methods in \code{MLSeq}. These names are directly used in \code{\link{classify}} function with arguement \code{method}.
#' See \url{http://topepo.github.io/caret/available-models.html} for a complete list of available methods in \code{caret}.
#' Run \code{printAvailableMethods()} to print detailed information about classification methods (prints to R Console).
#'
#' @note Available methods in \code{MLSeq} will be regularly updated. Some of the methods might be removed as well as some others
#' took its place in \code{MLSeq}. Please check the available methods before fitting the model. This function is inspired
#' from the function \code{getModelInfo()} in \code{caret} and some of the code chunks and help texts are used here.
#'
#' @param model a character string indicating the name of classification model. If NULL, all the available methods from \code{MLSeq}
#' is returned. Otherwise, the methods which are complete or partial matches to requested string is returned. See \code{regex} for
#' details.
#' @param regex a logical: should a regular expressions be used? If FALSE, a simple match is conducted against the whole name of the model.
#' @param \dots options to pass to \code{\link[base:grep]{grepl}}.
#'
#' @return a requested or complete character vector of available methods.
#'
#' @seealso \code{\link{classify}}, \code{\link[caret:modelLookup]{getModelInfo}}, \code{\link[caret]{train}}
#'
#' @name Available-classifiers
#' @rdname Available-classifiers
#' @aliases Available-classifiers printAvailableMethods
#'
#' @export
availableMethods <- function(model = NULL, regex = TRUE, ...){
  # This function returns the available classification methods in MLSeq. Name of the methods are identical to
  # methods in caret package.
  # Args:
  #   model: a character string indicating the name of classification model.
  #   regex: a logical. should a regular expressions be used? If FALSE, a simple match is conducted
  #          against the whole name of the model.
  #   ...: options to pass to grepl(...)

  loaded <- try({load(system.file("extdata", "availableModels.RData", package = "MLSeq"))})   # Use this for Bioconductor version

  # Local files:
  # Uncomment this line for loading dataset from source folder rather than installed package.
  # This is useful when debugging source files.
  # loaded <- try({load(paste(getwd(), "/inst/extdata/availableModels.RData", sep = ""))})  # iMAC

  if (class(loaded) == "try-error"){
    available_models <- NULL
    rm(loaded)
  } else {
    rm(loaded)
  }

  if (!is.null(available_models)){
    if (!is.null(model)){
      keepers <- if (regex){
        grepl(model, available_models, ...)
      } else {
        which(model == available_models)
      }
      models <- available_models[keepers]

      if (length(models) == 0) {
        stop("This model is not available in MLSeq.")
      } else {
        return(models)
      }
    } else {
      return(available_models)
    }

  } else {
    warning("Available models may not be loaded properly from \"MLSeq\" package.")
    return(NULL)
  }
}


#' @rdname Available-classifiers
#' @export
printAvailableMethods <- function(){
cat("
--------------------|------------------------------------------------------------------------------------
  Abbrev.           |   Description (B.C.: Binary Classification)
--------------------|------------------------------------------------------------------------------------
  amdai             |   Adaptive Mixture Discriminant Analysis
  AdaBag            |   Bagged AdaBoost
  treebag           |   Bagged CART
  bagFDA            |   Bagged Flexible Discriminant Analysis
  bayesglm          |   Bayesian Generalized Linear Model
  gamboost          |   Boosted Generalized Additive Model
  glmboost          |   Boosted Generalized Linear Model
  BstLm             |   Boosted Linear Model
  LogitBoost        |   Boosted Logistic Regression
  bstSm             |   Boosted Smoothing Spline
  blackboost        |   Boosted Tree  (B.C.)
  bstTree           |   Boosted Tree
  C5.0              |   C5.0
  rpart             |   CART
  rpart1SE          |   CART
  rpart2            |   CART
  rpartScore        |   CART or Ordinal Responses
  cforest           |   Conditional Inference Random Forest
  ctree             |   Conditional Inference Tree
  ctree2            |   Conditional Inference Tree
  C5.0Cost          |   Cost-Sensitive C5.0
  rpartCost         |   Cost-Sensitive CART
  deepboost         |   DeepBoost
  dda               |   Diagonal Discriminant Analysis
  dwdPoly           |   Distance Weighted Discrimination with Polynomial Kernel (B.C.)
  dwdRadial         |   Distance Weighted Discrimination with Radial Basis Function Kernel (B.C.)
  fda               |   Flexible Discriminant Analysis
  gam               |   Generalized Additive Model using Splines
  glm               |   Generalized Linear Model (B.C.)
  gpls              |   Generalized Partial Least Squares
  glmnet            |   glmnet
  protoclass        |   Greedy Prototype Selection
  hda               |   Heteroscedastic Discriminant Analysis
  hdda              |   High Dimensional Discriminant Analysis
  hdrda             |   High-Dimensional Regularized Discriminant Analysis
  kknn              |   k-Nearest Neighbors
  knn               |   k-Nearest Neighbors
  svmLinearWeights2 |   L2 Regularized Linear Support Vector Machines with Class Weights
  svmLinear3        |   L2 Regularized Support Vector Machine (dual) with Linear Kernel
  lvq               |   Learning Vector Quantization
  lda               |   Linear Discriminant Analysis
  lda2              |   Linear Discriminant Analysis
  stepLDA           |   Linear Discriminant Analysis with Stepwise Feature Selection
  dwdLinear         |   Linear Distance Weighted Discrimination (B.C.)
  loclda            |   Localized Linear Discriminant Analysis
  Mlda              |   Maximum Uncertainty Linear Discriminant Analysis
  mda               |   Mixture Discriminant Analysis
  avNNet            |   Model Averaged Neural Network
  mlp               |   Multi-Layer Perceptron
  mlpWeightDecay    |   Multi-Layer Perceptron
  mlpWeightDecayML  |   Multi-Layer Perceptron, multiple layers
  mlpML             |   Multi-Layer Perceptron, with multiple layers
  earth             |   Multivariate Adaptive Regression Spline
  gcvEarth          |   Multivariate Adaptive Regression Splines
  nb                |   Naive Bayes
  nbDiscrete        |   Naive Bayes
  pam               |   Nearest Shrunken Centroids
  nnet              |   Neural Network
  pcaNNet           |   Neural Networks with Feature Extraction
  ORFlog            |   Oblique Random Forest
  ORFpls            |   Oblique Random Forest (Partial least square as node model)
  ORFridge          |   Oblique Random Forest (Ridge regression as node model)
  ORFsvm            |   Oblique Random Forest (linear SVM as node model)
  pls               |   Partial Least Squares
  pda               |   Penalized Discriminant Analysis
  PenalizedLDA      |   Penalized Linear Discriminant Analysis
  plr               |   Penalized Logistic Regression
  multinom          |   Penalized Multinomial Regression
  qda               |   Quadratic Discriminant Analysis
  stepQDA           |   Quadratic Discriminant Analysis with Stepwise Feature Selection
  rbf               |   Radial Basis Function Network
  rf                |   Random Forest
  rda               |   Regularized Discriminant Analysis
  rlda              |   Regularized Linear Discriminant Analysis
  RRF               |   Regularized Random Forest
  Linda             |   Robust Linear Discriminant Analysis
  rmda              |   Robust Mixture Discriminant Analysis
  QdaCov            |   Robust Quadratic Discriminant Analysis
  rrlda             |   Robust Regularized Linear Discriminant Analysis
  bdk               |   Self-Organizing Map
  sdwd              |   Sparse Distance Weighted Discrimination
  sparseLDA         |   Sparse Linear Discriminant Analysis
  spls              |   Sparse Partial Least Squares
  gbm               |   Stochastic Gradient Boosting
  svmLinear         |   Support Vector Machines with Linear Kernel
  svmPoly           |   Support Vector Machines with Polynomial Kernel
  svmRadial         |   Support Vector Machines with Radial Basis Function Kernel
  voomDLDA          |   Voom-based Diagonal Linear Discriminant Analysis
  voomDQDA          |   Voom-based Diagonal Quadratic Discriminant Analysis
  voomNSC           |   Voom-based Nearest Shrunken Centroids
  PLDA              |   Poisson Linear Discriminant Analysis
  PLDA2             |   Poisson Linear Discriminant Analysis with Power Transformation
  NBLDA             |   Negative Binomial Linear Discriminant Analysis
---------------------------------------------------------------------------------------------------------
")
}


# This function gives the name of the model from caret's abbreviation.
# Burada verilen tanimlamalar bir dosya ile availableMethods() fonksiyonuna benzer sekilde kaydedilebilir.
modelDescription <- function(abbrv = NULL){

  abbreviations <- availableMethods()

  loaded <- try({load(system.file("extdata", "modelDescriptions.RData", package = "MLSeq"))})   # Use this for Bioconductor version

  # Local files:
  # Uncomment this line for loading dataset from source folder rather than installed package.
  # This is useful when debugging source files.
  # loaded <- try({load(paste(getwd(), "/inst/extdata/modelDescriptions.RData", sep = ""))})  # iMAC

  if (class(loaded) == "try-error"){
    descriptions <- NULL
    rm(loaded)
  } else {
    rm(loaded)
  }

  return(descriptions[which(abbrv == abbreviations)])
}


#' @importFrom caret createFolds
foldIndex <- function(data = NULL, n = NULL, nFolds = 5, repeats = 2){

  if (!is.null(data)){
    n = nrow(data)
  }

  indIn <- indOut <- list()

  for (j in 1:repeats){
    tmpIn = createFolds(y = 1:n, k = nFolds, list = TRUE, returnTrain = TRUE)
    tmpOut = lapply(tmpIn, function(x)c(1:n)[-x])

    indIn = c(indIn, tmpIn)
    indOut = c(indOut, tmpOut)
  }

  nms = paste(rep(paste("Fold", 1:nFolds, sep = ""), repeats),
              rep(paste(".Rep", 1:repeats, sep = ""), c(rep(nFolds, repeats))), sep = "")
  names(indIn) <- names(indOut) <- nms
  return(list(indexIn = indIn, indexOut = indOut))
}
