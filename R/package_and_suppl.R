###### Documentation for package and supplementary files ######

#' @title Machine learning interface for RNA-Seq data
#'
#' @description This package applies machine learning methods, such as Support Vector Machines (SVM), Random Forest (RF),
#' Classification and Regression Trees (CART), Linear Discriminant Analysis (LDA) and more to RNA-Seq data. MLSeq combines
#' well-known differential expression algorithms from bioconductor packages with functions from a famous package \code{caret},
#' which has comprehensive machine learning algorithms for classification and regression tasks. Although \code{caret} has 200+
#' classification/regression algorithm built-in, approximately 85 classification algorithms are used in \code{MLSeq} for classifying
#' gene-expression data. See \code{availableMethods()} for further information.
#'
#' @seealso \code{\link{availableMethods}}, \code{\link[caret:modelLookup]{getModelInfo}}
#'
#' \tabular{ll}{
#'   Package: \tab MLSeq\cr
#'   Type: \tab Package\cr
#'   License: \tab GPL (>= 2)\cr
#' }
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Bernd Klaus, Ahmet Ozturk and Ahmet Ergun Karaagaoglu
#'
#' -----------------
#'
#' Maintainers:
#'
#' Gokmen Zararsiz, \email{gokmenzararsiz@erciyes.edu.tr}
#'
#' Dincer Goksuluk \email{dincer.goksuluk@hacettepe.edu.tr}
#'
#' Selcuk Korkmaz \email{selcukorkmaz@hotmail.com}
#'
#' @docType package
#' @name MLSeq-package
#' @rdname MLSeq-package
#'
#' @keywords package
NULL


#' @title Cervical cancer data
#'
#' @description Cervical cancer data measures the expressions of 714 miRNAs of human samples. There are 29 tumor and 29 non-tumor cervical
#' samples and these two groups are treated as two separete classes.
#'
#' @format A data frame with 58 observations on the following 715 variables.
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2880020/#supplementary-material-sec}
#'
#' @references
#' Witten, D., et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical biomarker
#' analysis in a collection of cervical tumours and matched controls. BMC Biology, 8:58
#'
#' @docType data
#' @name cervical
#' @rdname cervical
#'
#' @examples
#' \dontrun{
#' data(cervical)
#' }

NULL
