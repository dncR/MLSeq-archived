######### COMPARISONS ##########
library(DESeq2)
library(edgeR)
library(doParallel)
library(caret)
library(sSeq)
library(MLSeq)

source("R/all_classes.R")
source("R/all_generics.R")
source("R/voomFunctions.R")
source("R/classify.R")
source("R/helper_functions.R")
source("R/package_and_suppl.R")
source("R/predict.R")
source("R/PLDA_NBLDA_Functions.R")
source("R/methods.R")

# # PoiClaClu
# source("R/R-PoiClaClu/PoiClaClu_Classify.cv.R")
# source("R/R-PoiClaClu/PoiClaClu_Classify.R")
# source("R/R-PoiClaClu/PoiClaClu_CountDataSet.R")
# source("R/R-PoiClaClu/PoiClaClu_FindBestTransform.R")
# source("R/R-PoiClaClu/PoiClaClu_Functions.R")
# source("R/R-PoiClaClu/PoiClaClu_NullModel.R")
# source("R/R-PoiClaClu/PoiClaClu_NullModelTest.R")
# source("R/R-PoiClaClu/PoiClaClu_PoissonDistance.R")
#


registerDoParallel(4)

data(cervical)

# a subset of cervical data with first 150 features.
data <- cervical[c(1:250), ]
# data <- cervical

# defining sample classes.
class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))

n <- ncol(data)  # number of samples
p <- nrow(data)  # number of features

# number of samples for test set (30% test, 70% train).
nTest <- ceiling(n*0.3)
ind <- sample(n, nTest, FALSE)

# train set
data.train <- data[ ,-ind]
data.train <- as.matrix(data.train + 1)

treatmentGroup <- sample(factor(rep(c("Treatment", "Control"), c(n - nTest))), n - nTest, FALSE)
classtr <- data.frame(condition = class[-ind, ], treat = treatmentGroup)

# train set in S4 class
data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
                                       colData = classtr, formula(~ 1))

######## BURADA KALDIM ########
# test set
data.test <- data[ ,ind]
data.test <- as.matrix(data.test + 1)
treatmentGroupTest <- sample(factor(rep(c(1, 2), nTest)), nTest, FALSE)
classts <- data.frame(condition = class[ind, ], treat = treatmentGroupTest)

# test set (trivial sample)
# data.test <- data[ ,ind]
# data.test <- as.matrix(data.test + 1)
# data.test <- as.matrix(data.test[ ,1])
# colnames(data.test) <- "N1"
#
# classts <- data.frame(condition=class[ind, ])
# classts <- data.frame(condition = classts[1, ])

data.testS4 <- DESeqDataSetFromMatrix(countData = data.test,
                                      colData = classts, formula(~ 1))

##### 1. Continous Classifiers  #########

## Number of repeats (repeats) might change model accuracies ##
# Classification and Regression Tree (CART) Classification
model <- classify(data = data.trainS4, method = "rf", class.labels = "condition",
                  ref = "T", preProcessing = "deseq-vst",
                  control = trainControl(method = "repeatedcv", number = 5,
                                         repeats = 2, classProbs = TRUE), tuneLength = 15)
model

# predicted classes of test samples for CART method (class probabilities)
pred.continous = predict(model, data.testS4, type = "raw")
pred.continous

confusionMatrix(table(Predicted = pred.continous, Actual = classts$condition), positive = "T")


##### 2. voom-based Classifiers #########
vmodel <- classify(data = data.trainS4, normalize = "deseq", method = "voomNSC",
                   class.labels = "condition", ref = "T", control = voomControl(number = 5, repeats = 10, tuneLength = 30))
vmodel

## voom.train için predict fonksiyonu uyarlanacak.
pred.voom <- predict(vmodel, data.testS4)
pred.voom

confusionMatrix(table(Predicted = pred.voom, Actual = classts$condition), positive = "T")

##### 3. Discrete Classifier #######
##### 3.1 PLDA/PLDA2  ######
pri <- NULL
pmodel <- classify(data = data.trainS4, method = "PLDA", ref = "T",
                   class.labels = "condition",normalize = "deseq",
                   control = discreteControl(number = 5, repeats = 10,
                                             tuneLength = 30, parallel = TRUE))
pmodel
(pp <- plot(pmodel))

pred.PLDA <- predict(object = pmodel, data.testS4)
pred.PLDA

confusionMatrix(table(Predicted = pred.PLDA, Actual = classts$condition), positive = "T")

##### 3.2 NBLDA  #######
NBmodel <- classify(data = data.trainS4, method = "NBLDA", ref = "T",
                    class.labels = "condition", normalize = "deseq",
                    control = discreteControl(number = 5, repeats = 10, tuneLength = 10, parallel = TRUE))


pred.NBLDA <- predict(NBmodel, data.testS4)
pred.NBLDA

confusionMatrix(table(Predicted = pred.NBLDA, Actual = classts$condition), positive = "T")

