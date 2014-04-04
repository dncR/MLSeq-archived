
## ----echo=FALSE, message=FALSE-------------------------------------------
require(knitr)
opts_chunk$set(cache = TRUE, dev = "pdf")


## ----warning=FALSE, message=FALSE----------------------------------------
library(MLSeq)


## ----chunk1--------------------------------------------------------------
filepath = system.file("extdata/cervical.txt", package = "MLSeq")
filepath


## ----chunk2, message=FALSE-----------------------------------------------
cervical = read.table(filepath, header=TRUE)


## ----chunk3--------------------------------------------------------------
head(cervical[,1:5])


## ----chunk4--------------------------------------------------------------
class(cervical)
dim(cervical)


## ----chunk5, echo=FALSE--------------------------------------------------
options(width = 63)


## ----chunk6--------------------------------------------------------------
class = data.frame(condition = factor(rep(c("N","T"), c(29,29))))
as.factor(class[,1])


## ----chunk7--------------------------------------------------------------
data = cervical[c(1:150),]


## ----chunk8--------------------------------------------------------------
nTest = ceiling(ncol(data)*0.2)  
set.seed(12345) 
ind = sample(ncol(data), nTest, FALSE)


## ----chunk9--------------------------------------------------------------
data.train = data[,-ind]    
data.train = as.matrix(data.train + 1)
classtr = data.frame(condition = class[-ind,])


## ----chunk10-------------------------------------------------------------
data.test = data[,ind]
data.test = as.matrix(data.test + 1)
classts = data.frame(condition = class[ind,])


## ----chunk11-------------------------------------------------------------
dim(data.train)
dim(data.test)


## ----chunk12, message=FALSE----------------------------------------------
data.trainS4 = DESeqDataSetFromMatrix(countData = data.train,
colData = classtr, formula(~condition))
data.trainS4 = DESeq(data.trainS4, fitType="local")
data.trainS4


## ----chunk13, message=FALSE----------------------------------------------
data.testS4 = DESeqDataSetFromMatrix(countData = data.test,
colData = classts, formula(~condition))
data.testS4 = DESeq(data.testS4, fitType = "local") 
data.testS4


## ----chunk14, message=FALSE, tidy.opts=list(width.cutoff=65)-------------
svm = classify(data = data.trainS4, method = "svm", normalize = "deseq",
               deseqTransform = "vst", cv = 5, rpt = 3, ref = "T")
svm


## ----chunk15-------------------------------------------------------------
getSlots("MLSeq")
trained(svm)


## ----chunk16, message=FALSE, tidy.opts=list(width.cutoff=65)-------------
rf = classify(data = data.trainS4, method = "randomforest", normalize = "deseq", 
              deseqTransform = "vst", cv = 5, rpt = 3, ref = "T")
rf


## ----chunk17, message=FALSE----------------------------------------------
pred.svm = predictClassify(svm, data.testS4)
pred.svm


## ----chunk18, message=FALSE----------------------------------------------
pred.rf = predictClassify(rf, data.testS4) 
pred.rf


## ----chunk19-------------------------------------------------------------
table(pred.svm, relevel(data.testS4$condition, 2))
table(pred.rf, relevel(data.testS4$condition, 2))


## ----chunk20-------------------------------------------------------------
sessionInfo()


