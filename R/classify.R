classify <-
function(data, method=c("svm","bagsvm","randomforest","cart"), normalize=c("deseq","none","tmm"),
deseqTransform=c("vst","voom"), cv=5, rpt=10, B=100, ref=NULL, ...){
    
    if(!is.null(ref)){
    if(!is.character(ref))stop("Reference class should be \"character\"")
    }
    if (is.null(ref)){ref=levels(data$condition)[1]}
    if (class(data)[1] != "DESeqDataSet") {stop("Data should be a \"DESeqDataSet Object\" of S4 class.")}
    if (is.null(method)) {stop("Classification method is not specified.")}
    
    
    method = match.arg(method)
    normalize = match.arg(normalize)
    deseqTransform = match.arg(deseqTransform)
    
    if ((normalize == "tmm" & deseqTransform == "vst")) {
        deseqTransform = "voom"
        warning("\"vst\" transformation can be applied with \"deseq\" normalization. \"voom\" transformation is used.")
    }
    
    if (normalize == "none") {
        deseqTransform = "NULL"
        warning("\"deseqTransform\" method is not used for normalize=\"none\" option")
    }
    
    conditions = as.factor(data$condition)
    conditions = relevel(conditions,which(levels(conditions)==ref))
    counts = counts(data)
    
    org.classes = conditions
    org.class.levels = levels(org.classes)
    
    ctrl <- trainControl(method = "repeatedcv", number=cv, repeats = rpt)
    
    #### Data preperation steps: ##
    if (normalize == "none"){
        counts = t(counts)
        conditions = conditions
        dataexp = data.frame(counts, conditions)
    }
    
    if (normalize == "tmm"){
        counts = counts(data)
        y <- DGEList(counts=counts, genes=rownames(counts))
        y <- calcNormFactors(y)
        design <- model.matrix(~conditions)
        v <- voom(y,design,plot=FALSE)
        
        dataexp = data.frame(t(v$E), conditions)
        counts = dataexp[,-length(dataexp)]
        conditions = as.factor(dataexp[, length(dataexp)])
    }
    
    
    if (normalize == "deseq"){
        data = estimateSizeFactors(data) #deseq size factor estimation
        
        if (deseqTransform == "vst") {
            data = estimateDispersions(data, fitType = "local")
            datavst = getVarianceStabilizedData(data) #vst transformation
            datavst <- as.matrix(datavst)
            cond=data.frame(conditions, row.names=colnames(datavst))
            datavst = ExpressionSet(datavst, AnnotatedDataFrame(cond))
            dataexp = data.frame(t(exprs(datavst)),conditions) #normalized gene expression data
            counts = dataexp[,-length(dataexp)]
            conditions = as.factor(dataexp[, length(dataexp)])
        }
        
        if (deseqTransform == "voom") {
            counts = counts(data, normalized=TRUE)
            y <- DGEList(counts=counts, genes=rownames(counts))
            design <- model.matrix(~conditions)
            v <- voom(y,design,plot=FALSE)
            dataexp = data.frame(t(v$E),conditions)
            counts = dataexp[,-length(dataexp)]
            conditions = as.factor(dataexp[, length(dataexp)])
        }
    }
    
    ######  Classification part (methods)
    if (method == "svm"){
        train <- train(counts, conditions, method = "svmRadial",
        trControl = ctrl, ...)
    }
    
    if (method == "bagsvm"){
        train <- train(counts, conditions, method="bag", B = B,
        bagControl = bagControl(fit = svmBag$fit,
        predict = svmBag$pred,
        aggregate = svmBag$aggregate),
        trControl = ctrl, ...
        )
    }
    
    if (method == "randomforest"){
        train <- train(counts, conditions, method = "rf",
        trControl = ctrl, ...)
    }
    
    if (method == "cart"){
        train <- train(counts, conditions, method = "rpart",
        trControl = ctrl, ...)
    }
    
    train.pred = predict(train)
    tbl.trn = table(Predicted=train.pred, Actual=org.classes)
    confM = confusionMatrix(tbl.trn,reference=ref)
    
    result = new("MLSeq", confusionMat=confM, trained=train, method=method, normalization=normalize,
    deseqTransform=deseqTransform, ref=ref)
    
    result
}
