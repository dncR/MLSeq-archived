predictClassify <-
function(model, test.data){
    
    if (class(model) != "MLSeq") stop("Model should be a \"MLSeq\" object.")
    if (class(test.data)[1] != "DESeqDataSet") {stop("Data should be a \"DESeqDataSet\" object.")}
    
    VOOM = FALSE
    
    
    normalize = normalization(model)
    deseqTransform = deseqTransform(model)
    
    counts = counts(test.data)
    conditions = test.data$condition
    
    if (normalize != "none"){
        if (deseqTransform == "voom" & length(levels(conditions))<=1) {
            warning("Voom transformation can be applied only to factors with 2 or more levels. \"vst\" transformation is performed with DESeq's \"blind\" dispersion estimation method.")
            VOOM = TRUE
            deseqTransform = "vst"
            normalize = "deseq"
            vst.method = "blind"
        }
    }
    
    #### Data preperation steps: ##
    if (normalize == "none"){
        counts = t(counts)
    }
    
    if (normalize == "tmm" & !VOOM){
        counts = counts(test.data)
        y <- DGEList(counts=counts, genes=rownames(counts))
        y <- calcNormFactors(y)
        design <- model.matrix(~conditions)
        v <- voom(y,design,plot=FALSE)
        counts = data.frame(t(v$E))
    }
    
    
    if (normalize == "deseq"){
        test.data = estimateSizeFactors(test.data) #deseq size factor estimation
        
        if (deseqTransform == "vst") {
            test.data = estimateDispersions(test.data, fitType = "local")
            test.datavst = getVarianceStabilizedData(test.data) #vst transformation
            test.datavst <- as.matrix(test.datavst)
            conds=data.frame(conditions, row.names=colnames(test.datavst))
            test.datavst = ExpressionSet(test.datavst, AnnotatedDataFrame(conds))
            counts = data.frame(t(exprs(test.datavst))) #normalized gene expression data
        }
        
        if (deseqTransform == "voom" & !VOOM) {
            counts = counts(test.data, normalized=TRUE)
            y <- DGEList(counts=counts, genes=rownames(counts))
            design <- model.matrix(~conditions)
            v <- voom(y,design,plot=FALSE)
            counts = data.frame(t(v$E))
        }
    }
    
    test.pred = predict(trained(model), counts)
    test.pred
}

