setOldClass(c("confusionMatrix","train"))

# setClass for main class "MLSeq"
setClass("MLSeq",
			slots = c(method = "character",
					  deseqTransform="character",
					  normalization = "character", confusionMat="confusionMatrix",
					  trained="train", ref="character"),
			prototype = prototype(confusionMat=structure(list(), class="confusionMatrix"),
            trained=structure(list(), class="train")
		))


# generic function for class "MLSeq"
setGeneric("MLSeq", function(object) standardGeneric("MLSeq"))


#method for class MLSeq
setMethod("show",
signature = "MLSeq",
definition = function(object) {
    if (dim(confusionMat(object)$table)[1] == 2){
      cat("\n", sep = " ")
      cat("  An object of class ", class(object), "\n\n", sep = " ")
      cat("            Method  : ", method(object), "\n\n")
      cat("       Accuracy(%)  : ", round(confusionMat(object)$overall[1],4)*100, "\n")
      cat("    Sensitivity(%)  : ", round(confusionMat(object)$byClass[1],4)*100, "\n")
      cat("    Specificity(%)  : ", round(confusionMat(object)$byClass[2],4)*100, "\n\n")
      cat("  Reference Class   : ", ref(object), "\n\n")
    }
    else {    
      cat("\n", sep = " ")
      cat("  An object of class ", class(object), "\n\n", sep = " ")
      cat("            Method  : ", method(object), "\n\n")
      cat("       Accuracy(%)  : ", round(confusionMat(object)$overall[1],4)*100, "\n")
      cat("    Sensitivity(%)  : ", round(confusionMat(object)$byClass[1,1],4)*100, "\n")
      cat("    Specificity(%)  : ", round(confusionMat(object)$byClass[1,2],4)*100, "\n\n")
      cat("  Reference Class   : ", ref(object), "\n\n")
      invisible(NULL)
    }
})


# validity for class "MLSeq"
setValidity( "MLSeq", function( object ) {
    
    if (!(method(object)  %in% c("svm", "bagsvm", "randomforest", "cart")))
    return("Error: 'method' slot must be in one of the following methods: \"svm\", \"bagsvm\", \"randomforest\", \"cart\" ")
    
    if (!(normalization(object)  %in% c("deseq", "none", "tmm")))
    return("Error: 'normalization' slot must be in one of the following: \"deseq\", \"none\", \"tmm\" ")
    
    if (!(deseqTransform(object)  %in% c("vst", "voom", "NULL")))
    return("Error: 'deseqTransform' slot must be in one of the following: \"vst\", \"voom\" ")
    
    if (!is.character(ref(object)))
    return("Error: 'ref' slot must be a character ")
    
    if ((normalization(object) == "tmm" & deseqTransform(object) == "vst"))
    return("Warning: \"vst\" transformation can be applied only with \"deseq\" normalization. \"voom\" transformation is used. ")
       
    TRUE
} )





