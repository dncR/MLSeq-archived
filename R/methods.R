# method for slot "method"
setMethod("method", "MLSeq", function(object) object@method)

# method for slot "deseqTransform"
setMethod("deseqTransform", "MLSeq", function(object) object@deseqTransform)

# method for slot "normalization"
setMethod("normalization", "MLSeq", function(object) object@normalization)

# method for slot "confusionMat"
setMethod("confusionMat", signature = "MLSeq", function(object) object@confusionMat)

# method for slot "trained"
setMethod("trained", signature = "MLSeq", function(object) object@trained)

# method for slot "ref"
setMethod("ref", "MLSeq", function(object) object@ref)

