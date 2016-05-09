SampleGroup <- setClass("SampleGroup", slots = c(description = "character",
                                                 values = "vector"))
                                                
RPPArray <- setClass("RPPArray",
                   slots = c(title = "character",
                             PMT = "integer",
                             blocksPerRow = "integer",
                             antibody = "character",
                             slideId = "integer",
                             barcode = "character",
                             positionData = "data.frame",
                             FG = "vector",
                             BG = "vector",
                             Signal = "vector",
                             spotData = "data.frame",
                             phenoData = "data.frame",
                             sampleGroupA = "SampleGroup",
                             sampleGroupB = "SampleGroup",
                             sampleGroupC = "SampleGroup",
                             plateData = "data.frame",
                             vshift = "vector",
                             hshift = "vector",
                             log = "logical",
                             shifted = "logical",
                             normalized = "logical",
                             normalized.method = "character"
                             ))

setMethod("plot", "RPPArray", function(x, fill="Signal"){
  
  plot.data <- cbind(x@positionData, FG=x@FG, BG=x@BG, Signal=x@Signal, x@spotData)
  
  if(!(fill %in% colnames(plot.data))) stop("invalid fill property selected")
  
  rppa.plot.heatmap(plot.data, fill=fill)
})

setMethod("log2", "RPPArray", function(x){
  
  if(x@log) stop("RPPArray object has already been log transformed.")
  if(x@normalized) stop("RPPArray object has already been normalized. Apply log transformation prior to normalization.")
  else{
   x@Signal <- log2(x@Signal)
   x@FG <- log2(x@FG)
   x@BG <- log2(x@BG)
  }
})
