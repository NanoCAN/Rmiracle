rppa.normalize.depos <- function(spots){
  library(dplyr)
  spots$Deposition <- as.factor(spots$Deposition)
  print(head(spots))
  
  normalize.lm <- function(x){
    correctionFactor <- coef(lm(concentrations~Deposition, x))
    #replace intercept with 1 for 1st deposition
    #intercept <- correctionFactor[1]
    #correctionFactor[1] <- 1
    
    result <- x %>% group_by(Deposition, add=T) %>% mutate(concentrations=concentrations/correctionFactor[which(levels(x$Deposition)==unique(Deposition))],
                                                           upper=upper/correctionFactor[which(levels(x$Deposition)==unique(Deposition))],
                                                           lower=lower/correctionFactor[which(levels(x$Deposition)==unique(Deposition))])
    return(result)
  }
  
  spots %>% group_by(A,B) %>% do(normalize.lm(.))
}