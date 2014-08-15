rppa.hshift <- function(spots){
  if(!is.null(attr(spots, "hshifted")))
  {
    cat("This slide has already been hshifted! Do you really want to continue? (yes/no)")
    answer <- readline()
    
    if(answer != "yes") return()
  }
  
  library(dplyr)
  spots <- spots %.% group_by(Block, Row) %.% arrange(Column) %.% mutate_each(funs(rppa.shift.vector(., min(hshift))), Signal, FG, BG, Flag, Diameter)
  return(spots)
}

rppa.vshift <- function(spots){
  if(!is.null(attr(spots, "vshifted"))){
    
    cat("This slide has already been vshifted! Do you really want to continue?")
    answer <- readline()
    
    if(answer != "yes") return()
  }
  
  library(dplyr)
  spots <- spots %.% group_by(Block, Column) %.% arrange(Row) %.% mutate_each(funs(rppa.shift.vector(., min(vshift))), Signal, FG, BG, Flag, Diameter)
  return(spots)
}


rppa.shift.vector <-  function (v, by) 
{
  if (by == 0) return(v)
  
  else if (by >= 0) return(c(NA + 1:by, v[1:(length(v) - by)]))
  
  else return(c(v[(-by + 1):(length(v) - by)]))
}
