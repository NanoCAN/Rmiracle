rppa.shift.vector <-  function (v, by) 
{
  if (by == 0) return(v)
  
  else if (by >= 0) return(c(NA + 1:by, v[1:(length(v) - by)]))
  
  else return(c(v[(-by + 1):(length(v) - by)]))
}

rppa.hshift <- function(spots, blocks=NA, cols=NA, by=NA)
{
  if(!is.null(attr(spots, "hshifted")))
  {
    cat("This slide has already been hshifted! Do you really want to continue? (yes/no)")
    answer <- readline()
  
    if(answer != "yes") return()
  }
  if(is.na(by) && sum(spots$hshift) == 0){
    message("nothing to do")
    return(spots)
  }
  
  spots <- rppa.shift(spots, blocks, cols, by, "h") 
  
  attr(spots, "hshifted") <- TRUE
  
  return(spots)
}

rppa.vshift <- function(spots, blocks=NA, rows=NA, by=NA)
{
  if(!is.null(attr(spots, "vshifted"))){
    
    cat("This slide has already been vshifted! Do you really want to continue?")
    answer <- readline()
    
    if(answer != "yes") return()
  }
  
  if(is.na(by) && sum(spots$vshift) == 0){
    message("nothing to do")
    return(spots)
  }
  spots <- rppa.shift(spots, blocks, rows, by, "v") 

  attr(spots, "vshifted") <- TRUE
  return(spots)
}

rppa.shift <- function(spots, blocks, selected.subset, by, direction="v"){
  if(is.na(blocks[1])) range <- min(spots$Block):max(spots$Block)
  else range <- blocks
  
  for(b in range)
  {
    blockB <- subset(spots, Block==b);
    
    if(direction == "v")
    {
      splitVar <- "Column"
      if(!is.na(by)) by <- by
      else if(sum(blockB$vshift) == 0) next  
      else by <- blockB$vshift[1]
    } else if(direction == "h"){
      splitVar <- "Row"
      if(!is.na(by)) by <- by
      else if(sum(blockB$hshift) == 0) next
      else by <- blockB$hshift[1]
    }
    else stop("direction can only be h for horizontal or v for vertical")
    
    if(by > 0){
      spots[spots$Block==b,] <- unsplit(
        lapply(split(blockB, blockB[,splitVar]), function(x, selected.subset, by)
        {
          if(direction=="v") x <- x[with(x, order(Row)),]  
          else if(direction=="h") x <- x[with(x, order(Column)),]  
          
          for(field in c("Signal", "FG", "BG", "Flag", "Diameter")){         
            if(!is.na(selected.subset)) x[selected.subset,field] <-  rppa.shift.vector(x[selected.subset, field], by)       
            else x[,field] <-  rppa.shift.vector(x[, field], by)       
          }
              
          
          return(x);
        }, selected.subset=selected.subset, by=by)  
        , blockB[,splitVar])
    }
  }
  return(spots)
}