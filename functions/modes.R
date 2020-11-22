# Calculates the mode(s) of a vector
#x <-  c(8,2,7,1,2,9,8,2,10,9,8)
modes <- function(x, na.rm = TRUE){
  if(na.rm == TRUE){
	x <- x[!is.na(x)]
	x <- sort(x)
  	names(table(x))[table(x)==max(table(x))]
	}else{
	x <- sort(x)
        names(table(x))[table(x)==max(table(x))]
	}
}
