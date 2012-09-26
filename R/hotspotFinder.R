hotspotFinder <-
function(peak,numgene){
	
	
if((dim(peak)[2] != 9) & (dim(peak)[2] !=3 )){
  stop("The peak data.frame need to be 3 or 9 columns.")
}
	   
out2 <- split(peak, peak[,2])
lapply(out2,dim)
out2 <- out2[which(lapply(out2,nrow) >= numgene)]
return(out2)
}

