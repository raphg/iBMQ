
calculateThreshold <-
function(prob, threshold)
{ 
	
if(threshold > 1| threshold < 0){
  stop("The threshold need to be a numerical value between 0 and 1")
}
	
probvec=as.vector(as.matrix(prob))
probvecs=sort(probvec, decreasing = TRUE)
	
i <- 0
ok <- 0
while (ok < threshold) {
  i <- i+1
  ok <- mean(1-probvecs[1:i])
}
	
cutoff <- probvecs[i]
return(cutoff)

}
