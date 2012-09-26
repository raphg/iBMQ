eqtlFinder <-
function(prob,threshold){


if(threshold>1|threshold<0){
  stop("The threshold need to be a numerical value between 0 and 1")
}
	
gene <- 0
snp <- 0
value <- 0
	
for (i in 1:dim(prob)[1]){
		
lis <- which(prob[i,] >= threshold)
	
if (length(lis)>=1){
  for (j in 1:length(lis)){
	gene <- c(gene,colnames(prob)[lis[j]])
	snp <- c(snp,rownames(prob)[i])
	value <- c(value,prob[i,lis[j]])
	}
  }
}
	
res <- data.frame(as.character(gene[-1]),as.character(snp[-1]),as.numeric(value[-1]))
colnames(res) <- c('Gene','SNP','PPA')
return(res)
	
}

