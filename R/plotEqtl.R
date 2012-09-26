plotEqtl <-
function(peak){
	

peak=peak[complete.cases(peak),]
	
if(dim(peak)[2]!=9){
  stop("The peak data.frame need to have 9 columns.")
}

t1 <- unique(peak[,4])
t2 <- unique(peak[,7])
numiChr <- max(c(length(t1), length(t2)))
maxChr <- max(c(peak[,6], peak[,8]))
numChr <- max(c(peak[,4], peak[,7]))
	
if(numiChr == 1){
  chr2 <- c(0,maxChr)
}
	
else{
chr <- c(0,rep(maxChr,numChr))
chr2 <- chr
  for (i in 2:length(chr))
	chr2[i] <- sum(as.numeric(chr[1:i]))
}
	
Xvect <- numeric(dim(peak)[1])
Yvect <- numeric(dim(peak)[1])

for(i in 1:dim(peak)[1]){
  Yvect[i] <- chr2[(as.numeric(peak[i,4]))]+as.numeric(peak[i,5])
  Xvect[i] <- chr2[(as.numeric(peak[i,7]))]+as.numeric(peak[i,8])
}
	
	
if(numiChr == 1){
  plot(Xvect, Yvect,ylim=c(0,chr2[length(chr2)]), xlim=c(0,chr2[length(chr2)]), cex=0.5, pch=1, xlab="SNP", ylab="GENE", xaxt="n" ,yaxt="n", main="Plot of eQTLs in relation to the positions of SNPs and genes")
}
else{if(numiChr > 1){
  plot(Xvect, Yvect,ylim=c(0,chr2[length(chr2)]), xlim=c(0,chr2[length(chr2)]), cex=0.5, pch=1, xlab="SNP", ylab="GENE", xaxt="n" ,yaxt="n", main="Plot of eQTLs in relation to the positions of SNPs and genes")
  abline(h=chr2,col="grey")
  abline(v=chr2,col="grey")
}}



}

