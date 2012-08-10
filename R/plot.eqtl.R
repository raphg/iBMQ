plot.eqtl <-
function(prob,threshold,posSNP,posGENE,chr){
prob=t(prob)

	
if(length(chr)<=1)
{
		nsnp=max(max(posSNP[,2]),max(posGENE[,2]))
		nlength=max(max(posSNP[,3]),max(posGENE[,3]))
		
		chr=rep(nlength,nsnp)
		chr2=chr
		for (i in 2:length(chr))
		chr2[i]=sum(as.numeric(chr[1:i]))
}
else {if(!is.null(chr))
{
chr2=chr
for (i in 2:length(chr))
chr2[i]=sum(as.numeric(chr[1:i]))
}}

	
snpVector=numeric(dim(posSNP)[1])

for (i in 1:dim(posSNP)[1]){
if (as.numeric(posSNP[i,2])==1){
snpVector[i]=as.numeric(posSNP[i,2])}
else{
snpVector[i]=chr2[(as.numeric(posSNP[i,2])-1)]+as.numeric(posSNP[i,3])
}
}

names(snpVector)=posSNP[,1]

geneVector=numeric(dim(posGENE)[1])

for (i in 1:dim(posGENE)[1]){
if (as.numeric(posGENE[i,2])==1){
geneVector[i]=as.numeric(posGENE[i,3])}
else{
geneVector[i]=chr2[(as.numeric(posGENE[i,2])-1)]+as.numeric(posGENE[i,3])
}
}

names(geneVector)=as.character(posGENE[,1])

colnames(prob)=names(snpVector)
rownames(prob)=names(geneVector)

Xvect=0
Yvect=0

for(i in 1:dim(prob)[1]){
lis=which(prob[i,]>=threshold)
if (length(lis)>=1){
for (j in 1:length(lis)){
Xvect=c(Xvect,geneVector[rownames(prob)[i]])
Yvect=c(Yvect,snpVector[colnames(prob)[lis[j]]])}}

}

plot(Yvect, Xvect,ylim=c(0,chr2[length(chr2)]), xlim=c(0,chr2[length(chr2)]), cex=0.5, pch=1, xlab="SNP", ylab="GENE", xaxt="n" ,yaxt="n", main="Plot of eQTLs in relation to the positions of SNPs and genes")
abline(h=chr2,col="grey")
abline(v=chr2,col="grey")



}

