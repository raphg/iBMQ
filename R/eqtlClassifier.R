eqtlClassifier <-
function(peak,posSNP,posGENE,max){


if(dim(peak)[2] != 3){
  stop("The peak data.frame need 3 colums.")
}
	
if(dim(posSNP)[2] != 3){
  stop("The SNP position data.frame need 3 colums.")
}
	
if(dim(posGENE)[2] != 4){
  stop("The gene position data.frame need 4 colums.")
}
	
#if(is.numerical(posGENE[,2]) == FALSE){
#	stop("The gene position data.frame need to be numerical")
#}

#if(is.numerical(posSNP[,2]) == FALSE){
#  stop("The snp position data.frame need to be numerical")
#}
	
	res1<-0
	res2<-0
	res3<-0
	res4<-0
	res5<-0
	res6<-""
	
for (i in 1:dim(peak)[1]) {
  Peak_gene <- as.character(peak[i,1])
  bool <- toupper(as.character(posGENE[,1])) %in% toupper(Peak_gene)
  if (any(bool)) {
			genechr <- posGENE[,2][bool][1]
            genestart <- posGENE[,3][bool][1]
            genestop <- posGENE[,4][bool][1]
        }
        else {
            cat("gene:", Peak_gene, "\t: position is missing\n")
			genechr <- "NA"
            genestart <- "NA"
            genestop <- "NA"
        }
		
		Peak_snp <- as.character(peak[i,2])
        bool <- toupper(as.character(posSNP[,1])) %in% toupper(Peak_snp)
        if (any(bool)) {
			chrsnp <- posSNP[,2][bool][1]
            pos <- posSNP[,3][bool][1]
           
        }
        else {
            cat("SNP",Peak_snp, "\t: position is missing\n")
            chrsnp <- "NA"
            pos <- "NA"
        }
		
		if(genechr!="NA" & chrsnp!="NA"){
		if(genechr==chrsnp)
		{
			if(abs(as.numeric(genestart)-as.numeric(pos))<=max){type="cis"}
			else{type="trans"}
		}
		else{type="trans"}
		
		}
				else{type="NA"}
		
		res1 <- c(res1,genechr)
		res2 <- c(res2,genestart)
		res3 <- c(res3,genestop)
		res4 <- c(res4,chrsnp)
		res5 <- c(res5,pos)
		res6 <- c(res6,type)

	}
	res <- data.frame(as.numeric(res1),as.numeric(res2),as.numeric(res3),as.numeric(res4),as.numeric(res5),res6)
	
	colnames(res) <- c("GeneChrm","GeneStart","GeneEnd","MarkerChrm","MarkerPosition","Type")
	
	return(cbind(peak,res[-1,]))
}