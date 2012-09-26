eqtlMcmc <-
function(snp,expr,n.iter,burn.in,n.sweep,nproc, constC = TRUE, write.output = TRUE, RIS=TRUE)
{
	
#snp=snp[complete.cases(snp),]
#expr=expr[complete.cases(expr),]
cat("Incomplete data will be remove","\n")
	
n.pheno <- length(as.data.frame(expr))
n.snp <- length(as.data.frame(snp))
	
	

if(!is.null(expr)){
  n.indiv <- length(as.data.frame(expr)[[1]])
  names.indiv <- rownames(expr)
}

if(!is.null(snp)){
  n.indiv1 <- length(as.data.frame(snp)[[1]])
  names.indiv1 <- rownames(snp)
  names.snp<- colnames(snp)
  snp <- unlist(snp)
}

if(!is.null(expr) &! is.null(snp)){
  if(n.indiv != n.indiv1){
    cat("indivs number in expr is",n.indiv,"\n")
    cat("indivs number in snp is",n.indiv1,"\n")
    stop("there is a problem with indiv numbers")
  }
  if(any(names.indiv != names.indiv1)){
    cat("indivs name in expr are",names.indiv,"\n")
    cat("indivs name in snp are",names.indiv1,"\n")
    stop("there is a problem with indivs names")
  }
}



if(RIS == TRUE){
  if(!all((snp == 1)| (snp == 0))){
    stop("The SNPs value need to be 0 and 1")
  }
  snp[snp==0] <- -0.5
  snp[snp==1] <- 0.5
}
	
else{
  if(!all((snp == 1)| (snp == 2) | (snp == 3))){
    stop("The SNPs value need to be 1,2,3")}
#snp[snp==0] <- -0.5
#snp[snp==1] <- 0
#snp[snp==2] <- 0.5
}
	
pheno <- unlist(expr)
#call "c_qtl_mcmc" of C code
start.time <- Sys.time()
outProb <- double(length = n.snp*n.pheno)


c.function="c_qtl_main_parallel_sparse"
#res <- .C("c_qtl_mcmc_wjg_new",as.double(pheno),as.integer(n.indiv),as.integer(n.pheno),as.double(snp),as.integer(n.snp),as.integer(n.iter),as.integer(burn.in),as.integer(n.sweep))
eps <- 10*(.Machine$double.eps)
nmax <- 500
res <- .C(c.function,as.double(pheno),as.integer(n.indiv),
		as.integer(n.pheno),as.double(snp),as.integer(n.snp),
		as.integer(n.iter),as.integer(burn.in),as.integer(n.sweep),
		as.double(outProb), as.integer(nproc), as.integer(nmax), as.double(eps),
		as.integer(write.output), as.integer(!constC))
end.time <- Sys.time()
cat("running MCMC takes ")
cat(as.character(round(difftime(end.time,start.time,units="min"),digits=2)))
cat(" minutes.\n")
dim(res[[9]]) <- c(n.snp, n.pheno)
out<-res[[9]]
colnames(out)<-colnames(expr)
rownames(out)<-names.snp
return(out)
}

