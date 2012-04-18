eqtl.mcmc <-
function(snp,expr,n.iter,burn.in,n.sweep,nproc, constC = TRUE, write.output = TRUE)
{
n.pheno <- length(expr)
n.snp <- length(snp)

if(!is.null(expr))
{
n.indiv <- length(expr[[1]])
names.indiv <- rownames(expr)
}

if(!is.null(snp))
{
n.indiv1 <- length(snp[[1]])-2
names.indiv1 <- rownames(snp)[-(1:2)]
chr <- snp$chr
pos <- snp$pos
snp <- unlist(snp[-(1:2),1:n.snp])
}

if(!is.null(expr)&!is.null(snp))
{
if(n.indiv!=n.indiv1)
{
cat("indivs number in expr is",n.indiv,"\n")
cat("indivs number in snp is",n.indiv1,"\n")
stop("there is a problem with indiv numbers")
}
if(any(names.indiv!=names.indiv1))
{
cat("indivs name in expr are",names.indiv,"\n")
cat("indivs name in snp are",names.indiv1,"\n")
stop("there is a problem with indivs names")
}
}

pheno <- unlist(expr)
#call "c_qtl_mcmc" of C code
start.time <- Sys.time()
outProb = double(length = n.snp*n.pheno)

c.function="c_qtl_main_parallel_sparse"
#res <- .C("c_qtl_mcmc_wjg_new",as.double(pheno),as.integer(n.indiv),as.integer(n.pheno),as.double(snp),as.integer(n.snp),as.integer(n.iter),as.integer(burn.in),as.integer(n.sweep))
eps = 10*(.Machine$double.eps)
nmax = 500
res <- .C(c.function,as.double(pheno),as.integer(n.indiv),
		as.integer(n.pheno),as.double(snp),as.integer(n.snp),
		as.integer(n.iter),as.integer(burn.in),as.integer(n.sweep),
		as.double(outProb), as.integer(nproc), as.integer(nmax), as.double(eps),
		as.integer(write.output), as.integer(!constC))
end.time <- Sys.time()
cat("running MCMC takes ")
cat(as.character(round(difftime(end.time,start.time,units="min"),digits=2)))
cat(" minutes.\n")
dim(res[[9]]) = c(n.snp, n.pheno)
return(res[[9]])
}

