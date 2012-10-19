\name{eqtlMcmc}
\alias{eqtlMcmc}
\title{
Bayesian Multiple eQTL mapping using MCMC
}

\description{

Compute MCMC algorithm to produce Posterior Probability of Association for 
  eQTL mapping.   }
\usage{
eqtlMcmc(snp, expr, n.iter, burn.in, n.sweep, mc.cores,  constC = TRUE,
 write.output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ SnpSet class object}
  \item{expr}{ ExpressionSet class object}
  \item{n.iter}{Number of iterations to be saved}
  \item{burn.in}{The initial burn-in period}
  \item{n.sweep}{The sweeping number which must be a positive number}
  \item{mc.cores}{The number of cores you would like to use for 
  	parallel processing. Can be set be set via `options(cores=4)', if not set, 
	the code will automatically detect the number of cores.} 
  \item{constC}{Should variables C be fixed? If TRUE, fixed to 
  	number of SNPS. If FALSE, C parameters are given independent
  	inverse gamma priors.}
  \item{write.output}{Write chain iterations to file? If TRUE, 
  	output for variables will be written to files created in 
  	the working directory.}
  \item{RIS}{If TRUE, the genotype need to be 0 and 1. If FALSE the genotype need to be 1,2 and 3.}
}

\value{
A matrix with Posterior Probability of Association. Rows correspond to snps
from original snp data object, columns correspond to genes from expr data object.
}
