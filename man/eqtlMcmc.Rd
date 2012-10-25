\name{eqtlMcmc}
\alias{eqtlMcmc}
\title{
Bayesian Multiple eQTL mapping using MCMC
}

\description{

Compute the MCMC algorithm to produce Posterior Probability of Association values for 
  eQTL mapping.   }
\usage{
eqtlMcmc(snp, expr, n.iter, burn.in, n.sweep, mc.cores,
 write.output = TRUE, RIS = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ SnpSet class object}
  \item{expr}{ ExpressionSet class object}
  \item{n.iter}{Number of samples to be saved from the Markov Chain}
  \item{burn.in}{Number of burn-in iterations for the Markov Chain}
  \item{n.sweep}{Number of iterations between samples of the Markov Chain (AKA thinning interval)}
  \item{mc.cores}{The number of cores you would like to use for 
  	parallel processing. Can be set be set via `options(cores=4)', if not set, 
	the code will automatically detect the number of cores.} 
  \item{write.output}{Write chain iterations to file. If TRUE, 
  	output for variables will be written to files created in 
  	the working directory.}
  \item{RIS}{If TRUE, the genotype needs to be either 0 and 1. If FALSE the genotype need to be either 1,2 and 3.}
}

\value{
A matrix with Posterior Probability of Association values. Rows correspond to snps
from original snp data objects, columns correspond to genes from expr data objects.
}
