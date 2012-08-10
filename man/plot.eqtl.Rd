\name{plot.eqtl}
\alias{plot.eqtl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the genome-wide eQTL mapping
}
\description{
The function to visualize the genome-wide eQTL mapping.
}
\usage{
plot.eqtl(prob, threshold, posSNP, posGENE, chr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob}{matrix or data frame that contains Posterior Probability of          Association (output of eqtl.mcmc function).
}
  \item{threshold}{threshold  of significant for a association
}
  \item{posSNP}{A data frame specifying the genomic locations of the SNP with columns: SNP name, the chromosome, the start location in base pair.}
  \item{posGENE}{A data frame specifying the genomic locations of the gene/probe with columns: gene name, the chromosome, the start location in base pair and the stop location in base pair.}
  \item{chr}{A vector of integer specifying the size of the chromosomes in base pair in order of the chromosomes.}
}
\examples{#To calculate the threshold leading to 10\% FDR#PPA= eqtl.mcmc(snp, gene, n.iter=1000,burn.in=1000,n.sweep=10,nproc=2,write.output=FALSE)ThresholdPPA =calculateThreshold(PPA,0.1)
data(snppos)
data(genepos)
mouseChr=c(1.97e+08, 1.82e+08, 1.60e+08, 1.56e+08, 1.53e+08, 1.50e+08,
1.53e+08, 1.32e+08,1.24e+08, 1.30e+08,1.22e+08,1.21e+08,1.20e+08, 1.25e+08,1.03e+08, 9.80e+07, 9.50e+07, 9.10e+07, 6.10e+07)
plot.eqtl(PPA,cutoff,snppos, genepos,mouseChr)
}