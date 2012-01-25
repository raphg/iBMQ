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
