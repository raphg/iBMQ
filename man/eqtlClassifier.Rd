\name{eqtlClassifier}
\alias{eqtlClassifier}
\title{
Classifying the eQTLs
}
\description{
It is customary to distinguish two kinds of eQTLs: 1) cis-eQTLs (where the eQTL is on the same locus as the expressed gene); and 2) trans-eQTLs (where the eQTL is on a locus other than that of the expressed gene). 
The eqtlClassifier allows us to classify the eQTLs as eithercis-eQTL or trans-eQTL according to their position in the genome.
}
\usage{
eqtlClassifier(peak, posSNP, posGENE, max)
}

\arguments{
  \item{peak}{
A data.frame of significant eQTLs (the output of the findEqtl function)
}
  \item{posSNP}{
A  data frame specifying the genomic locations of genomic markers (i.e. SNPs).
}
  \item{posGENE}{
A data frame specifying the genomic locations of the genes (or probes).
}
  \item{max}{
A cutoff value (in base pair) corresponding to the threshold where a eQTLis considered to be cis-eQTL. A numerical value.}
}

\value{
The output of the eqtlClassifier is a data frame where the first column contains the names of each gene, the second column contains the names of markers and the third column contains the PPA value for each significant eQTL. The fourth column contains the number of the chromosome to which the gene belongs, the fifth column contains the start position of the geneand the sixth column contains the end position of the gene. The seventh column contains the number of the chromosome to which the marker belongs, the eighth column contains position of the marker and the ninth column contains a descriptor of the type of eQTL (either cis or trans). Please note that in order to ascertain that an eQTL is either cis or trans, the positions of the markers and the gene need to the given to the function. If one of the values is missing the type of eQTL will be "NA".
}


