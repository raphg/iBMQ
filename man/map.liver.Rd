\name{map.liver}
\alias{map.liver}
\docType{data}
\title{
SNP position data frame 
}
\description{
A data frame specifying the genomiclocations of each SNP with following columns: SNP name, chromosome number, SNP location(in base pair).
}
\format{
  A data frame with 290 observations on the following 3 variables.
\describe{
    \item{\code{V1}}{a character vector}
    \item{\code{V2}}{a numeric vector}
    \item{\code{V3}}{a numeric vector}
  }
  }
\source{
 This F2 cross data set containing genotypic and phenotypic information for60 mice was obtained from the lab of Alan Attie Lab at the University of Wisconsin-Madison.These data are also available at GEO (accession number GSE3330). Only the 5000 most variableexpression traits out of 45,265 transcripts from the liver were used for the current example.
}

\examples{
data(map.liver)
}
\keyword{datasets}
