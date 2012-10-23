\name{hotspotFinder}
\alias{hotspotFinder}

\title{
hotspotFinder
}
\description{
One main advantage of our method is its increased sensitivity for finding trans-eQTL hotspots (corresponding to situations where a single SNP is linked to the expression of several genes across the genome).
}
\usage{
hotspotFinder(peak, numgene)
}
\arguments{
  \item{peak}{
A data frame (3 columns) corresponding to the output of the eqtlFinder function or the data frame (9 columns) corresponding to the output of the eqtlClassifier function. 

}
  \item{numgene}{
The minimum of gene to detect.
}
}

\value{
The output of this function is a list, where each element is a marker. For each markerthere is a data frame with all the eQTLs linked to this marker.
}
