\name{eqtlFinder}
\alias{eqtlFinder}
\title{
eqtlFinder
}
\description{
We can calculate how many eQTLs have PPA above the cutoff with the eqtlFinder function.
}
\usage{
find.eqtl(prob, threshold)
}

\arguments{
  \item{prob}{matrix or data frame that contains Posterior   Probability of Association  (output of eqtl.mcmc function)
}
  \item{threshold}{Threshold to use for significant QTL. This value can be the output of the calculateThreshold function.It must be a numerical value between 0 and 1.
}
}

\value{
The output of the eqtlFinderis a data frame where the first column contains the names of each gene, the second columncontains the names of corresponding markers and the third column contains the PPA valuefor each significant eQTL.

}

