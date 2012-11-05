\name{eqtlFinder}
\alias{eqtlFinder}
\title{
eqtlFinder
}
\description{
We can calculate how many eQTLs have PPA above the cutoff with the eqtlFinder function.
}
\usage{
eqtlFinder(prob, threshold)
}

\arguments{
  \item{prob}{matrix or data frame that contains the Posterior   Probability of Association  values (output of eqtlMcmc function)
}
  \item{threshold}{Threshold to be used to determine which QTLs are significant. This value can be the output of the calculateThreshold function.It must be a numerical value between 0 and 1.
}
}

\value{
The output of the eqtlFinder is a data frame where the first column contains the names of each gene, the second column contains the names of corresponding markers and the third column contains the PPA value for each significant eQTL.
}

\examples{
data(PPA.liver)
cutoff.liver <- calculateThreshold(PPA.liver, 0.2)
eqtl.liver <- eqtlFinder(PPA.liver, cutoff.liver)

}

