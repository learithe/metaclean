% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metagenomic_functions.R
\name{RemoveTaxa}
\alias{RemoveTaxa}
\title{RemoveTaxa Function}
\usage{
RemoveTaxa(physeq, badTaxa, prune = FALSE)
}
\arguments{
\item{physeq}{phyloseq object to prune}

\item{badTaxa}{list of taxa (OTU IDs) to remove}

\item{prune}{remove samples that have no remaining taxa (default FALSE)}
}
\description{
Function to remove unwanted taxa from a physeq object
}
\examples{
ps_noblanks <- RemoveTaxa(ps, blank_taxa, prune=TRUE)
}
\keyword{badtaxa,}
\keyword{removetaxa}
