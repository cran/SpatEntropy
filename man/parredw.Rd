% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneill_leibovici_entropy.R
\name{parredw}
\alias{parredw}
\title{Parresol and Edwards' entropy.}
\usage{
parredw(data)
}
\arguments{
\item{data}{A data matrix or vector, can be numeric, factor, character, ...}
}
\value{
a list of two elements:
\itemize{
\item \code{probabilities} - a table with the estimated probabilities (relative frequencies) for all couple categories;
\item \code{parredw} - Parresol and Edwards' entropy.
}
}
\description{
Compute Parresol and Edwards' entropy, following Parresol and Edwards (2014),
starting from data. References can be found at \code{SpatEntropy}.
}
\details{
This index is based on the transformed variable \eqn{Z} identifying couples of realizations
of the variable of interest. A distance of interest is fixed: Parresol and Edwards' entropy is
thought for areas sharing a border, as O'Neill's entropy. All contiguous couples
of realizations of the variable of interest are counted
and their relative frequencies are used to compute the index, which is the opposite of O'Neill's entropy.
The function is able to work with grids containing missing data, specified as NA values.
All NAs are ignored in the computation and only couples of non-NA observations are considered.
}
\examples{
#numeric data, square grid
data=matrix(sample(1:5, 100, replace=TRUE), nrow=10)
parredw(data)
#plot data
plot(as.im(data, W=square(nrow(data))),
     col=gray(seq(1,0,l=length(unique(c(data))))),
     main="", ribbon=TRUE)

#character data, rectangular grid
data=matrix(sample(c("a","b","c"), 300, replace=TRUE), nrow=30)
parredw(data)
#plot data
plot(as.im(data, W=owin(xrange=c(0,ncol(data)), yrange=c(0,nrow(data)))),
     col=terrain.colors(length(unique(c(data)))),
     main="", ribbon=TRUE)

}
