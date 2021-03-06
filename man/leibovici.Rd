% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneill_leibovici_entropy.R
\name{leibovici}
\alias{leibovici}
\title{Leibovici's entropy.}
\usage{
leibovici(data, cell.size = 1, ccdist = cell.size, verbose = F)
}
\arguments{
\item{data}{If data are lattice, a data matrix, which can be numeric, factor, character, ...
If the dataset is a point pattern, \code{data} is a \code{ppp} object.}

\item{cell.size}{A single number. If data are lattice, the length of the side of each pixel. Default to 1. Ignored if data are points.}

\item{ccdist}{A single number. The chosen distance for selecting couples of pixels/points within the observation area. Default to \code{cell.size}.}

\item{verbose}{Logical. If \code{TRUE} an output is printed in order to follow the progress of the work (recommended for large dataset).
Default set to \code{FALSE}.}
}
\value{
a list of three elements:
\itemize{
\item  \code{probabilities} - a table with the estimated probabilities (relative frequencies) for all couple categories;
\item \code{leib} - Leibovici's entropy;
\item \code{rel.leib} - Leibovici's relative entropy.
}
}
\description{
This function computes Leibovici's entropy according to a chosen distance \eqn{d}
(with O'Neill's entropy as a special case)
following Leibovici (2009), see also Altieri et al (2017). References can be found at \code{SpatEntropy}.
}
\details{
This index is based on the transformed variable \eqn{Z} identifying couples of realizations
of the variable of interest. A distance of interest is fixed, which in the case of O'Neill's
entropy is the contiguity, i.e. sharing a border for lattice data. Then, all couples
of realizations of the variable of interest lying at a distance smaller or equal to the distance of interest
are counted, and their relative frequencies are used to compute the index with the traditional Shannon's
formula.
#'\deqn{H_L(Z)=\sum p(z_r|d) \log(1/p(z_r|d))}
where \eqn{z_r=(x_i, x_{i'})} is a generic couple of realizations of the study variable \eqn{X}.
The conditioning on \eqn{d} means that only couples within a predefined distance are considered.
The maximum value for Leibovici's entropy is
\eqn{\log(I^2)} where \eqn{I} is the number of categories of the study variable \eqn{X}.
The relative version of Leibovici's entropy is obtained by dividing the entropy value by its maximum, and is useful for comparison across datasets with
a different number of categories.
The function is able to work with grids containing missing data, specified as NA values.
All NAs are ignored in the computation and only couples of non-NA observations are considered.
}
\examples{
#random grid data - high entropy
data=matrix(sample(c("a","b","c"), 400, replace=TRUE), nrow=20)
leibovici(data, cell.size=1, ccdist=2)
#plot data
plot(as.im(data, W=square(nrow(data))),
     col=gray(seq(1,0,l=length(unique(c(data))))),
     main="", ribbon=TRUE)

#compact grid data - low entropy
data=matrix(sort(sample(c("a","b","c"), 400, replace=TRUE)), nrow=20)
leibovici(data, cell.size=1, ccdist=1.5)
#plot data
plot(as.im(data, W=square(nrow(data))),
     col=heat.colors(length(unique(c(data)))),
     main="", ribbon=TRUE)

#point data
data=ppp(x=runif(400), y=runif(400), window=square(1), marks=sample(1:4, 400, replace=TRUE))
leibovici(data, ccdist=0.1)
#plot data
plot(data)
#nicer plot
data.cat=data; marks(data.cat)=as.character(marks(data))
plot(data.cat, cols=1:length(unique(marks(data.cat))), main="", pch=16)

}
