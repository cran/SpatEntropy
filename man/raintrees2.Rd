% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raintrees2}
\alias{raintrees2}
\title{Rainforest tree data 2.}
\format{
A \code{ppp} object (see package \code{spatstat}) with 5639 points, containing:
\describe{
\item{window}{An object of type \code{owin} (see package \code{spatstat}), the 1000x500 metres observation area}
\item{x, y}{Numeric vectors with points' coordinates}
\item{marks}{A character vector matching the tree species to the data points}
}
}
\source{
http://www.ctfs.si.edu
}
\usage{
raintrees2
}
\description{
A marked point pattern dataset about four rainforest tree species:
Astronium graveolens, Beilschmiedia pendula, Heisteria concinna and Inga sapindoides.
}
\details{
This dataset documents the presence of tree species over Barro Colorado Island, Panama.
Barro Colorado Island has been the focus of intensive research on lowland tropical
rainforest since 1923 (http://www.ctfs.si.edu). Research identified several tree species
over a rectangular observation window of size 1000x500 metres; the tree species
constitute the point data categorical mark. This dataset presents 4 species with
different spatial configurations: Astronium graveolens, Beilschmiedia pendula, Heisteria concinna and Inga sapindoides. The overall dataset has a total number of 5639 points.
}
\examples{
data(raintrees2)
#plot(raintrees2, main="", pch=16, cols=1:4)

#shannon's entropy of the four trees
shannon(raintrees2)

#shannon's entropy of Z (tree pairs)
shannonZ(raintrees2)

#leibovici's entropy
raintrees2$window #to check size and unit of measurement
#example run on a subset of the data to speed up computations
subdata=raintrees2[owin(c(0,200),c(0,100))]; plot(subdata)
outp=leibovici(subdata, ccdist=10, verbose=TRUE)
#do not worry about warnings like "data contain duplicated points": since
#coordinates are rounded to the first decimal place, it looks like
#some trees are overlapping when they are just very close.
#Entropy computation works properly anyway

#altieri's entropy
#example run on a subset of the data to speed up computations
subdata=raintrees2[owin(c(0,200),c(0,100))]; plot(subdata)
outp=altieri(subdata, distbreak=c(1,2,5,10), verbose=TRUE)

#batty's entropy
#on all points, with a random partition in 10 sub-areas
batty.ent=batty(unmark(raintrees2), partition=10)
#plot with partition
#plot(unmark(raintrees2), pch=16, cex=0.6, main="")
#plot(batty.ent$area.tess, add=TRUE, border=2, lwd=2)
#on a specific tree species, with a random partition in 6 sub-areas
unique(marks(raintrees2)) #to check the species' names
#plot(split.ppp(raintrees2), main="") #to plot by species
batty.ent=batty(raintrees2, category=levels(marks(raintrees2))[1], partition=6)
#plot with partition
#plot(split.ppp(raintrees2)[[1]], pch=16, cex=0.6, main="")
#plot(batty.ent$area.tess, add=TRUE, border=2)

#batty's entropy with a partition based on the covariate,
#exploiting spatstat functions
data(raintreesCOV)
#plot(raintreesCOV$grad, main="", col=gray(seq(1,0,l=100)))
data=split.ppp(raintrees2)[[1]]
#plot(data, add=TRUE, pch=16, cex=0.6, main="")
#discretize the covariate
slopecut=cut(raintreesCOV$grad,
                breaks = quantile(raintreesCOV$grad, probs = (0:4)/4),
                labels = 1:4)
maskv=tiles=list()
for(ii in 1:nlevels(slopecut))
{
maskv[[ii]]=as.logical(c(slopecut$v)==levels(slopecut)[ii])
tiles[[ii]]=owin(xrange=data$window$xrange,
                 yrange=data$window$yrange,
                 mask=matrix(maskv[[ii]],nrow(slopecut$v)))
}
slopetess=list(tiles=tiles, n=nlevels(slopecut))
#plot(slopecut, main = "", col=gray(seq(1,0.4,l=4)))
#plot(data, add=TRUE, pch=16, cex=0.6, main="", col=1)

batty(data, partition=slopetess)

#karlstrom and ceccato's entropy
#on a specific tree species, with a random partition in 6 sub-areas
unique(marks(raintrees2)) #to check the species' names
#plot(split.ppp(raintrees2), main="") #to plot by species
KC.ent=karlstrom(raintrees2, category=levels(marks(raintrees2))[2], partition=6, neigh=3)
#plot with partition
#plot(split.ppp(raintrees2)[[2]], pch=16, cex=0.6, main="")
#plot(KC.ent$area.tess, add=TRUE, border=2)
}
\keyword{datasets}
