###########################################
#THIS FILE CONTAINS
#1) data description for rainforest tree data 1 and 2 (point patterns)
#2) data description for tree data covariates (images)
#3) data description for turin urban data (lattice)
#4) data description for bologna urban data (lattice)
###########################################


###########################################
#1) tree

#'Rainforest tree data 1.
#'
#'A marked point pattern dataset about four rainforest tree species:
#'Acalypha diversifolia, Chamguava schippii,
#'Inga pezizifera and Rinorea sylvatica
#'
#'This dataset documents the presence of tree species over Barro Colorado Island, Panama.
#'Barro Colorado Island has been the focus of intensive research on lowland tropical
#'rainforest since 1923 (http://www.ctfs.si.edu). Research identified several tree species
#'over a rectangular observation window of size 1000x500 metres; the tree species
#'constitute the point data categorical mark. This dataset presents 4 species with
#'different spatial configurations: Acalypha diversifolia, Chamguava schippii,
#'Inga pezizifera and Rinorea sylvatica. The overall dataset has a total number of 7251 points.
#'The dataset is analyzed with spatial entropy measures in Altieri et al (2018) (references can be
#'found at \code{SpatEntropy}).
#'
#'@format A \code{ppp} object (see package \code{spatstat}) with 7251 points, containing:
#'\describe{
#'   \item{window}{An object of type \code{owin} (see package \code{spatstat}), the 1000x500 metres observation area}
#'   \item{x, y}{Numeric vectors with points' coordinates}
#'   \item{marks}{A character vector matching the tree species to the data points}
#' }
#' @source http://www.ctfs.si.edu
#'
#' @examples
#' data(raintrees)
#' #plot(raintrees, main="", pch=16, cols=1:4)
#'
#' #shannon's entropy of the four trees
#' shannon(raintrees)
#'
#' #shannon's entropy of Z (tree pairs)
#' shannonZ(raintrees)
#'
#' #leibovici's entropy
#' raintrees$window #to check size and unit of measurement
#' #example run on a subset of the data to speed up computations
#' subdata=raintrees[owin(c(0,200),c(0,100))]; plot(subdata)
#' outp=leibovici(subdata, ccdist=10, verbose=TRUE)
#' #do not worry about warnings like "data contain duplicated points": since
#' #coordinates are rounded to the first decimal place, it looks like
#' #some trees are overlapping when they are just very close.
#' #Entropy computation works properly anyway
#'
#' #altieri's entropy
#' #example run on a subset of the data to speed up computations
#' subdata=raintrees[owin(c(0,200),c(0,100))]; plot(subdata)
#' outp=altieri(subdata, distbreak=c(1,2,5,10), verbose=TRUE)
#'
#' #batty's entropy
#' #on all points, with a random partition in 10 sub-areas
#' batty.ent=batty(unmark(raintrees), partition=10)
#' #plot with partition
#' #plot(unmark(raintrees), pch=16, cex=0.6, main="")
#' #plot(batty.ent$area.tess, add=TRUE, border=2, lwd=2)
#' #on a specific tree species, with a random partition in 6 sub-areas
#' unique(marks(raintrees)) #to check the species' names
#' #plot(split.ppp(raintrees), main="") #to plot by species
#' batty.ent=batty(raintrees, category="cha2sc", partition=6)
#' #plot with partition
#' #plot(split.ppp(raintrees)$cha2sc, pch=16, cex=0.6, main="")
#' #plot(batty.ent$area.tess, add=TRUE, border=2)
#'
#' #batty's entropy with a partition based on the covariate,
#' #exploiting spatstat functions
#' data(raintreesCOV)
#' #plot(raintreesCOV$grad, main="", col=gray(seq(1,0,l=100)))
#' data=split.ppp(raintrees)$acaldi
#' #plot(data, add=TRUE, pch=16, cex=0.6, main="")
#' #discretize the covariate
#' slopecut=cut(raintreesCOV$grad,
#'                 breaks = quantile(raintreesCOV$grad, probs = (0:4)/4),
#'                 labels = 1:4)
#' maskv=tiles=list()
#' for(ii in 1:nlevels(slopecut))
#' {
#' maskv[[ii]]=as.logical(c(slopecut$v)==levels(slopecut)[ii])
#' tiles[[ii]]=owin(xrange=data$window$xrange,
#'                  yrange=data$window$yrange,
#'                  mask=matrix(maskv[[ii]],nrow(slopecut$v)))
#' }
#' slopetess=list(tiles=tiles, n=nlevels(slopecut))
#' #plot(slopecut, main = "", col=gray(seq(1,0.4,l=4)))
#' #plot(data, add=TRUE, pch=16, cex=0.6, main="", col=1)
#'
#' batty(data, partition=slopetess)
#'
#' #karlstrom and ceccato's entropy
#' #on a specific tree species, with a random partition in 6 sub-areas
#' unique(marks(raintrees)) #to check the species' names
#' #plot(split.ppp(raintrees), main="") #to plot by species
#' KC.ent=karlstrom(raintrees, category="rinosy", partition=6, neigh=3)
#' #plot with partition
#' #plot(split.ppp(raintrees)$rinosy, pch=16, cex=0.6, main="")
#' #plot(KC.ent$area.tess, add=TRUE, border=2)

"raintrees"
###########################################

#'Rainforest tree data 2.
#'
#'A marked point pattern dataset about four rainforest tree species:
#'Astronium graveolens, Beilschmiedia pendula, Heisteria concinna and Inga sapindoides.
#'
#'This dataset documents the presence of tree species over Barro Colorado Island, Panama.
#'Barro Colorado Island has been the focus of intensive research on lowland tropical
#'rainforest since 1923 (http://www.ctfs.si.edu). Research identified several tree species
#'over a rectangular observation window of size 1000x500 metres; the tree species
#'constitute the point data categorical mark. This dataset presents 4 species with
#'different spatial configurations: Astronium graveolens, Beilschmiedia pendula, Heisteria concinna and Inga sapindoides. The overall dataset has a total number of 5639 points.
#'
#'@format A \code{ppp} object (see package \code{spatstat}) with 5639 points, containing:
#'\describe{
#'   \item{window}{An object of type \code{owin} (see package \code{spatstat}), the 1000x500 metres observation area}
#'   \item{x, y}{Numeric vectors with points' coordinates}
#'   \item{marks}{A character vector matching the tree species to the data points}
#' }
#' @source http://www.ctfs.si.edu
#'
#' @examples
#' data(raintrees2)
#' #plot(raintrees2, main="", pch=16, cols=1:4)
#'
#' #shannon's entropy of the four trees
#' shannon(raintrees2)
#'
#' #shannon's entropy of Z (tree pairs)
#' shannonZ(raintrees2)
#'
#' #leibovici's entropy
#' raintrees2$window #to check size and unit of measurement
#' #example run on a subset of the data to speed up computations
#' subdata=raintrees2[owin(c(0,200),c(0,100))]; plot(subdata)
#' outp=leibovici(subdata, ccdist=10, verbose=TRUE)
#' #do not worry about warnings like "data contain duplicated points": since
#' #coordinates are rounded to the first decimal place, it looks like
#' #some trees are overlapping when they are just very close.
#' #Entropy computation works properly anyway
#'
#' #altieri's entropy
#' #example run on a subset of the data to speed up computations
#' subdata=raintrees2[owin(c(0,200),c(0,100))]; plot(subdata)
#' outp=altieri(subdata, distbreak=c(1,2,5,10), verbose=TRUE)
#'
#' #batty's entropy
#' #on all points, with a random partition in 10 sub-areas
#' batty.ent=batty(unmark(raintrees2), partition=10)
#' #plot with partition
#' #plot(unmark(raintrees2), pch=16, cex=0.6, main="")
#' #plot(batty.ent$area.tess, add=TRUE, border=2, lwd=2)
#' #on a specific tree species, with a random partition in 6 sub-areas
#' unique(marks(raintrees2)) #to check the species' names
#' #plot(split.ppp(raintrees2), main="") #to plot by species
#' batty.ent=batty(raintrees2, category=levels(marks(raintrees2))[1], partition=6)
#' #plot with partition
#' #plot(split.ppp(raintrees2)[[1]], pch=16, cex=0.6, main="")
#' #plot(batty.ent$area.tess, add=TRUE, border=2)
#'
#' #batty's entropy with a partition based on the covariate,
#' #exploiting spatstat functions
#' data(raintreesCOV)
#' #plot(raintreesCOV$grad, main="", col=gray(seq(1,0,l=100)))
#' data=split.ppp(raintrees2)[[1]]
#' #plot(data, add=TRUE, pch=16, cex=0.6, main="")
#' #discretize the covariate
#' slopecut=cut(raintreesCOV$grad,
#'                 breaks = quantile(raintreesCOV$grad, probs = (0:4)/4),
#'                 labels = 1:4)
#' maskv=tiles=list()
#' for(ii in 1:nlevels(slopecut))
#' {
#' maskv[[ii]]=as.logical(c(slopecut$v)==levels(slopecut)[ii])
#' tiles[[ii]]=owin(xrange=data$window$xrange,
#'                  yrange=data$window$yrange,
#'                  mask=matrix(maskv[[ii]],nrow(slopecut$v)))
#' }
#' slopetess=list(tiles=tiles, n=nlevels(slopecut))
#' #plot(slopecut, main = "", col=gray(seq(1,0.4,l=4)))
#' #plot(data, add=TRUE, pch=16, cex=0.6, main="", col=1)
#'
#' batty(data, partition=slopetess)
#'
#' #karlstrom and ceccato's entropy
#' #on a specific tree species, with a random partition in 6 sub-areas
#' unique(marks(raintrees2)) #to check the species' names
#' #plot(split.ppp(raintrees2), main="") #to plot by species
#' KC.ent=karlstrom(raintrees2, category=levels(marks(raintrees2))[2], partition=6, neigh=3)
#' #plot with partition
#' #plot(split.ppp(raintrees2)[[2]], pch=16, cex=0.6, main="")
#' #plot(KC.ent$area.tess, add=TRUE, border=2)

"raintrees2"
###########################################



###########################################
#2) tree COV

#'Covariates for the rainforest tree data.
#'
#'A list of two pixel images with covariates altitude and soil slope for the rainforest tree data 1 and 2,
#'i.e. [raintrees] and [raintrees2]
#'
#'For details of the point datasets, see [raintrees] and [raintrees2]. This accompanying dataset
#'gives information about the elevation in the study region.
#'It is a list containing two pixel images, elev (elevation in metres)
#'and grad (norm of elevation gradient). These pixel images are objects of class \code{im}.
#'Covariate values are continuous. Once discretized as wished, they can turn into categorical datasets
#'for the computation of all entropy measures. Moreover, they can be used
#'to build sensible sub-areas for Batty's and Karlstrom and Ceccato's entropies
#'(see the examples).
#'
#'@format A \code{list} of two elements:
#'\describe{
#'   \item{elev}{An object of type \code{im} (see package \code{spatstat}), the soil elevation}
#'   \item{grad}{An object of type \code{im} (see package \code{spatstat}), the soil slope (gradient of elevation)}
#' }
#' @source http://www.ctfs.si.edu
#'
#' @examples
#' data(raintreesCOV)
#' plot(raintreesCOV, main="")
#'


"raintreesCOV"
###########################################

###########################################
#3) turin

#'Turin urban data.
#'
#'A lattice dataset with Turin's Urban Morphological Zones (UMZ, see EEA, 2011).
#'
#'This raster/pixel/lattice dataset comes from the EU CORINE Land Cover project (EEA, 2011) and is dated 2011.
#'It is the result of classifying the original
#'land cover data into urbanised and non-urbanised zones, known as 'Urban
#'Morphological Zones' (UMZ, see EEA, 2011). UMZ data are useful to identify shapes
#'and patterns of urban areas, and thus to detect what is known as urban sprawl.
#'Turin's metropolitan area is extracted from the European dataset and is
#'composed by the municipality of Turin and the surrounding municipalities: Beinasco, Venaria Reale,
#'San Mauro Torinese, Grugliasco, Borgaro Torinese, Collegno, Pecetto Torinese, Pino Torinese,
#'Moncalieri, Nichelino, Settimo Torinese, Baldissero Torinese, Rivoli, Orbassano.
#'The dataset is made of 111x113 pixels of size 250x250 metres.
#'
#'@format A `matrix` with 111 rows and 113 columns. Values are either 0 (non-urban) or 1 (urban). Pixels
#'outside the administrative borders are classified as NA.
#' @source EEA (2011). Corine land cover 2000 raster data. Technical Report, downloadable at
#' http://www.eea.europa.eu/data-and-maps/ data/corine-land-cover-2000-raster-1
#' @examples
#' data(turin)
#' #plot(as.im(turin), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#'
#' #shannon's entropy
#' shannon(turin)
#'
#' #shannon's entropy of Z (urban/non-urban pairs)
#' shannonZ(turin)
#'
#' #oneill's entropy
#' oneill(turin)
#'
#' #leibovici's entropy only on Collegno's municipality
#' data(turinTess)
#' cell.size=250; ncl=ncol(turin); nrw=nrow(turin)
#' coords=expand.grid(rev(seq(cell.size/2, (nrw*cell.size-cell.size/2), l=nrw)),
#'                    seq(cell.size/2, (ncl*cell.size-cell.size/2), l=ncl))
#' data.pp=ppp(x=coords[which(!is.na(c(turin))),2],
#'             y=coords[which(!is.na(c(turin))),1],
#'             window=owin(xrange=c(0, ncl*cell.size), yrange=c(0,nrw*cell.size)),
#'             marks=c(turin)[which(!is.na(c(turin)))])
#' data=data.pp[turinTess$tiles[[which(turinTess$names=="Collegno")]]]
#' #plot(data, pch=16, cex=0.4)
#' outp=leibovici(data, cell.size=250, ccdist=400, verbose=TRUE)
#'
#' #altieri's entropy only on Collegno's municipality
#' outp=altieri(data, cell.size=250, distbreak=c(cell.size, 2*cell.size), verbose=TRUE)
#'
#' #batty's entropy
#' #on all points, with a random partition in 10 sub-areas
#' batty.ent=batty(turin, cell.size=250, partition=10)
#' #plot with partition
#' data(turinW)
#' #plot(as.im(turin, W=turinW), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' #plot(batty.ent$area.tess, add=TRUE, border=2)
#'
#' #batty's entropy with a partition based on the administrative areas
#' data(turinTess)
#' batty.ent=batty(turin, cell.size=250, partition=turinTess)
#' #plot(as.im(turin, W=turinW), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' #for(i in 1:turinTess$n) plot(turinTess$tiles[[i]], add=TRUE, border=2)
#'
#' #karlstrom and ceccato's entropy
#' data(turinW)
#' KC.ent=karlstrom(turin, cell.size=250, partition=15, neigh=3)
#' #plot with partition
#' #plot(as.im(turin, W=turinW), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' #plot(KC.ent$area.tess, add=TRUE, border=2)
#'
#' #karlstrom and ceccato's entropy with a partition based on the administrative areas
#' data(turinTess)
#' KC.ent=karlstrom(turin, cell.size=250, partition=turinTess, neigh=5000, method="distance")
#' #plot(as.im(turin, W=turinW), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' #for(i in 1:turinTess$n) plot(turinTess$tiles[[i]], add=TRUE, border=2)
#'

"turin"

###########################################

#'Observation window for Turin urban data.
#'
#'An \code{owin} object with the city border for the Turin dataset.
#'
#'This observation window is an \code{owin} object created as a binary mask. See \code{?owin} for details.
#'Examples on the usefulness of the window can be found at the topic [turin].
#'
#'@format An \code{owin} object. Units are given in metres; the basic image unit is a 250x250 metres pixels.
#' @source EEA (2011). Corine land cover 2000 raster data. Technical Report, downloadable at
#' http://www.eea.europa.eu/data-and-maps/ data/corine-land-cover-2000-raster-1
#'
#'@examples
#' data(turinW)
#' plot(turinW, col=c("red", "white"), main="")
#' plot(as.im(turin, W=turinW), main="", col=gray(c(0.8,0)), ribbon=FALSE, add=TRUE)
#'
#' #see examples under the topic "turin"

"turinW"

###########################################

#'Municipalities' administrative borders for Turin urban data.
#'
#'City borders of all municipalities included in the Turin dataset, in the format of polyognal windows \code{owin} objects.
#'
#'The object contains a list of 15 observation windows created as \code{owin} objects based on the coordinates of the border polygons, for each municipality.
#'See \code{?owin} for details.
#'The object also contains the names of the municipalities, in Italian.
#'Examples on the usefulness of the administrative borders can be found at the topic [turin].
#'
#'@format A \code{list} of three:
#'\itemize{
#'   \item `tiles` a \code{list} of 15, each element is a \code{owin} object with the
#'   administrative border of one municipality of Turin's dataset
#'   \item `n` the number of municipalities
#'   \item `names` the names of the 15 municipalities, in the same order as the windows
#'}
#' @source EEA (2011). Corine land cover 2000 raster data. Technical Report, downloadable at
#' http://www.eea.europa.eu/data-and-maps/ data/corine-land-cover-2000-raster-1
#'
#'@examples
#' data(turin); data(turinW); data(turinTess)
#' plot(turinW, col=c("black", "white"), main="")
#' plot(turinTess$tiles[[1]],border=2, add=TRUE, lwd=2)
#' for(ll in 2:turinTess$n) plot(turinTess$tiles[[ll]],border=2, add=TRUE, lwd=2)
#'
#' plot(as.im(turin, W=turinW), main="", col=gray(c(0.85,0.4)), ribbon=FALSE)
#' plot(turinTess$tiles[[1]],border=1, add=TRUE, lwd=2)
#' for(ll in 2:turinTess$n) plot(turinTess$tiles[[ll]],border=1, add=TRUE, lwd=2)
#'
#' #see examples under the topic "turin"

"turinTess"

###########################################
#4) bologna

#'Bologna urban data.
#'
#'A lattice dataset with Bologna's Urban Morphological Zones (UMZ, see EEA, 2011).
#'
#'This raster/pixel/lattice dataset comes from the EU CORINE Land Cover project (EEA, 2011) and is dated 2011.
#'It is the result of classifying the original
#'land cover data into urbanised and non-urbanised zones, known as 'Urban
#'Morphological Zones' (UMZ, see EEA, 2011). UMZ data are useful to identify shapes
#'and patterns of urban areas, and thus to detect what is known as urban sprawl.
#'Bologna's metropolitan area is extracted from the European dataset and is
#'composed by the municipality of Bologna and the surrounding municipalities: .
#'The dataset is made of 135x124 pixels of size 250x250 metres.
#'
#'@format A `matrix` with 135 rows and 124 columns. Values are either 0 (non-urban) or 1 (urban). Pixels
#'outside the administrative borders are classified as NA.
#' @source EEA (2011). Corine land cover 2000 raster data. Technical Report, downloadable at
#' http://www.eea.europa.eu/data-and-maps/ data/corine-land-cover-2000-raster-1
#' @examples
#' data(bologna)
#' #plot(as.im(bologna), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#'
#' #shannon's entropy
#' shannon(bologna)
#'
#' #shannon's entropy of Z (urban/non-urban pairs)
#' shannonZ(bologna)
#'
#' #oneill's entropy
#' oneill(bologna)
#'
#' #leibovici's entropy on a subset of the window
#' bolsub=bologna[30:70,45:85]
#' plot(as.im(bolsub), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' leibovici(bolsub, cell.size=250, ccdist=400, verbose=TRUE)
#'
#' #altieri's entropy
#' bolsub=bologna[30:70,45:85]
#' plot(as.im(bolsub), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' altieri(bolsub, cell.size=250, distbreak=c(250, 500), verbose=TRUE)
#'
#' #batty's entropy
#' #on all points, with a random partition in 10 sub-areas
#' batty.ent=batty(bologna, cell.size=250, partition=10, win=bolognaW)
#' #plot with partition
#' data(bolognaW)
#' #plot(as.im(bologna, W=bolognaW), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' #plot(batty.ent$area.tess, add=TRUE, border=2)
#'
#' #batty's entropy with a partition based on the administrative areas
#' data(bolognaTess)
#' batty.ent=batty(bologna, cell.size=250, partition=bolognaTess, win=bolognaW)
#' #plot(as.im(bologna, W=bolognaW), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' #for(i in 1:bolognaTess$n) plot(bolognaTess$tiles[[i]], add=TRUE, border=2)
#'
#' #karlstrom and ceccato's entropy
#' data(bolognaW)
#' KC.ent=karlstrom(bologna, cell.size=250, partition=15, win=bolognaW, neigh=3)
#' #plot with partition
#' #plot(as.im(bologna, W=bolognaW), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' #plot(KC.ent$area.tess, add=TRUE, border=2)
#'
#' #karlstrom and ceccato's entropy with a partition based on the administrative
#' #areas
#' data(bolognaTess)
#' KC.ent=karlstrom(bologna, cell.size=250, partition=bolognaTess, win=bolognaW,
#' neigh=10000, method="distance")
#' #plot(as.im(bologna, W=bolognaW), main="", col=gray(c(0.8,0)), ribbon=FALSE)
#' #for(i in 1:bolognaTess$n) plot(bolognaTess$tiles[[i]], add=TRUE, border=2)
#'

"bologna"

###########################################

#'Observation window for Bologna urban data.
#'
#'An \code{owin} object with the rectangle circumscribing the city border for the Bologna dataset.
#'
#'This observation window is an \code{owin} object with margins given by the CORINE project coordinates.
#'See \code{?owin} for details.
#'Examples on the usefulness of the window can be found at the topic [bologna].
#'
#'@format An \code{owin} object. The basic spatial unit is a 250x250 metres pixels.
#' @source EEA (2011). Corine land cover 2000 raster data. Technical Report, downloadable at
#' http://www.eea.europa.eu/data-and-maps/ data/corine-land-cover-2000-raster-1
#'
#'@examples
#' data(bolognaW)
#' plot(bolognaW, main="")
#' plot(as.im(bologna, W=bolognaW), main="", col=gray(c(0.8,0)), ribbon=FALSE, add=TRUE)
#'
#' #see examples under the topic "bologna"

"bolognaW"

###########################################

#'Municipalities' administrative borders for Bologna urban data.
#'
#'City borders of all municipalities included in the Bologna dataset, in the format of polyognal windows \code{owin} objects.
#'
#'The object contains a list of 11 observation windows created as \code{owin} objects based on the coordinates of the border polygons, for each municipality.
#'See \code{?owin} for details.
#'The object also contains the names of the municipalities, in Italian.
#'Examples on the usefulness of the administrative borders can be found at the topic [bologna].
#'
#'@format A \code{list} of three:
#'\itemize{
#'   \item `tiles` a \code{list} of 11, each element is a \code{owin} object with the
#'   administrative border of one municipality of Bologna's dataset
#'   \item `n` the number of municipalities
#'   \item `names` the names of the 11 municipalities, in the same order as the windows
#'}
#' @source EEA (2011). Corine land cover 2000 raster data. Technical Report, downloadable at
#' http://www.eea.europa.eu/data-and-maps/ data/corine-land-cover-2000-raster-1
#'
#'@examples
#' data(bologna); data(bolognaW); data(bolognaTess)
#' plot(bolognaW, main="")
#' plot(bolognaTess$tiles[[1]],border=2, add=TRUE, lwd=2)
#' for(ll in 2:bolognaTess$n) plot(bolognaTess$tiles[[ll]],border=2, add=TRUE, lwd=2)
#'
#' plot(as.im(bologna, W=bolognaW), main="", col=gray(c(0.85,0.4)), ribbon=FALSE)
#' plot(bolognaTess$tiles[[1]],border=1, add=TRUE, lwd=2)
#' for(ll in 2:bolognaTess$n) plot(bolognaTess$tiles[[ll]],border=1, add=TRUE, lwd=2)
#'
#' #see examples under the topic "bologna"

"bolognaTess"


