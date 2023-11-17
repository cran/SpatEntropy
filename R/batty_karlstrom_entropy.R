###########################################
#THIS FILE CONTAINS
#1) function for partitioning the observation window in subareas
#2) function for computing Batty's entropy
#3) function for computing Karlstrom and Ceccato's entropy
###########################################

###########################################
#1) build area partition

#'Area partition.
#'
#'This function partitions the observation area in a number of sub-areas,
#'and assigns the data points/pixels to the areas. This function is useful either
#'when a random partition wants to be created, or when the user wants to set the
#'area's centroids and is happy with an area tessellation in Voronoi polygons according
#'to the defined centroids.
#'
#'The function is preliminary to the computation of Batty's or Karlstrom and Ceccato's entropy.
#'An event of interest (in the form of a point or binary areal dataset) occurs
#'over an observation area divided into sub-areas. If the partition is random,
#'this function generates the sub-areas by randomly drawing the areas' centroids
#'over the observation window. Then, data points/pixels are assigned to the area with
#'the closest centroid. When data are pixels, each pixel is assigned to an area according to
#'the coordinates of its own centroid.
#'The function also works for non-binary datasets and marked ppp objects.
#'
#' @param data If data are lattice, a data matrix, which can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is a \code{ppp} object.
#' @param   cell.size A single number. If data are lattice, the length of the side of each pixel.
#'   Default to 1. Ignored if data are points.
#' @param G An integer if sub-areas are randomly generated, determining the number \eqn{G} of sub-areas.
#'          Alternatively, a 2-column matrix with the sub-areas centroids' coordinates.
#' @param win Optional, the observation area given as a \code{owin} object. If
#'            data are a point pattern \code{ppp} object, this argument is ignored
#'            and the observation area is extracted from the object. If data are
#'            given as a matrix and the area is not specified, the default is a
#'            rectangle with x range from 0 to the number of columns of the data, and
#'            y range from 0 to the number of rows of the data.
#' @param plotout Logical. Default to \code{TRUE}, produces an informative plot as part of the function output.
#'
#' @return A list with elements:
#'\itemize{
#'   \item `G.pp` a point pattern containing the \eqn{G} areas' centroids
#'   \item `data.assign` a four column matrix, with all pairs of data coordinates and data values
#'   matched to one of the \eqn{G} areas (numbered 1 to \eqn{G}). If the dataset is an unmarked ppp
#'   object, the data category column is a vector of 1s.
#'   }
#'   Moreover, a plot is produced showing the data and the area partition.
#'
#' @examples
#' #LATTICE DATA
#'
#' data=matrix(sort(sample(c("a","b","c"), 100, replace=TRUE)), nrow=10)
#' partition=areapart(data, G=5)
#' partition=areapart(data, G=5, cell.size=2)
#'
#' #providing a pre-fixed area partition
#' data=matrix(sort(sample(c("a","b","c"), 100, replace=TRUE)), nrow=10)
#' win=square(nrow(data))
#' GG=cbind(runif(5, win$xrange[1], win$xrange[2]),
#'          runif(5, win$yrange[1], win$yrange[2]))
#' partition=areapart(data, G=GG)
#'
#' #POINT DATA
#'
#' data=ppp(x=runif(100), y=runif(100), window=square(1))
#' partition=areapart(data, 10)
#'
#' #with marks
#' data=ppp(x=runif(100), y=runif(100), window=square(1),
#'          marks=(sample(c("a","b","c"), 100, replace=TRUE)))
#' GG=cbind(runif(10, data$window$xrange[1], data$window$xrange[2]),
#'          runif(10, data$window$yrange[1], data$window$yrange[2]))
#' partition=areapart(data, G=GG)
#'
#' @export

areapart=function(data, G, cell.size=1, win=NULL, plotout=T){

  if(!is.matrix(data) & !spatstat.geom::is.ppp(data))
    stop("For grid data, please provide the dataset as a matrix;
        for point pattern data, please provide the dataset as a ppp object")

  if(is.matrix(data))
  {
    ncl=ncol(data); nrw=nrow(data)

    if(length(cell.size)==1) {
      xsize=cell.size
      ysize=cell.size} else {
        xsize=cell.size[1]
        ysize=cell.size[2]}

    if(is.null(win))
      W=spatstat.geom::owin(xrange=c(0, ncl*xsize), yrange=c(0,nrw*ysize)) else
        W=win
      xx.c=seq(W$xrange[1]+xsize/2, W$xrange[1]+(ncl*xsize-xsize/2), length.out=ncl)
      yy.c=(seq(W$yrange[1]+ysize/2, W$yrange[1]+(nrw*ysize-ysize/2), length.out=nrw))
      cooords=expand.grid(yy.c, xx.c)
      suppressWarnings({data.pp=spatstat.geom::ppp(x=cooords[,2], y=cooords[,1], window=W)})
      spatstat.geom::marks(data.pp)=c(data)#[!is.na(c(data))]
    if(plotout==T)
    spatstat.geom::plot.im(spatstat.geom::as.im(data, W=W), main="",
         col=grDevices::gray(seq(0.9,0.1,l=length(unique(c(data)[!is.na(c(data))])))),
         ribbon=F)

  }

  if (spatstat.geom::is.ppp(data)) {
    W=data$window
    data.pp=data
    if(is.null(spatstat.geom::marks(data))) spatstat.geom::marks(data.pp)=rep(1, spatstat.geom::npoints(data))
    if(plotout==T)
    spatstat.geom::plot.ppp(data, pch=16, cex=.5, main="", legend=F)
  }

  if(length(G)==1){
    #areas' centroids
    part.pp=spatstat.random::runifpoint(G, W)
    part.coord=cbind(x=part.pp$x, y=part.pp$y, id=1:G)
  } else {
    if(min(G[,1])<W$xrange[1]|max(G[,1])>W$xrange[2]|
       min(G[,2])<W$yrange[1]|max(G[,2])>W$yrange[2])
      stop("The given coordinates for the area partition are outside the boundaries of the data observation window")
    part.pp=spatstat.geom::ppp(G[,1], G[,2], W)
    part.coord=cbind(x=part.pp$x, y=part.pp$y, id=1:nrow(G))
  }

  #match data to the nearest area centroid
  near.neigh=spatstat.geom::nncross(data.pp, part.pp)
  data.coord.area=data.frame(data.pp$x, data.pp$y, data.pp$marks, near.neigh$which)
  colnames(data.coord.area)=c("x", "y", "cat", "area")

  if(plotout==T)
    spatstat.geom::plot.tess(spatstat.geom::dirichlet(part.pp), add=TRUE, border=2)

  return(list(G.pp=part.pp, data.assign=data.coord.area))
}
###########################################

###########################################
#2) batty

#'Batty's entropy.
#'
#'This function computes Batty's spatial entropy, following Batty (1976), see also Altieri et al (2017 and following)
#'(references are under the topic \code{\link{SpatEntropy}}).
#'
#'Batty's spatial entropy measures the heterogeneity in the spatial distribution
#'of a phenomenon of interest, with regard to an area partition. It is high
#'when the phenomenon is equally intense over the sub-areas, and low when
#'it concentrates in one or few sub-areas. This function allows to compute Batty's entropy as
#'\deqn{H_B=\sum p_g \log(T_g/p_g)}
#'where \eqn{p_g} is the probability of occurrence of the phenomenon over sub-area \eqn{g},
#'and \eqn{T_g} is the sub-area size.
#'When data are categorical, the phenomenon of interest corresponds to
#'one category, which must be specified. If data are an unmarked
#'point pattern, a fake mark vector is be created with the same category for all points.
#'For comparison purposes, the relative version of Batty's entropy is also returned, i.e.
#'Batty's entropy divided by its maximum \eqn{\log(\sum T_g)}.
#'Note that when the total observation area is 1, then \eqn{\log(\sum T_g)=0}, therefore
#'in that case during the computation all \eqn{T_g}s are multiplied by 100 and a warning is produced.
#'The function is able to work with grids containing missing data, specified as NA values.
#'All NAs are ignored in the computation.
#'
#' @param data If data are lattice, a data matrix, which can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is a \code{ppp} object.
#' @param category A single value matching the data category of interest for computing Batty's entropy.
#'  Default to 1. If the dataset is an unmarked point pattern, this argument must not be changed from the default.
#'  In the plot, only data belonging to the selected category are displayed.
#' @param cell.size A single number or a vector of length two, only needed if data are lattice. It gives the length of the side of each pixel;
#' if the pixel is rectangular, the first number gives the horizontal side and the second number gives the vertical side. Default to 1. Ignored if data are points.
#' @param partition Input defining the partition into subareas. If an integer, it defines the
#' number of sub-areas that are randomly generated by [areapart]; if a two column matrix
#' with coordinates, they are the centroids of the subareas built by [areapart]. Alternatively,
#' it can be the output of [areapart], a \code{tess} object built by the user, a \code{list} object
#' with arguments \code{tiles}, i.e. a list of \code{owin} objects defining the partition, and \code{n}, the number of subareas.
#' Lastly, it can be an \code{im} object, i.e. a factor- or character-valued pixel image on the
#' same observation window as the data, so that the partition is defined according to the values
#' of the image.
#' The default option is \code{partition=areapart(data, G=10)}, which generates 10 random sub-areas.
#' @param win Optional, the observation area given as a \code{owin} object. If
#'            data are a point pattern \code{ppp} object, this argument is ignored
#'            and the observation area is extracted from the object. If data are
#'            given as a matrix, the area should be specified; the default is a
#'            rectangle with x range from 0 to the number of columns of the data, and
#'            y range from 0 to the number of rows of the data.
#' @param rescale Logical. Default to \code{TRUE}, checks whether the size of the observation area or of any of the
#' sub-areas is smaller than 1. If so, computational issues may arise due to negative logarithms, therefore the
#' function automatically performs a rescaling of all sub-areas, computes the entropy over the rescaled area and
#' then transforms it back to the entropy of the original dataset. In such case, a warning is produced to make the user
#' aware of such operation.
#' @param plotout Logical. Default to \code{TRUE}, produces an informative plot as part of the function output.
#'
#' @return A list of five elements:
#'\itemize{
#'   \item `batty` Batty's entropy
#'   \item `range` The theoretical range of Batty's entropy
#'   \item `rel.batty` Batty's entropy divided by \eqn{\log(\sum Tg)} for comparison across observation areas.
#'   \item `areas` a dataframe giving, for each sub-area of the partition, the absolute and relative frequency of
#'   the points/pixels of interest, the sub-area size and the intensity defined as \eqn{pg/Tg}
#'   \item `area.tess` a \code{tess} object with the area partition
#'}
#'Moreover, a plot is produced showing the data and the area partition.
#'
#' @examples
#' #LATTICE DATA
#'
#' data=matrix((sample(c("a","b","c"), 100, replace=TRUE)), nrow=10)
#' batty.entropy=batty(data, category="a")
#'
#' #POINT DATA
#'
#' #unmarked pp
#' data=ppp(x=runif(100, 0, 10), y=runif(100, 0, 10), window=square(10))
#' batty.entropy=batty(data)
#'
#' #smaller window so that some areas' size are smaller than 1
#' data=ppp(x=runif(100, 0, 3), y=runif(100, 0, 3), window=square(3))
#' batty.entropy=batty(data)
#'
#' #marked pp
#' data=ppp(x=runif(100, 0, 10), y=runif(100, 0, 10), window=square(10),
#'          marks=(sample(1:5, 100, replace=TRUE)))
#' plot(data) #see ?plot.ppp for options
#' #if you want to compute the entropy on all points
#' batty.entropy=batty(unmark(data))
#' #if you want to compute the entropy on a category, say 3
#' batty.entropy=batty(data, category=3)
#'
#' @export

batty=function(data, category=1, cell.size=1, partition=10, win=NULL,
               rescale=T, plotout=T){

  if(!is.matrix(data) & !spatstat.geom::is.ppp(data))
    stop("For grid data, please provide the dataset as a matrix;
        for point pattern data, please provide the dataset as a ppp object")

  #dichotomize dataset according to the category of interest
  if(spatstat.geom::is.ppp(data) & !spatstat.geom::is.marked(data) & category!=1)
    stop("Since data do not have different categories, please set category to the default 1")
  if(is.matrix(data)) datavec=c(data) else
    if(spatstat.geom::is.marked(data)) datavec=spatstat.geom::marks(data) else
      datavec=rep(1, spatstat.geom::npoints(data))
    if(is.factor(datavec)) datavec=as.character(datavec)
    if(length(which(unique(datavec)==category))==0)
      stop("In the function arguments, please select a category among the ones of the dataset.
      If the data is a marked point pattern and you wish to compute entropy over all points,
      please use unmark(data) as first argument.
           If the point pattern is unmarked, the argument category must be set to 1")
    datavec=as.numeric(datavec==category)
    #datavec[is.na(datavec)]=0

    #create data ppp object
    if(is.matrix(data))
    {
      ncl=ncol(data); nrw=nrow(data)

      if(length(cell.size)==1) {
        xsize=cell.size
        ysize=cell.size} else {
          xsize=cell.size[1]
          ysize=cell.size[2]}

      if(is.null(win))
        W=spatstat.geom::owin(xrange=c(0, ncl*xsize), yrange=c(0,nrw*ysize)) else
        W=win
      xx.c=seq(W$xrange[1]+xsize/2, W$xrange[1]+(ncl*xsize-xsize/2), length.out=ncl)
      yy.c=(seq(W$yrange[1]+ysize/2, W$yrange[1]+(nrw*ysize-ysize/2), length.out=nrw))
      cooords=expand.grid(yy.c, xx.c)
      suppressWarnings({data.pp=spatstat.geom::ppp(x=cooords[,2], y=cooords[,1], window=W)})
      spatstat.geom::marks(data.pp)=datavec#[!is.na(datavec)]
      #spatstat.geom::plot.im(spatstat.geom::as.im(matrix(datavec, nrow(data)), W=W), main="",
      #                       col=grDevices::gray(seq(0.9,0.1,l=length(unique(c(data)[!is.na(c(data))])))),
      #                       ribbon=F)
    }

    if (spatstat.geom::is.ppp(data)) {
      W=data$window
      data.pp=data
      spatstat.geom::marks(data.pp)=datavec[!is.na(datavec)]
      dataplot=spatstat.geom::unmark(data.pp[spatstat.geom::marks(data.pp)==1])
      #spatstat.geom::plot.ppp(dataplot, pch=16, cex=.5, main="")

    }

    if(spatstat.geom::area.owin(W)==1)
      stop("Batty's entropy cannot be computed on areas of total size = 1.
           You should change your measurement unit so that the window size is not 1.")

    #prepare the structure of object area partition as a tessellation
    if(is.numeric(partition) | is.matrix(partition))
    {
      if(is.matrix(data))
        areap=spatstat.geom::dirichlet(areapart(matrix(datavec,nrow(data)), G=partition,
                                       cell.size=cell.size, win=W, plotout=plotout)$G.pp)
      if(spatstat.geom::is.ppp(data))
        areap=spatstat.geom::dirichlet(areapart(dataplot, G=partition,
                                                cell.size=cell.size, win=W, plotout=plotout)$G.pp)
    } else

    if(is.list(partition) & ! spatstat.geom::is.im(partition)) {
          if(names(partition)[1]=="G.pp" & names(partition)[2]=="data.assign")
           { areap=spatstat.geom::dirichlet(partition$G.pp)
           if(plotout==T)
           spatstat.geom::plot.tess(areap, add=T, border=2)}
          if(names(partition)[1]=="tiles" & names(partition)[2]=="n")
            {areap=partition
            if(plotout==T){
             spatstat.geom::plot.owin(areap$tiles[[1]],border=2, add=T)
             for(ll in 2:areap$n) spatstat.geom::plot.owin(areap$tiles[[ll]],border=2, add=T)}}
        } else

    if(spatstat.geom::is.tess(partition)) {
            if(partition$window$xrange[1]!=W$xrange[1] | partition$window$xrange[2]!=W$xrange[2] |
               partition$window$yrange[1]!=W$yrange[1] | partition$window$yrange[2]!=W$yrange[2])
              stop("The given partition is not on the same observation window as the data")
            if(is.null(partition$tiles)) stop("If a tessellation is provided, it should contain tiles")
            areap=partition
            if(plotout==T)
            spatstat.geom::plot.tess(partition, add=T, border=2)
          } else

    if(spatstat.geom::is.im(partition))
          {
            partnames=names(table(partition$v, useNA="no"))
            vvv=length(partnames)
            if(vvv==length(datavec[!is.na(datavec)]))
              stop("You should provide an image with a limited number of categories for data partition")
            maskv = tiles = vector("list", length=vvv)
            for(ii in 1:vvv)
            {
              maskv[[ii]] = as.logical(c(partition$v)==partnames[ii])
              tiles[[ii]] = spatstat.geom::owin(xrange = W$xrange,
                                 yrange = W$yrange,
                                 mask = matrix(maskv[[ii]],
                                               nrow=partition$dim[1]))
            }
            areap = list(tiles = tiles, n = vvv, names = partnames)
            if(plotout==T){
            spatstat.geom::plot.im(partition, main = "",
                                   col = grDevices::gray(seq(0.9, 0.2, l = vvv)))
            if(spatstat.geom::is.ppp(data)) spatstat.geom::plot.ppp(dataplot, add=T, pch=16, cex=.5)
            if(is.matrix(data)) spatstat.geom::plot.ppp(data.pp[spatstat.geom::marks(data.pp)==1], add=T, pch=15, cex=.2)
       }} else stop("Please provide the area partition object in an accepted format, see ?batty.
                      If a tessellation is provided, it should contain tiles")

    n.G=areap$n
    tot.pG=sum(datavec, na.rm=T)
    absg=pg=Tg=lambdag=numeric(n.G)
    for(gg in 1:n.G)
    {
      subd=data.pp[areap$tiles[[gg]]]
      datatab=table(spatstat.geom::marks(subd))
      if(length(datatab[which(names(datatab)==1)])==1)
        {absg[gg]=datatab[which(names(datatab)==1)]
        pg[gg]=absg[gg]/tot.pG}
      Tg[gg]=spatstat.geom::area.owin(areap$tiles[[gg]])
      if(is.na(Tg[gg]))
        Tg[gg]=table(areap$tiles[[gg]]$m)[which(names(table(areap$tiles[[gg]]$m))==T)]*
        spatstat.geom::area.owin(data.pp$window)/(nrow(areap$tiles[[gg]]$m)*ncol(areap$tiles[[gg]]$m))
      lambdag[gg]=pg[gg]/Tg[gg]
    }
    G.count=data.frame(1:n.G, absg, pg, Tg, lambdag)
    colnames(G.count)=c("area.id", "abs.freq", "rel.freq", "area.size", "area.intens")

    #batty
    batty.terms=ifelse(G.count[,3]>0,G.count[,3]*log(Tg/G.count[,3]),0)
    if(min(Tg)<1)
    {
      if(rescale==F)
        warning("Results may be unreliable due to computational issues in taking logarithms
                of the sub-areas size, since there are areas with size < 1.
                We suggest to re-run the function with rescale = T") else
      {cc=1/min(Tg)+1e-02
      resc.Tg=Tg*cc
      batty.terms=ifelse(G.count[,3]>0,G.count[,3]*log(resc.Tg/G.count[,3]),0)
      batty.ent=sum(batty.terms)-log(cc)

      warning("Some sub-areas have size < 1, so they have been internally rescaled to avoid computational issues.
            The entropy in the output refers to the original area scale.")
      }
    } else
    {
      batty.terms=ifelse(G.count[,3]>0,G.count[,3]*log(Tg/G.count[,3]),0)
      batty.ent=sum(batty.terms)
    }

    batty.range=c(max(0,log(min(Tg))), log(sum(Tg)))
    names(batty.range)=c("Min", "Max")

    return(list(batty=batty.ent,
                range=batty.range,
                rel.batty=batty.ent/log(sum(Tg)),
                areas=G.count,
                area.tess=areap))
}



###########################################
#3) karlstrom

#'Karlstrom and Ceccato's entropy.
#'
#'This function computes Karlstrom and Ceccato's spatial entropy for a
#'chosen neighbourhood distance,
#'following Karlstrom and Ceccato (2002), see also Altieri et al (2017) and following works
#'(references are under the topic \code{\link{SpatEntropy}}).
#'
#'Karlstrom and Ceccato's spatial entropy measures the heterogeneity in the spatial distribution
#'of a phenomenon of interest, with regard to an area partition and accounting for the neighbourhood.
#'It is similar to Batty's entropy (see [batty]) discarding the sub-area size,
#'with the difference that the probability of occurrence of the phenomenon over area \eqn{g}
#'is actually a weighted sum of the neighbouring probabilities.
#'\deqn{H_{KC}=\sum p_g \log(1/ \tilde{p}_g)}
#'where \eqn{p_g} is the probability of occurrence of the phenomenon over sub-area \eqn{g},
#'and \eqn{\tilde{p}_g} is the averaged probability over the neighbouring areas (including the g-th area itself).
#'When data are categorical, the phenomenon of interest corresponds to
#'one category, which must be specified. If data are an unmarked
#'point pattern, a fake mark vector is be created with the same category for all points.
#'For comparison purposes, the relative version of Karlstrom and Ceccato's entropy is also returned, i.e.
#'Karlstrom and Ceccato's entropy divided by its maximum log(number of sub-areas).
#'The function is able to work with grids containing missing data, specified as NA values.
#'All NAs are ignored in the computation.
#'
#' @param data If data are lattice, a data matrix, which can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is a \code{ppp} object.
#' @param category A single value matching the data category of interest for computing Batty's entropy.
#'  Default to 1. If the dataset is an unmarked point pattern, this argument must not be changed from the default.
#' @param cell.size A single number or a vector of length two, only needed if data are lattice. It gives the length of the side of each pixel;
#' if the pixel is rectangular, the first number gives the horizontal side and the second number gives the vertical side. Default to 1. Ignored if data are points.
#' @param partition Input defining the partition into subareas. If an integer, it defines the
#' number of sub-areas that are randomly generated by [areapart]; if a two column matrix
#' with coordinates, they are the centroids of the subareas built by [areapart]. Alternatively,
#' it can be the output of [areapart], a \code{tess} object built by the user, a \code{list} object
#' with arguments \code{tiles}, i.e. a list of \code{owin} objects defining the partition, and \code{n}, the number of subareas.
#' Lastly, it can be an \code{im} object, i.e. a factor- or character-valued pixel image on the
#' same observation window as the data, so that the partition is defined according to the values
#' of the image.
#' The default option is \code{partition=areapart(data, G=10)}, which generates 10 random sub-areas.
#' @param win Optional, the observation area given as a \code{owin} object. If
#'            data are a point pattern \code{ppp} object, this argument is ignored
#'            and the observation area is extracted from the object. If data are
#'            given as a matrix and the area is not specified, the default is a
#'            rectangle with x range from 0 to the number of columns of the data, and
#'            y range from 0 to the number of rows of the data.
#' @param neigh A single number. It can be either the number of neighbours for each sub-area
#' (including the area itself). or the Euclidean distance to define which sub-areas are neighbours,
#' based on their centroids. Default to 4 neighbours.
#' @param method Character, it guides the interpretation of \code{neigh}. Either "number" (the default)
#' or "distance".
#' @param plotout Logical. Default to \code{TRUE}, produces an informative plot as part of the function output.
#'
#' @return A list of five elements:
#'\itemize{
#'   \item `karlstrom` Karlstrom and Ceccato's entropy
#'   \item `range` The theoretical range of Karlstrom and Ceccato's entropy
#'   \item `rel.karl` Karlstrom and Ceccato's entropy divided by \eqn{\log(G)} (number og sub-areas) for comparison across observation areas.
#'   \item `areas` a dataframe giving, for each sub-area, the absolute and relative frequency of
#'   the points/pixels of interest, the weighted probabilities of the neighbours and the sub-area size
#'   \item `area.tess` a \code{tess} object with the area partition
#'}
#'Moreover, a plot is produced showing the data and the area partition.
#'
#' @examples
#' #LATTICE DATA
#'
#' data=matrix((sample(c("a","b","c"), 100, replace=TRUE)), nrow=10)
#' KC.entropy=karlstrom(data, category="a")
#' KC.entropy=karlstrom(data, category="a", neigh=3.5, method="distance")
#' ##to plot
#' data.binary=matrix(as.numeric(data=="a"), nrow(data))
#' plot(as.im(data.binary, W=KC.entropy$area.tess$window), main="",
#'      col=grDevices::gray(seq(1,0,l=length(unique(c(data.binary))))), ribbon=FALSE)
#' plot(KC.entropy$area.tess, add=TRUE, border=2)
#'
#' #POINT DATA
#'
#' #unmarked pp
#' data=ppp(x=runif(100, 0, 10), y=runif(100, 0, 10), window=square(10))
#' KC.entropy=karlstrom(data)
#' ##to plot
#' plot(data)
#' plot(KC.entropy$area.tess, add=TRUE, border=2)
#'
#' #marked pp
#' data=ppp(x=runif(100, 0, 10), y=runif(100, 0, 10), window=square(10),
#'          marks=(sample(1:5, 100, replace=TRUE)))
#' #if you want to compute the entropy on all points
#' KC.entropy=karlstrom(unmark(data))
#' #if you want to compute the entropy on a category, say 3
#' KC.entropy=karlstrom(data, category=3)
#' ##to plot using the selected category
#' ind=which(spatstat.geom::marks(data)==3)
#' data.binary=unmark(data[ind])
#' plot(data.binary)
#' plot(KC.entropy$area.tess, add=TRUE, border=2)
#'
#' @export

karlstrom=function(data, category=1, cell.size=1, partition=10, win=NULL,
                   neigh=4, method="number", plotout=T){

  if(!is.matrix(data) & !spatstat.geom::is.ppp(data))
    stop("For grid data, please provide the dataset as a matrix;
        for point pattern data, please provide the dataset as a ppp object")

  #dichotomize dataset according to the category of interest
  if(spatstat.geom::is.ppp(data) & !spatstat.geom::is.marked(data) & category!=1)
    stop("Since data do not have different categories, please set category to the default 1")
  if(is.matrix(data)) datavec=c(data) else
    if(spatstat.geom::is.marked(data)) datavec=spatstat.geom::marks(data) else
      datavec=rep(1, spatstat.geom::npoints(data))
    if(is.factor(datavec)) datavec=as.character(datavec)
    if(length(which(unique(datavec)==category))==0)
      stop("In the function arguments, please select a category among the ones of the dataset.
      If the data is a marked point pattern and you wish to compute entropy over all points,
      please use unmark(data) as first argument.
           If the point pattern is unmarked, the argument category must be set to 1")
    datavec=as.numeric(datavec==category)
    #datavec[is.na(datavec)]=0

    #create data ppp object
    if(is.matrix(data))
    {
      ncl=ncol(data); nrw=nrow(data)

      if(length(cell.size)==1) {
        xsize=cell.size
        ysize=cell.size} else {
          xsize=cell.size[1]
          ysize=cell.size[2]}

      if(is.null(win))
        W=spatstat.geom::owin(xrange=c(0, ncl*xsize), yrange=c(0,nrw*ysize)) else
          W=win
        xx.c=seq(W$xrange[1]+xsize/2, W$xrange[1]+(ncl*xsize-xsize/2), length.out=ncl)
        yy.c=(seq(W$yrange[1]+ysize/2, W$yrange[1]+(nrw*ysize-ysize/2), length.out=nrw))
        cooords=expand.grid(yy.c, xx.c)
        suppressWarnings({data.pp=spatstat.geom::ppp(x=cooords[,2], y=cooords[,1], window=W)})
        spatstat.geom::marks(data.pp)=datavec#[!is.na(datavec)]
        #spatstat.geom::plot.im(spatstat.geom::as.im(matrix(datavec, nrow(data)), W=W), main="",
        #                       col=grDevices::gray(seq(0.9,0.1,l=length(unique(c(data)[!is.na(c(data))])))),
        #                       ribbon=F)
    }

    if (spatstat.geom::is.ppp(data)) {
      W=data$window
      data.pp=data
      spatstat.geom::marks(data.pp)=datavec[!is.na(datavec)]
      dataplot=spatstat.geom::unmark(data.pp[spatstat.geom::marks(data.pp)==1])
      #spatstat.geom::plot.ppp(dataplot, pch=16, cex=.5, main="")

    }

    #prepare the structure of object area partition as a tessellation
    if(is.numeric(partition) | is.matrix(partition))
    {
      if(is.matrix(data))
        areap=spatstat.geom::dirichlet(areapart(matrix(datavec,nrow(data)), G=partition,
                                                cell.size=cell.size, win=W, plotout=plotout)$G.pp)
      if(spatstat.geom::is.ppp(data))
        areap=spatstat.geom::dirichlet(areapart(dataplot, G=partition,
                                                cell.size=cell.size, win=W, plotout=plotout)$G.pp)
    } else

      if(is.list(partition) & !spatstat.geom::is.im(partition)) {
        if(names(partition)[1]=="G.pp" & names(partition)[2]=="data.assign")
        { areap=spatstat.geom::dirichlet(partition$G.pp)
        if(plotout==T)
        spatstat.geom::plot.tess(areap, add=T, border=2)}
        if(names(partition)[1]=="tiles" & names(partition)[2]=="n")
        {areap=partition
        if(plotout==T){
        spatstat.geom::plot.owin(areap$tiles[[1]],border=2, add=T)
        for(ll in 2:areap$n) spatstat.geom::plot.owin(areap$tiles[[ll]],border=2, add=T)}}
      } else

        if(spatstat.geom::is.tess(partition)) {
          if(partition$window$xrange[1]!=W$xrange[1] | partition$window$xrange[2]!=W$xrange[2] |
             partition$window$yrange[1]!=W$yrange[1] | partition$window$yrange[2]!=W$yrange[2])
            stop("The given partition is not on the same observation window as the data")
          if(is.null(partition$tiles)) stop("If a tessellation is provided, it should contain tiles")
          areap=partition
          if(plotout==T)
          spatstat.geom::plot.tess(partition, add=T, border=2)
        } else

          if(spatstat.geom::is.im(partition))
          {
            partnames=names(table(partition$v, useNA="no"))
            vvv=length(partnames)
            if(vvv==length(datavec[!is.na(datavec)]))
              stop("You should provide an image with a limited number of categories for data partition")
            maskv = tiles = vector("list", length=vvv)
            for(ii in 1:vvv)
            {
              maskv[[ii]] = as.logical(c(partition$v)==partnames[ii])
              tiles[[ii]] = spatstat.geom::owin(xrange = W$xrange,
                                 yrange = W$yrange,
                                 mask = matrix(maskv[[ii]],
                                               nrow=partition$dim[1]))
            }
            areap = list(tiles = tiles, n = vvv, names = partnames)
            if(plotout==T){
            spatstat.geom::plot.im(partition, main = "",
                                   col = grDevices::gray(seq(0.9, 0.2, l = vvv)))
            if(spatstat.geom::is.ppp(data)) spatstat.geom::plot.ppp(dataplot, add=T, pch=16, cex=.5)
            if(is.matrix(data)) spatstat.geom::plot.ppp(data.pp[spatstat.geom::marks(data.pp)==1], add=T, pch=15, cex=.2)
          }} else stop("Please provide the area partition object in an accepted format, see ?karlstrom
                      If a tessellation is provided, it should contain tiles")

    n.G=areap$n

    #neighbourhood
    enclW=spatstat.geom::owin(xrange=W$xrange, yrange=W$yrange)
    centroids=matrix(unlist(lapply(areap$tiles, spatstat.geom::centroid.owin)), byrow=T, ncol=2)
    end=spatstat.geom::ppp(centroids[,1], centroids[,2], window=enclW)
    maxdist=sqrt(diff(data.pp$window$xrange)^2+diff(data.pp$window$yrange)^2)
    mindist=min(spatstat.geom::nndist(end))
    if (method=="number" & neigh%%1!=0) stop("If method=number, neigh must be an integer")
    if (method=="number" & neigh>n.G) stop("The number of neighbours cannot exceed the number of sub-areas")
    if (method=="distance" & neigh>=maxdist)
      warning("The chosen neighbourhood distance is larger than the observation area.
        All areas will be neighbours of all other areas, i.e. all ptilde will be equal")
    if (method=="distance" & neigh<mindist)
      warning("The chosen neighbourhood distance is smaller than the minimum distance between areas.
        All areas will have 0 neighbours, i.e. all ptildeg will be equal to pg")

    neigh.indlist=vector("list", n.G)
    for(gg in 1:n.G)
    {
      start=spatstat.geom::ppp(centroids[gg,1], centroids[gg,2], window=enclW)
      cdist=spatstat.geom::crossdist(start,end)
      if (method=="number") neigh.indlist[[gg]]=which(cdist<=sort(cdist)[neigh]) else
        if (method=="distance")  neigh.indlist[[gg]]=which(cdist<=neigh) else
          stop("Method should be set to either number or distance. If method=number, neigh must be integer.")
    }

    tot.pG=sum(datavec, na.rm=T)
    absg=pg=Tg=ptildeg=numeric(n.G)
    for(gg in 1:n.G)
    {
      subd=data.pp[areap$tiles[[gg]]]
      datatab=table(spatstat.geom::marks(subd))
      if(length(datatab[which(names(datatab)==1)])==1)
        {absg[gg]=datatab[which(names(datatab)==1)]
        pg[gg]=absg[gg]/tot.pG}
      Tg[gg]=spatstat.geom::area.owin(areap$tiles[[gg]])
    }
    for(gg in 1:n.G)
      ptildeg[gg]=mean(pg[neigh.indlist[[gg]]])
    G.count=data.frame(1:n.G, absg, pg, ptildeg, Tg)
    colnames(G.count)=c("area.id", "abs.freq", "rel.freq", "neigh.mean","area.size")

    #karlstrom
    karl.terms=ifelse(G.count[,3]>0&G.count[,4]>0,G.count[,3]*log(1/G.count[,4]),0)
    karl.ent=sum(karl.terms)
    if(karl.ent>log(n.G)) karl.ent=log(n.G)-1e-05
    karl.range=c(0, log(n.G))
    names(karl.range)=c("Min", "Max")

    return(list(karlstrom=karl.ent,
                range=karl.range,
                rel.karl=karl.ent/log(n.G),
                areas=G.count,
                area.tess=areap))
}

#' @rdname karlstrom
#' @export
battyLISA <- karlstrom
