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
#'
#' @return A list with elements:
#'\itemize{
#'   \item `G.pp` a point pattern containing the \eqn{G} areas' centroids
#'   \item `data.assign` a four column matrix, with all pairs of data coordinates and data values
#'   matched to one of the \eqn{G} areas (numbered 1 to \eqn{G}). If the dataset is an unmarked ppp
#'   object, the data category column is a vector of 1s.
#'   }
#'
#' @examples
#' #LATTICE DATA
#'
#' data=matrix(sort(sample(c("a","b","c"), 100, replace=TRUE)), nrow=10)
#' partition=areapart(data, G=5)
#' partition=areapart(data, G=5, cell.size=2)
#' ##to plot
#' #data
#' plot(as.im(data, W=partition$G.pp$window), main="",
#'      col=gray(seq(1,0,l=length(unique(c(data)[!is.na(c(data))])))))
#' #partition centroids
#' plot(partition$G.pp, add=TRUE, pch=16, col=2)
#' #partition tessellation
#' plot(dirichlet(partition$G.pp), add=TRUE, border=2)
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
#' #to plot
#' plot(data)
#' plot(partition$G.pp, add=TRUE, col=2, pch=16)
#' plot(dirichlet(partition$G.pp), add=TRUE, border=2)
#'
#' #with marks
#' data=ppp(x=runif(100), y=runif(100), window=square(1),
#'          marks=(sample(c("a","b","c"), 100, replace=TRUE)))
#' GG=cbind(runif(10, data$window$xrange[1], data$window$xrange[2]),
#'          runif(10, data$window$yrange[1], data$window$yrange[2]))
#' partition=areapart(data, G=GG)
#'
#' @export

areapart=function(data, G, cell.size=1){

  if(!is.matrix(data) & !spatstat.geom::is.ppp(data))
    stop("For grid data, please provide the dataset as a matrix;
        for point pattern data, please provide the dataset as a ppp object")

  if(is.matrix(data))
  {
    ncl=ncol(data); nrw=nrow(data)
    W=spatstat.geom::owin(xrange=c(0, ncl*cell.size), yrange=c(0,nrw*cell.size))
    xx.c=seq(cell.size/2, (ncl*cell.size-cell.size/2), l=ncl)
    yy.c=rev(seq(cell.size/2, (nrw*cell.size-cell.size/2), l=nrw))
    coords=expand.grid(yy.c, xx.c)
    data.pp=spatstat.geom::ppp(x=coords[,2], y=coords[,1], window=W)
    spatstat.geom::marks(data.pp)=c(data)
  }

  if (spatstat.geom::is.ppp(data)) {
    W=data$window
    data.pp=data
    if(is.null(spatstat.geom::marks(data))) spatstat.geom::marks(data.pp)=rep(1, spatstat.geom::npoints(data))
  }


  if(length(G)==1){
    #areas' centroids
    part.pp=spatstat.core::runifpoint(G, W)
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
#' @param cell.size A single number. If data are lattice, the length of the side of each pixel.
#'   Default to 1. Ignored if data are points.
#' @param partition Input defining the partition into subareas. If an integer, it defines the
#' number of sub-areas that are randomly generated by [areapart]; if a two column matrix
#' with coordinates, they are the centroids of the subareas built by [areapart]. Alternatively,
#' it can be the output of [areapart], or a \code{tess} object built by the user.
#' The default option is \code{partition=areapart(data, G=10)}, which generates 10 random sub-areas.
#'
#' @return A list of four elements:
#'\itemize{
#'   \item `area.tess` a \code{tess} object with the area partition
#'   \item `areas.freq` a dataframe giving, for each sub-area, the absolute and relative frequency of
#'   the points/pixels of interest and the sub-area size
#'     for each sub-area
#'   \item `batty` Batty's entropy
#'   \item `rel.batty` Batty's entropy divided by \eqn{\log(\sum Tg)} for comparison across observation areas.
#'}
#'
#' @examples
#' #LATTICE DATA
#'
#' data=matrix((sample(c("a","b","c"), 100, replace=TRUE)), nrow=10)
#' batty.entropy=batty(data, category="a")
#' ##to plot
#' data.binary=matrix(as.numeric(data=="a"), nrow(data))
#' plot(as.im(data.binary, W=batty.entropy$areas.tess$window), main="",
#'      col=gray(seq(1,0,l=length(unique(c(data.binary))))), ribbon=FALSE)
#' plot(batty.entropy$areas.tess, add=TRUE, border=2)
#'
#' #POINT DATA
#'
#' #unmarked pp
#' data=ppp(x=runif(100, 0, 10), y=runif(100, 0, 10), window=square(10))
#' batty.entropy=batty(data)
#' ##to plot
#' plot(data)
#' plot(batty.entropy$areas.tess, add=TRUE, border=2)
#'
#' #marked pp
#' data=ppp(x=runif(100, 0, 10), y=runif(100, 0, 10), window=square(10),
#'          marks=(sample(1:5, 100, replace=TRUE)))
#' plot(data) #see ?plot.ppp for options
#' #if you want to compute the entropy on all points
#' batty.entropy=batty(unmark(data))
#' #if you want to compute the entropy on a category, say 3
#' batty.entropy=batty(data, category=3)
#' ##to plot using the selected category
#' ind=which(marks(data)==3)
#' data.binary=unmark(data[ind])
#' plot(data.binary)
#' plot(batty.entropy$areas.tess, add=TRUE, border=2)
#'
#' @export

batty=function(data, category=1, cell.size=1, partition=10){

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
      stop("Please choose a category that is present in the dataset.
           If the point pattern is unmarked, category must be set to 1")
    datavec=as.numeric(datavec==category)
    datavec[is.na(datavec)]=0

    #create data ppp object
    if(is.matrix(data))
    {
      ncl=ncol(data); nrw=nrow(data)
      W=spatstat.geom::owin(xrange=c(0, ncl*cell.size), yrange=c(0,nrw*cell.size))
      xx.c=seq(cell.size/2, (ncl*cell.size-cell.size/2), l=ncl)
      yy.c=rev(seq(cell.size/2, (nrw*cell.size-cell.size/2), l=nrw))
      coords=expand.grid(yy.c, xx.c)
      data.pp=spatstat.geom::ppp(x=coords[,2], y=coords[,1], window=W)
      spatstat.geom::marks(data.pp)=datavec
    }
    if (spatstat.geom::is.ppp(data)) {
      W=data$window
      data.pp=data
      spatstat.geom::marks(data.pp)=datavec
    }

    #prepare the structure of object area partition as a tessellation
    if(is.numeric(partition) | is.matrix(partition))
      areap=spatstat.geom::dirichlet(areapart(data, G=partition, cell.size=cell.size)$G.pp) else
        if(is.list(partition)) {
          if(names(partition)[1]=="G.pp" & names(partition)[2]=="data.assign")
            areap=spatstat.geom::dirichlet(partition$G.pp)
          if(names(partition)[1]=="tiles" & names(partition)[2]=="n")
            areap=partition
        } else
          if(spatstat.geom::is.tess(partition)) {
            if(partition$window$xrange[1]!=W$xrange[1] | partition$window$xrange[2]!=W$xrange[2] |
               partition$window$yrange[1]!=W$yrange[1] | partition$window$yrange[2]!=W$yrange[2])
              stop("The given partition is not on the same observation window as the data")
            if(is.null(partition$tiles)) stop("If a tessellation is provided, it should contain tiles")
            areap=partition
          } else stop("please provide the area partition object in an accepted format.
                      If a tessellation is provided, it should contain tiles")



    n.G=areap$n
    tot.pG=sum(datavec)
    pg=Tg=numeric(n.G)
    for(gg in 1:n.G)
    {
      subd=data.pp[areap$tiles[[gg]]]
      datatab=table(spatstat.geom::marks(subd))
      if(length(datatab[which(names(datatab)==1)])==1)
        pg[gg]=datatab[which(names(datatab)==1)]/tot.pG
      Tg[gg]=spatstat.geom::area.owin(areap$tiles[[gg]])
      if(is.na(Tg[gg]))
        Tg[gg]=table(areap$tiles[[gg]]$m)[which(names(table(areap$tiles[[gg]]$m))==T)]*
        spatstat.geom::area.owin(data.pp$window)/(nrow(areap$tiles[[gg]]$m)*ncol(areap$tiles[[gg]]$m))
    }
    G.count=data.frame(1:n.G, pg*tot.pG, pg, Tg)
    colnames(G.count)=c("area.id", "abs.freq", "rel.freq", "area.size")

    if(sum(Tg)==1)
    {
      Tg=Tg*100
      warning("The total observation area is 1, which returns problems in the computation of Batty's entropy, since the maximum is log(1)=0.
      For this reason, during the computation all areas are multiplied by 100." )
    }

    #batty
    batty.terms=ifelse(G.count[,2]>0,G.count[,3]*log(Tg/G.count[,3]),0)
    batty.ent=sum(batty.terms)

    return(list(areas.tess=areap, areas.freq=G.count, batty=batty.ent, rel.batty=batty.ent/log(sum(Tg))))
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
#' @param cell.size A single number. If data are lattice, the length of the side of each pixel.
#'   Default to 1. Ignored if data are points.
#' @param partition Input defining the partition into subareas. If an integer, it defines the
#' number of sub-areas that are randomly generated by [areapart]; if a two column matrix
#' with coordinates, they are the centroids of the subareas built by [areapart]. Alternatively,
#' it can be the output of [areapart], or a \code{tess} object built by the user.
#' The default option is \code{partition=areapart(data, G=10)}, which generates 10 random sub-areas.
#' @param neigh A single number. It can be either the number of neighbours for each sub-area
#' (including the area itself). or the Euclidean distance to define which sub-areas are neighbours,
#' based on their centroids. Default to 4 neighbours.
#' @param method Character, it guides the interpretation of \code{neigh}. Either "number" (the default)
#' or "distance".
#'
#' @return A list of four elements:
#'\itemize{
#'   \item `area.tess` a \code{tess} object with the area partition
#'   \item `areas.freq` a dataframe giving, for each sub-area, the absolute and relative frequency of
#'   the points/pixels of interest, the weighted probabilities of the neighbours and the sub-area size
#'     for each sub-area
#'   \item `karlstrom` Karlstrom and Ceccato's entropy
#'   \item `rel.karl` Karlstrom and Ceccato's entropy divided by \eqn{\log(G)} (number og sub-areas) for comparison across observation areas.
#'}
#'
#' @examples
#' #LATTICE DATA
#'
#' data=matrix((sample(c("a","b","c"), 100, replace=TRUE)), nrow=10)
#' KC.entropy=karlstrom(data, category="a")
#' KC.entropy=karlstrom(data, category="a", neigh=3.5, method="distance")
#' ##to plot
#' data.binary=matrix(as.numeric(data=="a"), nrow(data))
#' plot(as.im(data.binary, W=KC.entropy$areas.tess$window), main="",
#'      col=gray(seq(1,0,l=length(unique(c(data.binary))))), ribbon=FALSE)
#' plot(KC.entropy$areas.tess, add=TRUE, border=2)
#'
#' #POINT DATA
#'
#' #unmarked pp
#' data=ppp(x=runif(100, 0, 10), y=runif(100, 0, 10), window=square(10))
#' KC.entropy=karlstrom(data)
#' ##to plot
#' plot(data)
#' plot(KC.entropy$areas.tess, add=TRUE, border=2)
#'
#' #marked pp
#' data=ppp(x=runif(100, 0, 10), y=runif(100, 0, 10), window=square(10),
#'          marks=(sample(1:5, 100, replace=TRUE)))
#' #if you want to compute the entropy on all points
#' KC.entropy=karlstrom(unmark(data))
#' #if you want to compute the entropy on a category, say 3
#' KC.entropy=karlstrom(data, category=3)
#' ##to plot using the selected category
#' ind=which(marks(data)==3)
#' data.binary=unmark(data[ind])
#' plot(data.binary)
#' plot(KC.entropy$areas.tess, add=TRUE, border=2)
#'
#' @export

karlstrom=function(data, category=1, cell.size=1, partition=10, neigh=4, method="number"){

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
      stop("Please choose a category that is present in the dataset.
           If the point pattern is unmarked, category must be set to 1")
    datavec=as.numeric(datavec==category)
    datavec[is.na(datavec)]=0

    #create data ppp object
    if(is.matrix(data))
    {
      ncl=ncol(data); nrw=nrow(data)
      W=spatstat.geom::owin(xrange=c(0, ncl*cell.size), yrange=c(0,nrw*cell.size))
      xx.c=seq(cell.size/2, (ncl*cell.size-cell.size/2), l=ncl)
      yy.c=rev(seq(cell.size/2, (nrw*cell.size-cell.size/2), l=nrw))
      coords=expand.grid(yy.c, xx.c)
      data.pp=spatstat.geom::ppp(x=coords[,2], y=coords[,1], window=W)
      spatstat.geom::marks(data.pp)=datavec
    }
    if (spatstat.geom::is.ppp(data)) {
      W=data$window
      data.pp=data
      spatstat.geom::marks(data.pp)=datavec
    }

    #prepare the structure of object area partition as a tessellation
    if(is.numeric(partition) | is.matrix(partition))
      areap=spatstat.geom::dirichlet(areapart(data, G=partition, cell.size=cell.size)$G.pp) else
        if(is.list(partition)) {
          if(names(partition)[1]=="G.pp" & names(partition)[2]=="data.assign")
            areap=spatstat.geom::dirichlet(partition$G.pp)
          if(names(partition)[1]=="tiles" & names(partition)[2]=="n")
            areap=partition
        } else
          if(spatstat.geom::is.tess(partition)) {
            if(partition$window$xrange[1]!=W$xrange[1] | partition$window$xrange[2]!=W$xrange[2] |
               partition$window$yrange[1]!=W$yrange[1] | partition$window$yrange[2]!=W$yrange[2])
              stop("The given partition is not on the same observation window as the data")
            areap=partition
          } else stop("please provide the area partition object in an accepted format")

    n.G=areap$n

    #neighbourhood
    centroids=matrix(unlist(lapply(areap$tiles, spatstat.geom::centroid.owin)), byrow=T, ncol=2)
    end=spatstat.geom::ppp(centroids[,1], centroids[,2], window=W)
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
      start=spatstat.geom::ppp(centroids[gg,1], centroids[gg,2], window=W)
      cdist=spatstat.geom::crossdist(start,end)
      if (method=="number") neigh.indlist[[gg]]=which(cdist<=sort(cdist)[neigh]) else
        if (method=="distance")  neigh.indlist[[gg]]=which(cdist<=neigh) else
          stop("Method should be set to either number or distance. If method=number, neigh must be integer.")
    }

    tot.pG=sum(datavec)
    pg=Tg=ptildeg=numeric(n.G)
    for(gg in 1:n.G)
    {
      subd=data.pp[areap$tiles[[gg]]]
      datatab=table(spatstat.geom::marks(subd))
      if(length(datatab[which(names(datatab)==1)])==1)
        pg[gg]=datatab[which(names(datatab)==1)]/tot.pG
      Tg[gg]=spatstat.geom::area.owin(areap$tiles[[gg]])
    }
    for(gg in 1:n.G)
      ptildeg[gg]=mean(pg[neigh.indlist[[gg]]])
    G.count=data.frame(1:n.G, pg*tot.pG, pg, ptildeg, Tg)
    colnames(G.count)=c("area.id", "abs.freq", "rel.freq", "neigh.mean","area.size")

    #karlstrom
    karl.terms=ifelse(G.count[,3]>0&G.count[,4]>0,G.count[,3]*log(1/G.count[,4]),0)
    karl.ent=sum(karl.terms)
    if(karl.ent>log(n.G)) karl.ent=log(n.G)-1e-05

    return(list(areas.tess=areap, areas.freq=G.count,
                karlstrom=karl.ent, rel.karl=karl.ent/log(n.G)))
}
