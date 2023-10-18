###########################################
#THIS FILE CONTAINS
#1) function for o'neill's entropy
#2) function for leibovici's entropy
#2) function for the contagion index
#3) function for parresol and edwards' entropy
###########################################


###########################################
#1) oneill

#'O'Neill's entropy.
#'
#'This function computes O'Neill's entropy for a data matrix (see O'Neill et al, 1988).
#'
#'O'Neill's entropy index is based on the transformed variable \eqn{Z}, identifying couples of realizations
#'of the variable of interest:
#'\deqn{H_O(Z)=\sum p(z_r|C) \log(1/p(z_r|C))}
#'where \eqn{z_r=(x_i, x_{i'})} is a generic couple of realizations of the study variable \eqn{X}.
#'The conditioning on \eqn{C} for grid data means that only contiguous couples are considered, i.e.
#'couples of pixels sharing a border.
#'All contiguous couples of realizations of the variable of interest are counted
#'and their relative frequencies are used to compute the index. The maximum value for O'Neill's entropy is
#'\eqn{\log(I^2)} where \eqn{I} is the number of categories of \eqn{X}. The relative version of O'Neill's entropy
#'is obtained by dividing the entropy value by its maximum, and is useful for comparison across datasets with
#'a different number of categories.
#'The function is able to work with grids containing missing data, specified as NA values.
#'All NAs are ignored in the computation and only couples of non-NA observations are considered.
#'
#' @param data A data matrix, can be numeric, factor, character, ...
#' @param win Optional, an object of class \code{owin}, the observation window for data plotting
#'
#' @return a list of four elements:
#' \itemize{
#' \item `oneill` O'Neill's entropy
#' \item `range` the theoretical range of O'Neill's entropy, from 0 to \eqn{\log(I^2)}
#' \item`rel.oneill` O'Neill's relative entropy
#' \item `probabilities` a table with absolute frequencies and estimated probabilities (relative frequencies) for all couple categories
#' }
#' Moreover, a plot of the dataset is produced.
#'
#' @examples
#' #numeric data, square grid
#' data=matrix(sample(1:5, 100, replace=TRUE), nrow=10)
#' oneill(data)
#' #plot data
#' plot(as.im(data, W=square(nrow(data))),
#'      col=grDevices::gray(seq(1,0,length.out=length(unique(c(data))))),
#'      main="", ribbon=TRUE)
#'
#' #character data, rectangular grid
#' data=matrix(sample(c("a","b","c"), 300, replace=TRUE), nrow=30)
#' oneill(data)
#' #plot data
#' plot(as.im(data, W=owin(xrange=c(0,ncol(data)), yrange=c(0,nrow(data)))),
#'      col=terrain.colors(length(unique(c(data)))),
#'      main="", ribbon=TRUE)
#'
#' #data with missing values
#' data=matrix(sample(1:5, 100, replace=TRUE), nrow=10)
#' data=rbind(rep(NA, ncol(data)), data, rep(NA, ncol(data)))
#' oneill(data)
#'
#' @export

oneill=function(data, win=spatstat.geom::owin(xrange=c(0, ncol(data)), yrange=c(0, nrow(data))))
{
  if(!is.matrix(data))
    stop("This function works for grid data. Please provide the dataset as a matrix, or use the function -leibovici- instead")

  spatstat.geom::plot.im(spatstat.geom::as.im(matrix(c(data), nrow(data)), W=win), main="",
                         col=grDevices::gray(seq(0.9,0.1,l=length(unique(c(data)[!is.na(c(data))])))),
                         ribbon=T)
  Xcat=unique(c(data)[!is.na(c(data))])
  if(length(Xcat)==1) warning("Data must have at least two categories to compute entropy")

  datapairs=c()
  for(i in 1:(nrow(data)-1))
  {
    ind=which(!is.na(data[i,]) & !is.na(data[i+1,]))
    if(length(ind)>0)
      datapairs=c(datapairs, paste0(data[i,ind],"-",data[i+1,ind]))
  }
  for(j in 1:(ncol(data)-1))
  {
    ind=which(!is.na(data[,j]) & !is.na(data[,j+1]))
    if(length(ind)>0)
      datapairs=c(datapairs, paste0(data[ind,j], "-",data[ind,j+1]))
  }
  tabpairs=table(datapairs)

  probabilities=data.frame("couple"=names(tabpairs),
                           "abs.freq"=as.numeric(tabpairs),
                           "rel.freq"=as.numeric(prop.table(tabpairs)))


  terms=ifelse(probabilities$rel.freq!=0, probabilities$rel.freq*log(probabilities$rel.freq), 0)
  oneill=-sum(terms)

  oneill.range=c(0, log(length(Xcat)^2))
  names(oneill.range)=c("Min", "Max")

  return(list(oneill=oneill,
              range=oneill.range,
              rel.oneill=oneill/log(length(Xcat)^2),
              probabilities=probabilities))

}

###########################################
#2) leibovici

#'Leibovici's entropy.
#'
#'This function computes Leibovici's entropy according to a chosen distance \eqn{d}
#'(with O'Neill's entropy as a special case)
#'following Leibovici (2009), see also Altieri et al (2017). References can be found at \code{SpatEntropy}.
#'
#'This index is based on the transformed variable \eqn{Z} identifying couples of realizations
#'of the variable of interest. A distance of interest is fixed, which in the case of O'Neill's
#'entropy is the contiguity, i.e. sharing a border for lattice data. Then, all couples
#'of realizations of the variable of interest lying at a distance smaller or equal to the distance of interest
#'are counted, and their relative frequencies are used to compute the index with the traditional Shannon's
#'formula.
#'#'\deqn{H_L(Z)=\sum p(z_r|d) \log(1/p(z_r|d))}
#'where \eqn{z_r=(x_i, x_{i'})} is a generic couple of realizations of the study variable \eqn{X}.
#'The conditioning on \eqn{d} means that only couples within a predefined distance are considered.
#'The maximum value for Leibovici's entropy is
#'\eqn{\log(I^2)} where \eqn{I} is the number of categories of the study variable \eqn{X}.
#'The relative version of Leibovici's entropy is obtained by dividing the entropy value by its maximum, and is useful for comparison across datasets with
#'a different number of categories.
#'The function is able to work with grids containing missing data, specified as NA values.
#'All NAs are ignored in the computation and only couples of non-NA observations are considered.
#'
#' @param data If data are lattice, a data matrix, which can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is a \code{ppp} object.
#' @param cell.size A single number or a vector of length two, only needed if data are lattice. It gives the length of the side of each pixel;
#' if the pixel is rectangular, the first number gives the horizontal side and the second number gives the vertical side. Default to 1. Ignored if data are points.
#' @param ccdist A single number. The chosen distance for selecting couples of pixels/points within the observation area. Default to \code{cell.size[1]}.
#' @param win Optional, an object of class \code{owin}, the observation window for data plotting
#' @param verbose Logical. If \code{TRUE} an output is printed in order to follow the progress of the work (recommended for large dataset).
#'   Default set to \code{FALSE}.
#'
#' @return a list of four elements:
#' \itemize{
#' \item `leib` Leibovici's entropy
#' \item `range` the theoretical range of Leibovici's entropy, from 0 to \eqn{\log(I^2)}
#' \item `rel.leib` Leibovici's relative entropy
#' \item `probabilities` a table with absolute frequencies and estimated probabilities (relative frequencies) for all couple categories
#' }
#' Moreover, a plot of the dataset is produced. Over the plot, a random point is displayed as a red star, and a circle is plotted
#' around that point. The radius of the circle is set by \code{ccdist}, so that a visual idea is given about the
#' choice of the distance for building co-occurrences.
#'
#' @examples
#'#random grid data - high entropy
#' data=matrix(sample(c("a","b","c"), 400, replace=TRUE), nrow=20)
#' leibovici(data, cell.size=1, ccdist=2)
#' #plot data
#' plot(as.im(data, W=square(nrow(data))),
#'      col=grDevices::gray(seq(1,0,length.out=length(unique(c(data))))),
#'      main="", ribbon=TRUE)
#'
#'#compact grid data - low entropy
#' data=matrix(sort(sample(c("a","b","c"), 400, replace=TRUE)), nrow=20)
#' #Note: with sorted data, only some couple categories will be present
#' leibovici(data, cell.size=1, ccdist=1.5)
#' #plot data
#' plot(as.im(data, W=square(nrow(data))),
#'      col=heat.colors(length(unique(c(data)))),
#'      main="", ribbon=TRUE)
#'
#'#point data
#' data=ppp(x=runif(400), y=runif(400),
#' window=square(1), marks=sample(1:4, 400, replace=TRUE))
#' leibovici(data, ccdist=0.1)
#' #plot data
#' plot(data)
#'
#' @export

leibovici=function(data, cell.size=1, ccdist=cell.size[1], win=NULL, verbose=F)
{
  if(!is.matrix(data) & ! spatstat.geom::is.ppp(data))
    stop("For grid data, please provide the dataset as a matrix;
        for point pattern data, please provide the dataset as a marked ppp object")

  if(is.matrix(data)) {
    if(ccdist<min(cell.size)) stop("The distance of interest is too small for building any couple")
    if(ccdist>=sqrt((nrow(data)*cell.size[length(cell.size)])^2+(ncol(data)*cell.size[1])^2))
      stop("The chosen distance is equal or larger than the maximum distance over the observation area. Maybe you wish to compute the non-spatial Shannon's entropy of Z instead?")
  }

  if(spatstat.geom::is.ppp(data)) {
    if(ccdist<min(spatstat.geom::nndist(data))) stop("The distance of interest is too small for building any couple")
    if(ccdist>=sqrt(diff(data$window$xrange)^2+diff(data$window$yrange)^2))
      stop("The chosen distance is equal or larger than the maximum distance over the observation area. Maybe you wish to compute the non-spatial Shannon's entropy of Z instead?")
  }

  if(is.matrix(data) & length(cell.size)==1 & ccdist==cell.size[1]) {
    outp=oneill(data)
    return(list(probabilities=outp$probabilities, leib=outp$oneill, rel.leib=outp$rel.oneill))
  }

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
    datavec=c(data)
    ind=which(!is.na(datavec))
    centr.pp=spatstat.geom::ppp(x=cooords[ind,2], y=rev(cooords[ind,1]), window=W)
    datavec=datavec[ind]
    data.tab=table(datavec)
    Xcat=names(data.tab)
    if(length(Xcat)==1)
      warning("There is only one category in the dataset, so the entropy will be equal to 0")
    Xcat.proxy=1:length(Xcat)
    datavec.proxy=c()
    for(ii in 1:length(datavec)) datavec.proxy[ii]=Xcat.proxy[which(Xcat==datavec[ii])]
    spatstat.geom::marks(centr.pp)=datavec.proxy

    speedstep=100
    datapairs.list=vector("list", spatstat.geom::npoints(centr.pp)%/%speedstep+
                            as.numeric(spatstat.geom::npoints(centr.pp)%%speedstep>0))
    for(nn in 1:length(datapairs.list)) datapairs.list[[nn]]=rep(NA,1)

    cat("Computing the pairwise distances for all observations. This may take some time...", fill=T)
    for(ii in 1:(spatstat.geom::npoints(centr.pp)-1))
    {
      start=centr.pp[ii]
      possend=centr.pp[-ii]
      end=possend[which(possend$x>=start$x & possend$y<=start$y)]
      pdis=c(spatstat.geom::crossdist(start,end))
      ind=which(pdis<=ccdist)
      if(length(ind)>0)
        datapairs.list[[ii%/%speedstep+1]]=c(datapairs.list[[ii%/%speedstep+1]],
                                             paste0(spatstat.geom::marks(start),"-",spatstat.geom::marks(end[ind])))
      if(verbose) {if(ii%%100==0) print(paste0("Done for ", ii, "/", spatstat.geom::npoints(centr.pp), " (non-NA) observations"))}
    }
    for(nn in 1:length(datapairs.list))
      datapairs.list[[nn]]=datapairs.list[[nn]][-1]
    datapairs=c()
      for(nn in 1:length(datapairs.list))
        datapairs=c(datapairs, datapairs.list[[nn]])

    cat("Done", fill=T)
  }

  if (spatstat.geom::is.ppp(data)) {
    datavec=spatstat.geom::marks(data)
    if(is.factor(datavec)) datavec=as.character(datavec)
    if(is.null(spatstat.geom::marks(data))) datavec=rep(1, data$n)
    data.tab=table(datavec)
    Xcat=names(data.tab)
    if(length(Xcat)==1)
      warning("There is only one category in the dataset, so the entropy will be equal to 0")
    Xcat.proxy=1:length(Xcat)
    datavec.proxy=c()
    for(ii in 1:length(datavec)) datavec.proxy[ii]=Xcat.proxy[which(Xcat==datavec[ii])]
    spatstat.geom::marks(data)=datavec.proxy

    speedstep=100
    datapairs.list=vector("list", spatstat.geom::npoints(data)%/%speedstep+
                            as.numeric(spatstat.geom::npoints(data)%%speedstep>0))
    for(nn in 1:length(datapairs.list))
      datapairs.list[[nn]]=rep(NA,1)

    cat("Computing the pairwise distances for all observations. This may take some time...", fill=T)
    for(ii in 1:(spatstat.geom::npoints(data)-1))
    {
      start=data[ii]
      possend=data[-ii]
      end=possend[which(possend$x>=start$x & possend$y<=start$y)]
      pdis=c(spatstat.geom::crossdist(start,end))
      ind=which(pdis<=ccdist)
      if(length(ind)>0)
        datapairs.list[[ii%/%speedstep+1]]=c(datapairs.list[[ii%/%speedstep+1]],
                                             paste0(spatstat.geom::marks(start),"-",spatstat.geom::marks(end[ind])))

      if(verbose) {if(ii%%100==0) print(paste0("Done for ", ii, "/", spatstat.geom::npoints(data), " observations"))}
    }

    for(nn in 1:length(datapairs.list))
      datapairs.list[[nn]]=datapairs.list[[nn]][-1]

    datapairs=c()
      for(nn in 1:length(datapairs.list))
        datapairs=c(datapairs, datapairs.list[[nn]])

    cat("Done", fill=T)
  }

  tabpairs=table(datapairs)
  catnames=c()
  for(ii in 1:length(Xcat)) catnames=c(catnames,paste0(Xcat[ii],"-", Xcat))
  catnames.proxy=c()
  for(ii in 1:length(Xcat)) catnames.proxy=c(catnames.proxy,paste0(Xcat.proxy[ii], "-",Xcat.proxy))

  for(ii in 1:length(tabpairs)) names(tabpairs)[ii]=catnames[which(catnames.proxy==names(tabpairs)[ii])]
  probabilities=data.frame("couple"=names(tabpairs),
                           "abs.freq"=as.numeric(tabpairs),
                           "rel.freq"=as.numeric(prop.table(tabpairs)))

  terms=ifelse(probabilities$rel.freq!=0, probabilities$rel.freq*log(probabilities$rel.freq), 0)
  leib=-sum(terms)

  leib.range=c(0, log(length(Xcat)^2))
  names(leib.range)=c("Min", "Max")

  if(spatstat.geom::is.ppp(data)){
    data.cat=data; spatstat.geom::marks(data.cat)=as.character(spatstat.geom::marks(data))
    spatstat.geom::plot.ppp(data.cat, cols=(1:(length(unique(spatstat.geom::marks(data.cat)))+1))[-2],
                            main="", pch=16, cex=.5)
    spatstat.geom::plot.ppp(spatstat.geom::ppp(spatstat.geom::centroid.owin(data$window)$x,
                                spatstat.geom::centroid.owin(data$window)$y,
                                window=data$window), add=T, pch="*", cex=1.5, col=2)
    spatstat.geom::plot.owin(spatstat.geom::disc(radius = ccdist,
                                  centre = spatstat.geom::centroid.owin(data$window)),
                             add=T, border=2, lwd=2)
  }
  if(is.matrix(data)){
    if(is.null(win))
      win=spatstat.geom::owin(xrange=c(0, ncol(data)*xsize), yrange=c(0, nrow(data)*ysize))
    spatstat.geom::plot.im(spatstat.geom::as.im(matrix(c(data), nrow(data)), W=win), main="",
                           col=grDevices::gray(seq(0.9,0.1,l=length(unique(c(data)[!is.na(c(data))])))),
                           ribbon=T)
    spatstat.geom::plot.ppp(spatstat.geom::ppp(spatstat.geom::centroid.owin(win)$x,
                                spatstat.geom::centroid.owin(win)$y,
                                window=win), add=T, pch="*", cex=1.5, col=2)
    spatstat.geom::plot.owin(spatstat.geom::disc(radius = ccdist,
                                  centre = spatstat.geom::centroid.owin(win)),
                             add=T, border=2, lwd=2)

  }

  return(list(leib=leib,
              range=leib.range,
              rel.leib=leib/log(length(Xcat)^2),
              probabilities=probabilities))
}

###########################################

###########################################
#3) contagion

#'Li and Reynolds' relative contagion index.
#'
#'This function computes Li and Reynold's contagion index, following Li and Reynolds (1993),
#'starting from a data matrix. References can be found at \code{SpatEntropy}.
#'
#'This index is based on the transformed variable \eqn{Z} identifying couples of realizations
#'of the variable of interest. A distance of interest is fixed: the contagion index is
#'originally thought for areas sharing a border, as O'Neill's entropy. Then, all contiguous couples
#'of realizations of the variable of interest are counted
#'and their relative frequencies are used to compute the index, which is \eqn{1-NO} where \eqn{NO}
#'is the relative version of O'Neill's entropy, i.e. O'Neill's entropy divided by its maximum \eqn{\log(I^2)},
#'\eqn{I} being the number of categories of the variable under study. The relative contagion index ranges
#'from 0 (no contagion, maximum entropy) to 1 (maximum contagion).
#'The function is able to work with grids containing missing data, specified as NA values.
#'All NAs are ignored in the computation and only couples of non-NA observations are considered.
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#' @param win Optional, an object of class \code{owin}, the observation window for data plotting
#'
#' @return a list of two elements:
#' \itemize{
#' \item `contagion` Li and Reynold's relative contagion index
#'   \item `probabilities` a table with absolute frequencies and estimated probabilities (relative frequencies) for all couple categories
#' }
#' Moreover, a plot of the dataset is produced.
#'
#' @examples
#' #numeric data, square grid
#' data=matrix(sample(1:5, 100, replace=TRUE), nrow=10)
#' contagion(data)
#' #plot data
#' plot(as.im(data, W=square(nrow(data))),
#'      col=grDevices::gray(seq(1,0,length.out=length(unique(c(data))))),
#'      main="", ribbon=TRUE)
#'
#' #character data, rectangular grid
#' data=matrix(sample(c("a","b","c"), 300, replace=TRUE), nrow=30)
#' contagion(data)
#' #plot data
#' plot(as.im(data, W=owin(xrange=c(0,ncol(data)), yrange=c(0,nrow(data)))),
#'      col=terrain.colors(length(unique(c(data)))),
#'      main="", ribbon=TRUE)
#'
#' @export

contagion=function(data, win=spatstat.geom::owin(xrange=c(0, ncol(data)), yrange=c(0, nrow(data))))
{
  outp=oneill(data, win)
  return(list(contagion=1-outp$rel.oneill,
              probabilities=outp$probabilities))
}
###########################################

###########################################
#4) parresol

#'Parresol and Edwards' entropy.
#'
#'Compute Parresol and Edwards' entropy, following Parresol and Edwards (2014),
#'starting from data. References can be found at \code{SpatEntropy}.
#'
#'This index is based on the transformed variable \eqn{Z} identifying couples of realizations
#'of the variable of interest. A distance of interest is fixed: Parresol and Edwards' entropy is
#'thought for areas sharing a border, as O'Neill's entropy. All contiguous couples
#'of realizations of the variable of interest are counted
#'and their relative frequencies are used to compute the index, which is the opposite of O'Neill's entropy.
#'The function is able to work with grids containing missing data, specified as NA values.
#'All NAs are ignored in the computation and only couples of non-NA observations are considered.
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#' @param win Optional, an object of class \code{owin}, the observation window for data plotting
#'
#' @return a list of four elements:
#' \itemize{
#' \item `parredw` Parresol and Edwards' entropy
#' \item `range` the theoretical range of Parresol and Edwards' entropy, from \eqn{-\log(I^2)} to 0
#' \item`rel.parredw` Parresol and Edwards' relative entropy (with the same interpretation as O'Neill's relative entropy)
#' \item `probabilities` a table with absolute frequencies and estimated probabilities (relative frequencies) for all couple categories
#' }
#' Moreover, a plot of the dataset is produced.
#'
#' @examples
#' #numeric data, square grid
#' data=matrix(sample(1:5, 100, replace=TRUE), nrow=10)
#' parredw(data)
#' #plot data
#' plot(as.im(data, W=square(nrow(data))),
#'      col=grDevices::gray(seq(1,0,length.out=length(unique(c(data))))),
#'      main="", ribbon=TRUE)
#'
#' #character data, rectangular grid
#' data=matrix(sample(c("a","b","c"), 300, replace=TRUE), nrow=30)
#' parredw(data)
#' #plot data
#' plot(as.im(data, W=owin(xrange=c(0,ncol(data)), yrange=c(0,nrow(data)))),
#'      col=terrain.colors(length(unique(c(data)))),
#'      main="", ribbon=TRUE)
#'
#' @export

parredw=function(data, win=spatstat.geom::owin(xrange=c(0, ncol(data)), yrange=c(0, nrow(data))))
{
  outp=oneill(data, win)
  parr.range=-rev(outp$range)
  names(parr.range)=c("Min", "Max")
  return(list(parredw=-outp$oneill,
              range=parr.range,
              rel.parredw=as.numeric(-outp$oneill/parr.range[1]),
              probabilities=outp$probabilities))
}
###########################################
