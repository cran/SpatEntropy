###########################################
#THIS FILE CONTAINS
#1) function for building Altieri et al spatial entropy measure
###########################################


###########################################
#'Altieri's spatial entropy.
#'
#'This function computes spatial mutual information and spatial residual entropy as in Altieri et al (2017) and following works.
#'References can be found at \code{SpatEntropy}.
#'
#'The computation of Altieri's entropy starts from a point or areal dataset, for which
#'Shannon's entropy of the transformed variable \eqn{Z} (for details see \code{\link{shannonZ}})
#'\deqn{H(Z)=\sum p(z_r)\log(1/p(z_r))}
#'is computed using all
#'possible pairs within the observation area. Then, its two components
#'spatial mutual information
#'\deqn{SMI(Z,W)=\sum p(w_k) \sum p(z_r|w_k)\log(p(z_r|w_k)/p(z_r))}
#' and spatial residual entropy
#'\deqn{H(Z)_W=\sum p(w_k) \sum p(z_r|w_k)\log(1/p(z_r|w_k))}
#'are calculated
#'in order to account for the overall role of space in determining
#'the data heterogeneity. Besides, starting from a partition into distance
#'classes, a list of adjacency matrices is
#'built, which identifies what pairs of units
#'must be considered for each class. Spatial mutual information and spatial residual
#'entropy are split into local terms according to the chosen distance breaks, so that the role of space
#'can be investigated both in absolute and relative terms. In the function output, the relative partial terms
#'are returned so that they sum to 1 for each distance class: e.g. if the relative SPI terms is 0.3 and
#'the relative residual term is 0.7, the interpretation is that, at the specific distance class, 30% of
#'the entropy is due to the role of space as a source of heterogeneity.
#'The function is able to work with lattice data with missing data, as long as they are specified as NAs:
#'missing data are ignored in the computations.
#'The function is able to work with grids containing missing data, specified as NA values.
#'All NAs are ignored in the computation and only couples of non-NA observations are considered.
#'
#' @param data If data are lattice, a data matrix, which can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is a \code{ppp} object.
#' @param cell.size A single number. If data are lattice, the length of the side of each pixel.
#'   Default to 1. Ignored if data are points.
#' @param distbreak Numeric. The chosen distance breaks for selecting pairs of pixels/points within the observation area.
#'   The default option is \code{c(cell.size, 2cell.size)} for lattice data, and \code{c(mindist, 2mindist)} for point data,
#'   where \code{mindist} is one tenth of the maximum distance within the observation area.
#'   Only the internal breaks have to be specified, the first and last break are automatically added as 0 and the maximum distance within the observation area, respectively.
#' @param verbose Logical. If \code{TRUE} an output is printed in order to follow the progress of the work (recommended for large dataset).
#'   Default set to \code{FALSE}.
#'
#' @return A list with elements:
#'\itemize{
#'   \item `distance.breaks` a two column matrix with the lower and upper extreme of each distance class
#'   \item `SPI.terms` the spatial partial information terms
#'   \item `rel.SPI.terms` the relative version of spatial partial information terms (see the details)
#'   \item `RES.terms` the spatial partial residual entropies
#'   \item `rel.RES.terms` the relative version of spatial partial residual entropies (see the details)
#'   \item `SMI` the spatial mutual information
#'   \item `RES` the global residual entropy
#'   \item `ShannonZ` Shannon's entropy of \eqn{Z} in the same format as the output of [shannonZ()]
#'   \item `W.distribution` the spatial weights for each distance range
#'   \item `total.pairs` the total number of pairs over the area (realizations of \eqn{Z})
#'   \item `class.pairs` the number of pairs for each distance range.
#'   \item `cond.Z.distribution` a list with the conditional absolute and relative frequencies of \eqn{Z} for each distance range
#'}
#'
#' @examples
#' #lattice data
#' data=matrix(sample(1:5, 100, replace=TRUE), nrow=10)
#' outp=altieri(data)
#' outp=altieri(data, cell.size=2) #same result
#' outp=altieri(data, cell.size=2, distbreak=c(2, 5))
#' #plot data
#' plot(as.im(data, W=square(nrow(data))),
#'      col=gray(seq(1,0,l=length(unique(c(data))))),
#'      main="", ribbon=TRUE)
#'
#' #lattice data with missing values
#' data=matrix(sample(1:5, 100, replace=TRUE), nrow=10)
#' data=rbind(rep(NA, ncol(data)), data, rep(NA, ncol(data)))
#' outp=altieri(data)
#' #plot data
#' plot(as.im(data, W=square(nrow(data))),
#'      col=topo.colors(length(unique(c(data)[!is.na(c(data))]))),
#'      main="", ribbon=TRUE)
#'
#' #point data
#' data=ppp(x=runif(400), y=runif(400), window=square(1),
#'          marks=(sample(c("a","b","c"), 400, replace=TRUE)))
#' outp=altieri(data)
#' outp=altieri(data, verbose=TRUE)
#' #plot data
#' plot(data, cols=1:length(unique(marks(data))), main="", pch=16)
#' #check what happens for badly specified distance breaks
#' #outp=altieri(data, distbreak=c(1,1.4))
#' #outp=altieri(data, distbreak=c(1,2))
#'
#' @export

altieri=function(data, cell.size=1, distbreak="default", verbose=F)
{
  if(!is.matrix(data) & !spatstat.geom::is.ppp(data))
    stop("For grid data, please provide the dataset as a matrix;
        for point pattern data, please provide the dataset as a marked ppp object")

  if(is.matrix(data)) {
    if(length(distbreak)>1 & distbreak[1]<cell.size) stop("The first distance break is too small for building any couple")
    if(length(distbreak)>1 & distbreak[length(distbreak)]>=sqrt((nrow(data)*cell.size)^2+(ncol(data)*cell.size)^2))
      stop("The last distance break is equal or larger than the maximum distance over the observation area")
  }

  if(spatstat.geom::is.ppp(data)) {
    if(length(distbreak)>1 & distbreak[1]<min(spatstat.geom::nndist(data))) stop("The first distance break is too small for building any couple")
    if(length(distbreak)>1 & distbreak[length(distbreak)]>=sqrt(diff(data$window$xrange)^2+diff(data$window$yrange)^2))
      stop("The last distance break is equal or larger than the maximum distance over the observation area")
  }

  if(is.matrix(data))
  {
    ncl=ncol(data); nrw=nrow(data)
    W=spatstat.geom::owin(xrange=c(0, ncl*cell.size), yrange=c(0,nrw*cell.size))
    xx.c=seq(cell.size/2, (ncl*cell.size-cell.size/2), l=ncl)
    yy.c=rev(seq(cell.size/2, (nrw*cell.size-cell.size/2), l=nrw))
    coords=expand.grid(yy.c, xx.c)
    datavec=c(data)
    ind=which(!is.na(datavec))
    centr.pp=spatstat.geom::ppp(x=coords[ind,2], y=coords[ind,1], window=W)
    datavec=datavec[ind]
    data.tab=table(datavec)
    Xcat=names(data.tab)
    if(length(Xcat)==1)
      warning("There is only one category in the data, so all entropies will be equal to 0. You may stop the function to avoid wasting time.")
    Xcat.proxy=1:length(Xcat)
    datavec.proxy=c()
    for(ii in 1:length(datavec)) datavec.proxy[ii]=Xcat.proxy[which(Xcat==datavec[ii])]
    spatstat.geom::marks(centr.pp)=datavec.proxy

    if(distbreak[1]=="default")
      dbreak=c(0, cell.size, 2*cell.size, sqrt((nrow(data)*cell.size)^2+(ncol(data)*cell.size)^2)) else
        dbreak=c(0, distbreak, sqrt((nrow(data)*cell.size)^2+(ncol(data)*cell.size)^2))

    speedstep=100
    datapairs.list=vector("list", spatstat.geom::npoints(centr.pp)%/%speedstep+
                            as.numeric(spatstat.geom::npoints(centr.pp)%%speedstep>0))
    for(nn in 1:length(datapairs.list))
    {
      datapairs.list[[nn]]=vector("list", length(dbreak)-1)
      for(dd in 1:(length(dbreak)-1)) datapairs.list[[nn]][[dd]]=rep(NA,1)
    }

    cat("Computing the pairwise distances for all observations. This may take some time...", fill=T)
    for(ii in 1:(spatstat.geom::npoints(centr.pp)-1))
    {
      start=centr.pp[ii]
      end=centr.pp[(ii+1):spatstat.geom::npoints(centr.pp)]
      pdis=c(spatstat.geom::crossdist(start,end))

      for(dd in 1:(length(dbreak)-1))
      {
        ind=which(pdis>dbreak[dd] & pdis<=dbreak[dd+1])
        if(length(ind)>0)
          datapairs.list[[ii%/%speedstep+1]][[dd]]=c(datapairs.list[[ii%/%speedstep+1]][[dd]],
                                                     paste0(spatstat.geom::marks(start),spatstat.geom::marks(end[ind])))
      }
      if(verbose) {if(ii%%100==0) print(paste0("Done for ", ii, "/", spatstat.geom::npoints(centr.pp), " (non-NA) observations"))}
    }

    for(nn in 1:length(datapairs.list))
      for(dd in 1:(length(dbreak)-1)) datapairs.list[[nn]][[dd]]=datapairs.list[[nn]][[dd]][-1]

    datapairs=vector("list", length(dbreak)-1)
    for(dd in 1:(length(dbreak)-1))
    {
      datapairs[[dd]]=rep(NA,1)
      for(nn in 1:length(datapairs.list))
        datapairs[[dd]]=c(datapairs[[dd]], datapairs.list[[nn]][[dd]])
    }
    for(dd in 1:(length(dbreak)-1)) datapairs[[dd]]=datapairs[[dd]][-1]

    cat("All done", fill=T)

    shZ=shannonZ(data)
    names(shZ)=c("marginal distribution", "shannZ", "rel.shannZ")
  }

  if (spatstat.geom::is.ppp(data)) {
    datavec=spatstat.geom::marks(data)
    if(is.factor(datavec)) datavec=as.character(datavec)
    data.tab=table(datavec)
    Xcat=names(data.tab)
    if(length(Xcat)==1)
      warning("There is only one category in the data, so all entropies will be equal to 0. You may stop the function to avoid wasting time.")
    Xcat.proxy=1:length(Xcat)
    datavec.proxy=c()
    for(ii in 1:length(datavec)) datavec.proxy[ii]=Xcat.proxy[which(Xcat==datavec[ii])]
    spatstat.geom::marks(data)=datavec.proxy

    if(distbreak[1]=="default")
    {
      mindist=sqrt(diff(data$window$xrange)^2+diff(data$window$yrange)^2)/10
      dbreak=c(0, mindist, 2*mindist, sqrt(diff(data$window$xrange)^2+diff(data$window$yrange)^2))
    } else
    dbreak=c(0, distbreak, sqrt(diff(data$window$xrange)^2+diff(data$window$yrange)^2))

    speedstep=100
    datapairs.list=vector("list", spatstat.geom::npoints(data)%/%speedstep+
                            as.numeric(spatstat.geom::npoints(data)%%speedstep>0))
    for(nn in 1:length(datapairs.list))
    {
      datapairs.list[[nn]]=vector("list", length(dbreak)-1)
      for(dd in 1:(length(dbreak)-1)) datapairs.list[[nn]][[dd]]=rep(NA,1)
    }

    cat("Computing the pairwise distances for all observations. This may take some time...", fill=T)
    for(ii in 1:(spatstat.geom::npoints(data)-1))
    {
      start=data[ii]
      end=data[(ii+1):spatstat.geom::npoints(data)]
      pdis=c(spatstat.geom::crossdist(start,end))
      if(any(pdis==0)) pdis[which(pdis==0)]=1e-04

      for(dd in 1:(length(dbreak)-1))
      {
        ind=which(pdis>dbreak[dd] & pdis<=dbreak[dd+1])
        if(length(ind)>0)
          datapairs.list[[ii%/%speedstep+1]][[dd]]=c(datapairs.list[[ii%/%speedstep+1]][[dd]],
                            paste0(spatstat.geom::marks(start),"-", spatstat.geom::marks(end[ind])))
      }
      if(verbose) {if(ii%%100==0) print(paste0("Done for ", ii, "/", spatstat.geom::npoints(data), " observations"))}
    }

    for(nn in 1:length(datapairs.list))
      for(dd in 1:(length(dbreak)-1)) datapairs.list[[nn]][[dd]]=datapairs.list[[nn]][[dd]][-1]

    datapairs=vector("list", length(dbreak)-1)
    for(dd in 1:(length(dbreak)-1))
    {
      datapairs[[dd]]=rep(NA,1)
      for(nn in 1:length(datapairs.list))
        datapairs[[dd]]=c(datapairs[[dd]], datapairs.list[[nn]][[dd]])
    }
    for(dd in 1:(length(dbreak)-1)) datapairs[[dd]]=datapairs[[dd]][-1]

    cat("All done", fill=T)

    shZ=shannonZ(datavec)
    names(shZ)=c("marginal distribution", "shannZ", "rel.shannZ")
  }

  if(any(lapply(datapairs,length)==0))
  {
    inddrop=which(lapply(datapairs,length)==0)
    datapairs=datapairs[-inddrop]
    indbr=inddrop+1; indbr[which(inddrop==(length(dbreak)-1))]=inddrop
    dbreak=dbreak[-indbr]
    warning(paste0("Class ", inddrop, " dropped because it does not contain pairs. "))
  }
  catnames=c()
  for(ii in 1:length(Xcat)) catnames=c(catnames,paste0(Xcat[ii], "-", Xcat[ii]))
  if(length(Xcat)>1) for(ii in 1:(length(Xcat)-1)) catnames=c(catnames,paste0(Xcat[ii], "-", Xcat[(ii+1): length(Xcat)]))
  catnames.proxy=c()
  for(ii in 1:length(Xcat)) catnames.proxy=c(catnames.proxy,paste0(Xcat.proxy[ii], "-", Xcat.proxy[ii]))
  if(length(Xcat)>1) for(ii in 1:(length(Xcat)-1)) catnames.proxy=c(catnames.proxy,paste0(Xcat.proxy[ii], "-", Xcat.proxy[(ii+1): length(Xcat)]))

  probabilities=vector("list", length(dbreak)-1)
  names(probabilities)=paste0("conditional distribution for distance class ", 1:(length(dbreak)-1),
                              ": [", round(dbreak[1:(length(dbreak)-1)],2), ", ",
                              round(dbreak[2:length(dbreak)],2), "]")

  for(dd in 1:(length(dbreak)-1))
    probabilities[[dd]]=data.frame("pair"=catnames,
                                   "abs.freq"=numeric(length(catnames)),
                                   "rel.freq"=numeric(length(catnames)))

  tabpairs=lapply(datapairs,table)
  inv.tabpairs=tabpairs
  for(dd in 1:(length(dbreak)-1))
  {
    if(length(tabpairs[[dd]])==1) probabilities[[dd]]$abs.freq=tabpairs[[dd]] else
    {
    lll=lapply(strsplit(names(tabpairs[[dd]]),split="-"),rev)
    for(ii in 1:length(lll))
      names(inv.tabpairs[[dd]])[ii]=paste0(lll[[ii]][1], lll[[ii]][2])

    for(ii in 1:length(Xcat))
    {
      quale=which(names(tabpairs[[dd]])==catnames.proxy[ii])
      if(length(quale)==1) probabilities[[dd]]$abs.freq[ii]=tabpairs[[dd]][quale]
    }
    for(ii in (length(Xcat)+1):nrow(probabilities[[dd]]))
    {
      quale=which(names(tabpairs[[dd]])==catnames.proxy[ii]|names(inv.tabpairs[[dd]])==catnames.proxy[ii])
      if(length(quale)>0) probabilities[[dd]]$abs.freq[ii]=sum(tabpairs[[dd]][quale])
    }
    }
    probabilities[[dd]]$rel.freq=probabilities[[dd]]$abs.freq/sum(probabilities[[dd]]$abs.freq)
  }

  ##ingredients:
  #1) Z marginal frequencies
  P.zr=shZ$"marginal distribution"$rel.freq

  #2) classes wk, and W marginal frequencies (relative frequencies of pairs)
  wks=cbind(round(dbreak[1:(length(dbreak)-1)],2), round(dbreak[2:length(dbreak)],2))
  colnames(wks)=c("from", "to"); rownames(wks)=paste("class", 1:(length(dbreak)-1))
  QQ=choose(length(datavec),2)
  Qk=unlist(lapply(datapairs,length))
  names(Qk)=paste0("class [", round(dbreak[1:(length(dbreak)-1)],2), ", ", round(dbreak[2:length(dbreak)],2), "]")
  P.wk=Qk/QQ

  #3) partial mutual: sum pzr|wk * log(pzr|wk / pzr)
  #   partial resid:  sum pzr|wk * log(1 / pzr|wk)
  mut.partial=res.partial=numeric(length(dbreak)-1)
  names(mut.partial)=names(res.partial)=paste0("class [", round(dbreak[1:(length(dbreak)-1)],2), ", ",
                                               round(dbreak[2:length(dbreak)],2), "]")
  for(dd in 1:(length(dbreak)-1))
  {
    cond.probs=probabilities[[dd]]$rel.freq[probabilities[[dd]]$rel.freq>0]
    marg.probs=P.zr[probabilities[[dd]]$rel.freq>0]
    res.partial[dd]=sum(cond.probs*log(1/cond.probs))
    mut.partial[dd]=sum(cond.probs*log(cond.probs/marg.probs))
  }
  res.global=sum(P.wk*res.partial)
  mut.global=sum(P.wk*mut.partial)
  #sum(P.wk*(res.partial + mut.partial)) #check that this equals shZ$shannZ

  return(list("distance.breaks"=wks,
              "SPI.terms"=mut.partial,
              "rel.SPI.terms"=if(length(Xcat)>1)mut.partial/(mut.partial+res.partial) else NA,
              "RES.terms"=res.partial,
              "rel.RES.terms"=if(length(Xcat)>1)res.partial/(mut.partial+res.partial) else NA,
              "SMI"=mut.global, "RES"=res.global,
              "ShannonZ"=shZ,
              "W.distribution"=P.wk, "total.pairs"=QQ,
              "class.pairs"=Qk,
              "cond.Z.distribution"=probabilities))
}

