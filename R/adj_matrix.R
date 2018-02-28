###########################################
#THIS FILE CONTAINS
#1) function for euclidean distance between two coordinate matrices
#2) function for generating centroids' coordinates for lattice data
#3) function for building the adjacency matrix for a chosen distance range
#4) function for building an adjacency matrix list for multiple distance ranges
###########################################


###########################################
#1)euclidean distances between two coordinate matrices

#'Euclidean distance.
#'
#'\code{euclid_dist} computes the Euclidean distance between all point couples/pairs
#'identified by two coordinate matrices.
#'Useful alone, or together with [adj_mat()].
#'
#'`euclid_dist` needs two matrices listing the coordinates of two sets of points.
#'It computes the Euclidean distance between each point of the first matrix and
#'all points of the second matrix. The default option provides the Euclidean distance
#'between all points in a single set.
#'
#' @param coords1 A two column matrix containing starting coordinates.
#'                Provided by user or as output of [coords_pix()].
#' @param coords2 A two column matrix with containing coordinates, same as \code{coords1} by default.
#'
#' @return An upper-triangular distance matrix: a `nrow(coords1)` x `nrow(coords2)` matrix with Euclidean distances between points.
#'
#' @examples
#' euclid_dist(cbind(runif(10), runif(10)))
#'
#' @export

euclid_dist=function(coords1, coords2=coords1)
{
  distmat=matrix(,nrow(coords1),nrow(coords2))
  for (i in 1:nrow(coords1))
  {
    for(j in 1:nrow(coords2))
    {
      if (coords1[i,1]<coords2[j,1] |
          (coords1[i,1]==coords2[j,1] & coords1[i,2]>coords2[j,2])) {
        a=coords1[i,1]-coords2[j,1]
        b=coords1[i,2]-coords2[j,2]
        if(i<j) distmat[i,j]=sqrt(a^2+b^2) else distmat[j,i]=sqrt(a^2+b^2)}}
  }
  return(distmat)
}
###########################################


###########################################
#2) lattice coordinates generation

#'Pixel coordinates generation.
#'
#'\code{coords_pix} generates the coordinates of the pixel centroids,
#'given an obervation area and the pixel size.
#'The resulting coordinates may be used as arguments of [euclid_dist()].
#'
#'When the Euclidean distance is computed between areas, the most common choice is
#'to use the area centroid as the reference point. This function returns the centroids
#'of a regular lattice. You can either provide the number of pixels along the two cardinal
#'directions (using `ncol` for the \eqn{x} direction and `nrow` for the \eqn{y} direction), or alternatively
#'provide the pixel size along the two directions (rectangular areas are allowed).
#'
#' @param area The observation area, an \code{owin} object, see package `spatstat`
#' @param pixel.xsize A scalar, length of the pixel side along the \eqn{x} axis (unnecessary if \code{ncol} is provided)
#' @param pixel.ysize A scalar, length of the pixel side along the \eqn{y} axis (unnecessary if \code{nrow} is provided)
#' @param nrow An integer, number of pixels along the \eqn{y} axis (unnecessary if \code{pixel.ysize} is provided)
#' @param ncol An integer, number of pixels along the \eqn{x}  axis (unnecessary if \code{pixel.xsize} is provided)
#'
#' @return A two column matrix containing the \eqn{x}  and \eqn{y} coordinates of the pixel centroids.
#'
#' @examples
#' ccc=coords_pix(area=square(10), nrow=10, ncol=10)
#' plot(square(10)); points(ccc)
#'
#' ccc=coords_pix(area=square(10), pixel.xsize = 2, pixel.ysize = 2)
#' plot(square(10)); points(ccc, pch=16)
#'
#' @export


coords_pix=function(area, pixel.xsize=diff(area$xrange)/ncol, pixel.ysize=diff(area$yrange)/nrow,
                nrow=diff(area$yrange)/pixel.ysize, ncol=diff(area$xrange)/pixel.xsize)
{
  if(ncol!=round(ncol,0)|nrow!=round(nrow,0))
    cat("please check the pixel size, it must be set so that nrow and ncol are integer numbers")
  xx<-seq(pixel.xsize/2,(area$xrange[2]-pixel.xsize/2),l=ncol)
  yy<-seq(pixel.ysize/2,(area$yrange[2]-pixel.ysize/2),l=nrow)

  coords.mat<-cbind(x=rep(xx,each=nrow),y=rep(rev(yy),ncol))

  return(coords.mat)
}
###########################################


###########################################
#3)adjacency matrix for a chosen distance range

#'Adjacency matrix.
#'
#'\code{adj_mat} builds an upper-triangular adjacency matrix for the set of points/areas
#'in a chosen distance range.
#'
#'The adjacency matrix is a square matrix, with each row corresponding to a point/area, and with
#'1 values along the row marking the points/areas that are considered 'neighbours' or 'pairing/coupling'.
#'In the context of spatial entropy, an adjacency matrix may take two roles.
#'If Karlstrom and Ceccato's entropy is computed, the adjacency matrix identifies
#'what areas are neighbours, i.e. what areas enter the computation of the
#'local entropy of a specific area.
#'If a spatial entropy based on the transformed variable \eqn{Z} is computed (see \code{\link{shannonZ}}
#'for details on \eqn{Z}), the adjacency matrix
#'identifies what pairs/couples of points/areas must be considered for the computation, according to
#'the chosen distance range of interest.
#'
#' @param dist.mat An upper-triangular matrix of Euclidean distances, as returned by [euclid_dist()].
#' @param dd0 Numeric, minimum distance for the neighbourhood/couples/pairs, 0 by default.
#' @param dd1 Numeric, maximum distance for the neighbourhood/couples/pairs.
#'
#' @return An \eqn{n}x\eqn{n} upper-triangular adjacency matrix
#' (where \eqn{n} in the data vector length)
#' with value 1 if two units are neighbours or form a couple/pair, 0 otherwise.
#'
#' @examples
#' dist.mat=euclid_dist(cbind(rep(1:5, each=5), rep(1:5,5)))
#' plot(cbind(rep(1:5, each=5), rep(1:5,5)))
#' adj_mat(dist.mat, dd1=dist.mat[1,2]) #for the contiguity matrix
#' adj_mat(dist.mat, 1, 3)
#'
#' @export

adj_mat=function(dist.mat, dd0=0, dd1)
{
  if(dd0==0) {adj.mat=(dist.mat<=(dd1+1e-10))+0
  } else {
    adj.mat0=(dist.mat<=(dd0+1e-10))+0
    adj.mat1=(dist.mat<=(dd1+1e-10))+0
    adj.mat=adj.mat1-adj.mat0
  }
  return(adj.mat)
}
###########################################


###########################################
#4) adjacency matrix list

#'Adjacency list for spatial entropy.
#'
#'\code{adj_list} builds a list of adjacency matrices for the computation of Altieri's entropy,
#'one for each possible distance range within the observation area.
#'
#'This function is needed for the computation of Altieri's spatial entropy.
#'After defining the distance classes in which the observation window has to be partitioned,
#'for each class ad adjacency matrix is constructed. Each adjacency matrix identifies what
#'pairs of points/areas fall within the specific distance range is the basis
#'for the computation of the local term of Altieri's spatial entropy.
#'
#' @param dist.mat An upper-triangular matrix of Euclidean distances, as returned by [euclid_dist()].
#' @param dist.breaks A numeric vector with the breaks for the partition of the maximum distance into distance classes.
#'                    It must start at 0 and end at the maximum possible distance within the observation area.
#'
#' @return A list of length `length(dist.breaks)-1`; each element is an upper-triangular
#'   adjacency matrix as returned by [adj_mat()] according to the corresponding distance range.
#'
#' @examples
#' dist.breaks=c(0,1,2,5,10*sqrt(2))
#' dist.mat=euclid_dist(coords_pix(square(10), nrow=10, ncol=10))
#' my.adj.list=adj_list(dist.mat, dist.breaks)
#'
#' @export

adj_list=function(dist.mat, dist.breaks)
{
  dd.vec=rep(dist.breaks,each=2)
  dd.vec=dd.vec[-c(1, length(dd.vec))]
  DD=matrix(dd.vec, ncol=2, byrow=TRUE)
  n.dist=nrow(DD)

  A.list=vector("list", n.dist)
  for(dd in 1:n.dist)
    A.list[[dd]]=adj_mat(dist.mat, dd0=DD[dd,1], dd1=DD[dd,2])
  return(A.list)
}

###########################################
