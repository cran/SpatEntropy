###########################################
#THIS FILE CONTAINS
#1) function for building Altieri et al spatial entropy measure
###########################################


###########################################
#'Altieri's spatial entropy.
#'
#'This function computes spatial mutual information and spatial residual entropy as in Altieri et al (2017).
#'References can be found at \code{SpatEntropy}.
#'
#'The computation of Altieri's entropy starts from a point or areal dataset, for which
#'Shannon's entropy of the transformed variable \eqn{Z} (for details see \code{\link{shannonZ}})
#'\deqn{H(Z)=\sum p(z_r)\log(1/p(z_r))}
#'is computed using all
#'possible pairs within the observation area. Then, its two components
#'spatial mutual information
#'\deqn{MI(Z,W)=\sum p(w_k) \sum p(z_r|w_k)\log(p(z_r|w_k)/p(z_r))}
#' and spatial residual entropy
#'\deqn{H(Z)_W=\sum p(w_k) \sum p(z_r|w_k)\log(1/p(z_r|w_k))}
#'are calculated
#'in order to account for the overall role of space in determining
#'the data heterogeneity. Besides, starting from a partition into distance
#'classes, a list of adjacency matrices is
#'built using [adj_mat()], which identifies what pairs of units
#'must be considered for each class. Spatial mutual information and spatial residual
#'entropy are split into local terms according to the chosen distance breaks, so that the role of space
#'can be investigated.
#'
#'@param data A data matrix or vector, can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is the mark vector.
#' @param adj.list A list of adjacency matrices, provided by user or returned by [adj_list()].
#'   Each element of the list contains a binary adjacency matrix for building pairs at a specific distance,
#'   provided by user or returned by [adj_mat()].
#' @param shannZ The output of [shannonZ()]: Shannon's entropy of \eqn{Z}
#'            as well as the pairs frequency table.
#' @param missing.cat Optional, a vector with the names of all categories that are absent in `data`.
#'
#' @return A list with elements:
#'\itemize{
#'   \item `mut.global` the spatial mutual information
#'   \item `res.global` the global residual entropy
#'   \item `shannZ` Shannon's entropy of \eqn{Z}
#'   \item `mut.local` the partial information terms
#'   \item `res.local` the partial residual entropies
#'   \item `pwk` the spatial weights for each distance range
#'   \item `pzr.marg` the marginal probability distribution of \eqn{Z}
#'   \item `pzr.cond` a list with the conditional probability distribution of \eqn{Z} for each distance range
#'   \item `Q` the number of pairs (realizations of \eqn{Z})
#'   \item `Qk` the number of pairs for each distance range.
#'}
#'
#' @examples
#' data=matrix(sample(1:5, 25, replace=TRUE), nrow=5)
#' plot(as.im(data, W=square(5)))
#' dist.breaks=c(0,1,2,5,5*sqrt(2))
#' dist.mat=euclid_dist(coords_pix(square(5), nrow=5, ncol=5))
#' my.adj.list=adj_list(dist.mat, dist.breaks) #see adj_list
#' my.shZ=shannonZ(data)
#' spat_entropy(data=data, adj.list=my.adj.list, shannZ=my.shZ)
#'
#' @export

spat_entropy=function(data, adj.list, shannZ, missing.cat=NULL)
{
  ##ingredients:
  #1) Z marginal frequencies
  P.zr=shannZ$probabilities$proportion
  names(P.zr)=shannZ$probabilities$couple

  #2) W marginal frequencies and Z|wk conditional frequencies
  n.dist=length(adj.list)
  QQ=sum(shannZ$probabilities$abs.frequency)
  P.zr.cond.wk=vector("list", n.dist)
  P.wk=numeric(n.dist)
  for (dd in 1:n.dist)
  {
    output=pair_count(data, adj.list[[dd]], missing.cat)
    P.zr.cond.wk[[dd]]=output$probabilities$proportion
    names(P.zr.cond.wk[[dd]])=output$probabilities$couple
    P.wk[dd]=output$Qk/QQ
  }

  ##partial terms
  res.local=mut.local=numeric(n.dist)
  for(dd in 1:n.dist)
  {
    cond.probs=as.numeric(P.zr.cond.wk[[dd]][P.zr.cond.wk[[dd]]>0])
    marg.probs=as.numeric(P.zr[P.zr.cond.wk[[dd]]>0])
    res.local[dd]=sum(cond.probs*log(1/cond.probs))
    mut.local[dd]=sum(cond.probs*log(cond.probs/marg.probs))
  }
  res.global=sum(P.wk*res.local)
  mut.global=sum(P.wk*mut.local)

  ##output
  return(list(mut.global=mut.global, res.global=res.global, shannZ=shannZ$shannon.Z,
              mut.local=mut.local, res.local=res.local,
              pwk=P.wk, pzr.marg=P.zr, pzr.cond=P.zr.cond.wk,
              Q=QQ, Qk=P.wk*QQ))
}
###########################################

