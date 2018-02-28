###########################################
#THIS FILE CONTAINS
#1) SHANNON'S ENTROPY OF X
#2) SHANNON'S ENTROPY OF X^2
#3) SHANNON'S ENTROPY OF Z
#4) SHANNON'S ENTROPY OF Z^2
###########################################


###########################################
#1) SHANNON'S ENTROPY OF X

#'Shannon's entropy.
#'
#'This function computes Shannon's entropy of a variable \eqn{X} with a finite number of categories.
#'
#'Shannon's entropy measures the heterogeneity of a set of categorical data. It
#'is computed as \deqn{H(X)=\sum p(x_i) \log(1/p(x_i))} where \eqn{p(x_i)} is the
#'probability of occurrence of the \eqn{i}-th category, here estimated by its relative
#'frequency. This is both the non parametric and the maximum likelihood estimator.
#'Shannon's entropy varies between 0 and \eqn{\log(I)}, \eqn{I} being the
#'number of categories of the variable under study.
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is the mark vector.
#'
#' @return Estimated probabilities for all data categories, and Shannon's entropy.
#'
#' @examples
#' #NON SPATIAL DATA
#' shannonX(sample(1:5, 50, replace=TRUE))
#'
#' #POINT DATA
#' data.pp=runifpoint(100, win=square(10))
#' marks(data.pp)=sample(c("a","b","c"), 100, replace=TRUE)
#' shannonX(marks(data.pp))
#'
#' #LATTICE DATA
#' data.lat=matrix(sample(c("a","b","c"), 100, replace=TRUE), nrow=10)
#' shannonX(data.lat)
#'
#' @export


shannonX=function(data)
{
  if(!is.matrix(data)&!is.vector(data)) print("data must be a matrix or a vector")
  datavec=c(data)
  probs=table(datavec)/sum(table(datavec))
  sh=-sum(probs*log(probs))
  probs=as.data.frame(probs); names(probs)=c("category", "frequency")
  return(list(probabilities=probs, shannon=sh))
}
###########################################


###########################################
#2) SHANNON'S ENTROPY OF X^2

#'Shannon's entropy with a squared information function.
#'
#'This function computes Shannon's entropy of \eqn{X} with the square of the information function.
#'
#'This computes a version of Shannon's entropy (see [shannonX()]) where
#'the information function \eqn{\log(1/p(x_i))} is squared:
#'\deqn{H(X)_2=\sum p(x_i) \log(1/p(x_i))^2}
#'It is useful for estimating the
#'variance of the maximum likelihood estimator of Shannon's entropy given by
#'[shannonX()].
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is the mark vector.
#'
#' @return Estimated probabilities for all data categories, and Shannon's entropy of \eqn{X} with a squared
#' information function.
#'
#' @examples
#' #NON SPATIAL DATA
#' shannonX_sq(sample(1:5, 50, replace=TRUE))
#'
#' #POINT DATA
#' data.pp=runifpoint(100, win=square(10))
#' marks(data.pp)=sample(c("a","b","c"), 100, replace=TRUE)
#' shannonX_sq(marks(data.pp))
#'
#' #LATTICE DATA
#' data.lat=matrix(sample(c("a","b","c"), 100, replace=TRUE), nrow=10)
#' shannonX_sq(data.lat)
#'
#' @export

shannonX_sq=function(data)
{
  if(!is.matrix(data)&!is.vector(data)) print("data must be a matrix or a vector")
  datavec=c(data)
  probs=table(datavec)/sum(table(datavec))
  logprob.sq=as.numeric(log(1/probs)^2)
  sh=sum(probs*logprob.sq)
  probs=as.data.frame(probs)
  probs=data.frame(probs, logprob.sq); names(probs)=c("category", "frequency", "logprob.sq")
  return(list(probabilities=probs, shannon.square=sh))
}
###########################################


###########################################
#3) SHANNON'S ENTROPY OF Z

#'Shannon's entropy of the transformed variable \eqn{Z}.
#'
#'This function computes Shannon's entropy of variable \eqn{Z},
#'where \eqn{Z} identifies pairs of realizations of the variable of interest.
#'
#'Many spatial entropy indices are based on the trasformation \eqn{Z} of the study variable,
#'i.e. on pairs (unordered couples) of realizations of the variable of interest. 'Unordered couples'
#'means that the relative spatial location is irrelevant, i.e. that a couple
#'where category \eqn{i} occurs at the left of category \eqn{j} is identical to a couple
#'where category \eqn{j} occurs at the left of category \eqn{i}.
#'When all possible pairs occurring within the observation areas are considered,
#'Shannon's entropy of the variable \eqn{Z} may be computed as
#'\deqn{H(Z)=\sum p(z_r)\log(1/p(z_r))}
#'where \eqn{p(z_r)} is the probability of the \eqn{r}-th pair of realizations, here
#'estimated by its relative frequency.
#'Shannon's entropy of \eqn{Z} varies between 0 and \eqn{\log(R)}, \eqn{R} being the
#'number of possible pairs of categories of the variable under study.
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#'   If the dataset is a point pattern, `data` is the mark vector.
#' @param missing.cat Optional, a vector with the names of all categories that are absent in `data`.
#'
#' @return Estimated probabilities for all \eqn{Z} categories (data pairs), and Shannon's entropy of \eqn{Z}.
#'
#' @examples
#' #NON SPATIAL DATA
#' shannonZ(sample(1:5, 50, replace=TRUE))
#'
#' #POINT DATA
#' data.pp=runifpoint(100, win=square(10))
#' marks(data.pp)=sample(c("a","b","c"), 100, replace=TRUE)
#' shannonZ(marks(data.pp))
#'
#' #LATTICE DATA
#' data.lat=matrix(sample(c("a","b","c"), 100, replace=TRUE), nrow=10)
#' shannonZ(data.lat)
#'
#' #when categories are missing
#' shannonZ(data.lat, missing.cat=c("d", "e"))
#'
#' @export

shannonZ=function(data, missing.cat=NULL) #ncat= num categorie di X
{
  if(!is.matrix(data)&!is.vector(data)) print("data must be a matrix or a vector")
  datavec=c(data)
  print("building Z variable...")
  adj.mat=matrix(0,length(datavec), length(datavec))
  for(j in 1:(length(datavec)-1)) adj.mat[j, (j+1):length(datavec)]=1
  output=pair_count(datavec, adj.mat, missing.cat)

  print("computing entropy...")
  prop=output$probabilities$proportion[output$probabilities$proportion>0]
  localH=sum(prop*log(1/prop))

  print("done")
  return(list(probabilities=output$probabilities, shannon.Z=localH))
}
###########################################


###########################################
#4) SHANNON'S ENTROPY OF Z^2

#'Shannon's entropy of \eqn{Z} with a squared information function.
#'
#'This function computes Shannon's entropy of \eqn{Z} with the square of the information function.
#'
#'This computes a version of Shannon's entropy of \eqn{Z} (see [shannonZ()]) where
#'the information function \eqn{\log(1/p(z_r))} is squared:
#'\deqn{H(Z)_2=\sum p(z_r)\log(1/p(z_r))^2}
#'It is useful for estimating the
#'variance of the maximum likelihood estimator of Shannon's entropy given by
#'[shannonZ()].
#'
#' @param shannZ Output of [shannonZ()]
#'
#' @return Estimated probabilities for all \eqn{Z} categories (data pairs), and Shannon's entropy of \eqn{Z} with a squared information function.
#'
#' @examples
#' #NON SPATIAL DATA
#' shZ=shannonZ(sample(1:5, 50, replace=TRUE))
#' shannonZ_sq(shZ)
#'
#' #POINT DATA
#' data.pp=runifpoint(100, win=square(10))
#' marks(data.pp)=sample(c("a","b","c"), 100, replace=TRUE)
#' shZ=shannonZ(marks(data.pp))
#' shannonZ_sq(shZ)
#'
#' #LATTICE DATA
#' data.lat=matrix(sample(c("a","b","c"), 100, replace=TRUE), nrow=10)
#' shZ=shannonZ(data.lat)
#' shannonZ_sq(shZ)
#'
#' @export

shannonZ_sq=function(shannZ)
{
  probs=shannZ$probabilities$proportion[shannZ$probabilities$proportion>0]
  log.square=log(1/probs)^2
  sh=sum(probs*log.square)
  logprob.sq=numeric(nrow(shannZ$probabilities))
    logprob.sq[shannZ$probabilities$proportion>0]=log.square
  probabilities=cbind(shannZ$probabilities, logprob.sq)
  return(list(probabilities=probabilities, shannonZ.square=sh))
}
###########################################


