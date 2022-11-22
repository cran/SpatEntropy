###########################################
#THIS FILE CONTAINS
#1) SHANNON'S ENTROPY OF X
#2) SHANNON'S ENTROPY OF X^2
#3) SHANNON'S ENTROPY OF Z
#4) SHANNON'S ENTROPY OF Z^2
###########################################


###########################################
#1) SHANNON'S ENTROPY

#'Shannon's entropy.
#'
#'This function computes Shannon's entropy of a variable \eqn{X} with a finite number of categories. Shannon's entropy is a non-spatial measure.
#'
#'Shannon's entropy measures the heterogeneity of a set of categorical data. It
#'is computed as \deqn{H(X)=\sum p(x_i) \log(1/p(x_i))} where \eqn{p(x_i)} is the
#'probability of occurrence of the \eqn{i}-th category, here estimated, as usual, by its relative
#'frequency. This is both the non parametric and the maximum likelihood estimator for entropy.
#'Shannon's entropy varies between 0 and \eqn{\log(I)}, \eqn{I} being the
#'number of categories of the variable under study. The relative version of Shannon's entropy, i.e. the entropy divided by
#'\eqn{\log(I)}, is also computed, under the assumption that all data categories are present in the dataset.
#'The relative entropy is useful for comparison across datasets with differen \eqn{I}.
#'The function is able to work with lattice data with missing data, as long as they are specified as NAs:
#'missing data are ignored in the computations.
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#'   Alternatively, a marked \code{ppp} object.
#'
#' @return a list of three elements:
#' \itemize{
#' \item `probabilities` - a table with the estimated probabilities (relative frequencies) for all data categories;
#' \item `shann` - Shannon's entropy;
#' \item `rel.shann` - Shannon's relative entropy.
#' }
#'
#' @examples
#' #NON SPATIAL DATA
#' shannon(sample(1:5, 50, replace=TRUE))
#'
#' #POINT DATA
#' #requires marks with a finite number of categories
#' data.pp=runifpoint(100, win=square(10))
#' marks(data.pp)=sample(c("a","b","c"), 100, replace=TRUE)
#' shannon(marks(data.pp))
#'
#' #LATTICE DATA
#' data.lat=matrix(sample(c("a","b","c"), 100, replace=TRUE), nrow=10)
#' shannon(data.lat)
#'
#' @export


shannon=function(data)
{
  if(!is.matrix(data) & !is.vector(data) & ! spatstat.geom::is.ppp(data))
    stop("Please provide the dataset as a matrix, vector, or marked ppp object")
  if(is.matrix(data) | is.vector(data)) datavec=c(data) else
    if(spatstat.geom::is.ppp(data)){
      if(is.null(spatstat.geom::marks(data)))
      stop("Please provide marks for the point pattern: at least two categories should be present to compute Shannon's entropy.") else
      datavec=spatstat.geom::marks(data)
    }

  probs=prop.table(table(datavec))
  if(length(probs)==1) warning("There is only one category, so Shannon's entropy is 0")
  sh=-sum(probs*log(probs))
  probs=as.data.frame(probs); names(probs)=c("category", "frequency")
  return(list(probabilities=probs, shann=sh,
              rel.shann=ifelse(nrow(probs)>1,sh/log(nrow(probs)),0)))
}
###########################################


###########################################
#2) VARIANCE OF SHANNON'S ENTROPY

#'Estimated variance of Shannon's entropy.
#'
#'This function estimates the variance of Shannon's entropy of a variable \eqn{X}.
#'
#'[varshannon] estimates the
#'variance of the maximum likelihood estimator of Shannon's entropy given by
#'[shannon]. The variance is \deqn{V(H(X))=H(X)_2- H(X)^2}, where \eqn{H(X)_2} is
#'a version of Shannon's entropy (see [shannon]) where
#'the information function \eqn{\log(1/p(x_i))} is squared:
#'\deqn{H(X)_2=\sum p(x_i) \log(1/p(x_i))^2}.
#'The function is able to work with lattice data with missing data, as long as they are specified as NAs:
#'missing data are ignored in the computations.
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#'   Alternatively, a marked \code{ppp} object.
#'
#' @return the estimated variance of Shannon's entropy.
#'
#' @examples
#' #NON SPATIAL DATA
#' varshannon(sample(1:5, 50, replace=TRUE))
#'
#' #POINT DATA
#' data.pp=runifpoint(100, win=square(10))
#' marks(data.pp)=sample(c("a","b","c"), 100, replace=TRUE)
#' varshannon(marks(data.pp))
#'
#' #LATTICE DATA
#' data.lat=matrix(sample(c("a","b","c"), 100, replace=TRUE), nrow=10)
#' varshannon(data.lat)
#'
#' @export

varshannon=function(data)
{
  if(!is.matrix(data) & !is.vector(data) & ! spatstat.geom::is.ppp(data))
    stop("Please provide the dataset as a matrix, vector, or marked ppp object")
  if(is.matrix(data) | is.vector(data)) datavec=c(data) else
    if(spatstat.geom::is.ppp(data)){
      if(is.null(spatstat.geom::marks(data)))
        stop("Please provide marks for the point pattern: at least two categories should be present to compute Shannon's entropy.") else
          datavec=spatstat.geom::marks(data)
    }

  sh=shannon(datavec)$shann

  probs=prop.table(table(datavec))
  if(length(probs)==1) warning("There is only one category, so Shannon's entropy and its variance are 0")
  logprob.sq=as.numeric(log(1/probs)^2)
  sh.sq=sum(probs*logprob.sq)
  return(sh.sq-sh^2)
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
#'Shannon's entropy of \eqn{Z} varies between 0 and \eqn{\log(R)}, \eqn{R=binom(n+1,2)} (where \eqn{n} is the number of observations) being the
#'number of possible pairs of categories of the variable under study.
#'The function is able to work with lattice data with missing data, as long as they are specified as NAs:
#'missing data are ignored in the computations.
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#'   Alternatively, a marked \code{ppp} object.
#'
#' @return a list of three elements:
#' \itemize{
#'   \item `probabilities` - a table with the estimated probabilities (relative frequencies) for all \eqn{Z} categories (data pairs);
#' \item `shannZ` - Shannon's entropy of \eqn{Z};
#' \item `rel.shannZ` - Shannon's relative entropy of \eqn{Z}.
#' }
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
#' @export

shannonZ=function(data)
{
  if(!is.matrix(data) & !is.vector(data) & !is.factor(data) &
     !is.character(data)& ! spatstat.geom::is.ppp(data))
    stop("Please provide the dataset as a matrix, vector, or marked ppp object")

  if(is.matrix(data)) datavec=c(data) else
  if(is.vector(data) | is.factor(data) | is.character(data)) datavec=data else
  if(spatstat.geom::is.ppp(data)){
    if(is.null(spatstat.geom::marks(data)))
      stop("Please provide marks for the point pattern: at least two categories should be present to compute Shannon's entropy.") else
      datavec=spatstat.geom::marks(data)
    }

  datavec=datavec[!is.na(datavec)]
  ns=length(datavec)
  data.tab=table(datavec)
  Xcat=names(data.tab)
  Xcat.proxy=1:length(Xcat)
  data.tab.proxy=data.tab; names(data.tab.proxy)=Xcat.proxy
  catnames=c()
  for(ii in 1:length(Xcat)) catnames=c(catnames,paste0(Xcat[ii], "-", Xcat[ii]))
  if(length(Xcat)>1) for(ii in 1:(length(Xcat)-1)) catnames=c(catnames,paste0(Xcat[ii], "-",Xcat[(ii+1): length(Xcat)]))
  catnames.proxy=c()
  for(ii in 1:length(Xcat)) catnames.proxy=c(catnames.proxy,paste0(Xcat.proxy[ii], "-",Xcat.proxy[ii]))
  if(length(Xcat)>1)for(ii in 1:(length(Xcat)-1)) catnames.proxy=c(catnames.proxy,paste0(Xcat.proxy[ii], "-",Xcat.proxy[(ii+1): length(Xcat)]))

  probabilities=data.frame("pair"=catnames,
                           "abs.freq"=numeric(length(catnames)),
                           "rel.freq"=numeric(length(catnames)))
  if (length(data.tab)==1)
    {
    warning("there is only one category, so Shannon's entropy of Z is 0")
    probabilities$abs.freq=choose(data.tab, 2)
    }else {

    probabilities$abs.freq[1:length(Xcat)]=choose(data.tab, 2)

    for (ii in (length(Xcat)+1):length(catnames))
    {
      pair.el=unlist(strsplit(catnames.proxy[ii], split="-"))
      probabilities$abs.freq[ii]=as.numeric(choose(
        data.tab[names(data.tab.proxy)==pair.el[1]]+
        data.tab[names(data.tab.proxy)==pair.el[2]],
        2)-
        (choose(data.tab[names(data.tab.proxy)==pair.el[1]],2)+choose(data.tab[names(data.tab.proxy)==pair.el[2]],2)))
    }
    }

    probabilities$rel.freq=probabilities$abs.freq/sum(probabilities$abs.freq)
    shZ=-sum(probabilities$rel.freq[probabilities$rel.freq!=0]*
               log(probabilities$rel.freq[probabilities$rel.freq!=0]))
    return(list(probabilities=probabilities, shannZ=shZ,
                rel.shannZ=ifelse(nrow(probabilities)>1,shZ/log(nrow(probabilities)),0)))
}
###########################################


###########################################
#4) VARIANCE OF SHANNON'S ENTROPY OF Z

#'Estimated variance of Shannon's entropy of \eqn{Z}.
#'
#'This function estimates the variance of Shannon's entropy of \eqn{Z}, where \eqn{Z} identifies pairs of categories of the original study variable.
#'
#'[varshannonZ] estimates the
#'variance of the maximum likelihood estimator of Shannon's entropy of \eqn{Z} given by
#'[shannonZ]. The variance is \deqn{V(H(Z))=H(Z)_2- H(Z)^2}, where
#'\deqn{H(Z)_2=\sum p(z_r)\log(1/p(z_r))^2}.
#'The function is able to work with lattice data with missing data, as long as they are specified as NAs:
#'missing data are ignored in the computations.
#'
#' @param data A data matrix or vector, can be numeric, factor, character, ...
#'   Alternatively, a marked \code{ppp} object.
#'
#' @return the estimated variance of Shannon's entropy of \eqn{Z}.
#'
#' @examples
#' #NON SPATIAL DATA
#' data=sample(1:5, 50, replace=TRUE)
#' varshannonZ(data)
#'
#' #POINT DATA
#' data.pp=runifpoint(100, win=square(10))
#' marks(data.pp)=sample(c("a","b","c"), 100, replace=TRUE)
#' varshannonZ(marks(data.pp))
#'
#' #LATTICE DATA
#' data.lat=matrix(sample(c("a","b","c"), 100, replace=TRUE), nrow=10)
#' varshannonZ(data.lat)
#'
#' @export

varshannonZ=function(data)
{
  if(!is.matrix(data) & !is.vector(data) & ! spatstat.geom::is.ppp(data))
    stop("Please provide the dataset as a matrix, vector, or marked ppp object")
  if(is.matrix(data) | is.vector(data)) datavec=c(data) else
    if(spatstat.geom::is.ppp(data)){
      if(is.null(spatstat.geom::marks(data)))
        stop("Please provide marks for the point pattern: at least two categories should be present to compute Shannon's entropy.") else
          datavec=spatstat.geom::marks(data)
    }

  datavec=datavec[!is.na(datavec)]
  shZ=shannonZ(datavec)$shannZ
  ns=length(datavec)
  data.tab=table(datavec)
  Xcat=names(data.tab)
  catnames=c()
  for(ii in 1:length(Xcat)) catnames=c(catnames,paste0(Xcat[ii], "-", Xcat[ii]))
  for(ii in 1:(length(Xcat)-1)) catnames=c(catnames,paste0(Xcat[ii], "-", Xcat[(ii+1): length(Xcat)]))

  probabilities=data.frame("pair"=catnames,
                           "abs.freq"=numeric(length(catnames)),
                           "rel.freq"=numeric(length(catnames)))
  if (length(data.tab)==1)
  {
    warning("there is only one category, so Shannon's entropy of Z and its variance are both 0")
    probabilities$abs.freq=choose(data.tab, 2)
  }else {

    probabilities$abs.freq[1:length(Xcat)]=choose(data.tab, 2)

    for (ii in (length(Xcat)+1):length(catnames))
    {
      pair.el=unlist(strsplit(catnames[ii], split="-"))
      probabilities$abs.freq[ii]=as.numeric(choose(data.tab[names(data.tab)==pair.el[1]]+data.tab[names(data.tab)==pair.el[2]],2))-
        (choose(data.tab[names(data.tab)==pair.el[1]],2)+choose(data.tab[names(data.tab)==pair.el[2]],2))
    }
    probabilities$rel.freq=probabilities$abs.freq/sum(probabilities$abs.freq)
    logprob.sq=as.numeric(log(1/probabilities$rel.freq)^2)
    shZ.sq=sum(probabilities$rel.freq*logprob.sq)
    return(shZ.sq-shZ^2)
  }
}
###########################################


