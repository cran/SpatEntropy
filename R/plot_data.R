###########################################
#THIS FILE CONTAINS
#1) function for plotting lattice data
###########################################

###########################################
#1)

#'Plot lattice data.
#'
#'`plot_lattice` produces a gray scale plot of a matrix of
#'categorical data.
#'
#'This function allows to easily produce a gray scale map given a matrix of categorical data
#'and, optionally, the observation area. It ensures that data is displayed following the matrix order
#'(where position (1,1) corresponds to the top-left corner of plot), avoiding risks of row inversion or
#'transposition. A few #'options may be tuned: the extent of the gray scale, the title and the legend.
#'
#' @param data A matrix, can be numeric, factor, character.
#' @param win An \code{owin} object, the observation area (see package \code{spatstat}). Automatically
#'   created, if not provided, by setting the pixel size as 1.
#' @param gray.ext A vector of length two with the two extremes of the gray scale (between 0 and 1).
#' @param main Optional, a character string with the plot main title.
#' @param ribbon Logical, whether to display a ribbon showing the colour map.
#'
#' @return A gray scale plot of the categorical lattice dataset.
#'
#' @examples
#' data.lat=matrix(sample(c("a","b","c"), 100, replace=TRUE), nrow=10)
#' plot_lattice(data.lat)
#'
#' plot_lattice(data.lat, win=square(100))
#'
#' plot_lattice(data.lat, win=square(10), gray.ext=c(1,.4), ribbon=FALSE)
#'
#' @export

plot_lattice=function(data, win=spatstat::owin(xrange=c(1, ncol(data)),
                                               yrange=c(1, nrow(data))),
                      gray.ext=c(1,0), main="", ribbon=TRUE){
  color=grDevices::gray(seq(gray.ext[1],gray.ext[2],
                            length.out=length(unique(c(data)))))
  data2=data
  for(i in 1:nrow(data)) data2[i,]=data[nrow(data)-(i-1),]
  spatstat::plot.im(spatstat::as.im(data2, W=win), col=color, ribbon=ribbon, main=main)
}
###########################################

