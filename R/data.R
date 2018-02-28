###########################################
#THIS FILE CONTAINS
#1) data description for rainforest tree data (point pattern)
#2) data description for urban data (lattice)
###########################################


###########################################
#1) tree

#'Rainforest tree data.
#'
#'A marked point pattern dataset about four rainforest tree species.
#'
#'This dataset documents the presence of tree species over Barro Colorado Island, Panama.
#'Barro Colorado Island has been the focus of intensive research on lowland tropical
#'rainforest since 1923 (http://www.ctfs.si.edu). Research identified several tree species
#'over a rectangular observation window of size 1000x500 metres; the tree species
#'constitute the point data categorical mark. This dataset presents 4 species with
#'different spatial configurations: Acalypha diversifolia, Chamguava schippii,
#'Inga pezizifera and Rinorea sylvatica. The overall dataset has a total number of 7251 points.
#'The dataset is analyzed with spatial entropy measures in Altieri et al (2017) (references can be
#'found at \code{SpatEntropy}).
#'
#'@format A \code{ppp} object (see package \code{spatstat}) with 7251 points, containing:
#'\describe{
#'   \item{window}{An object of type \code{owin} (see package \code{spatstat}), the 1000x500 metres observation area}
#'   \item{x, y}{Numeric vectors with points' coordinates}
#'   \item{marks}{A character vector matching the tree species to the data points}
#' }
#' @source http://www.ctfs.si.edu

"data_rainforest"
###########################################


###########################################
#2) bologna

#'Bologna data.
#'
#'A lattice dataset with Bologna Urban Morphological Zones.
#'
#'This raster/pixel/lattice dataset comes from the EU CORINE Land Cover project (EEA, 2011) and is dated 2011.
#'It is the result of classifying the original
#'land cover data into urbanised and non-urbanised zones, known as 'Urban
#'Morphological Zones' (UMZ, see EEA, 2011). UMZ data are useful to identify shapes
#'and patterns of urban areas, and thus to detect what is known as urban sprawl.
#'The Bologna metropolitan area is extracted from the European dataset and is
#'composed by the municipality of Bologna and the surrounding municipalities.
#'The dataset is made of 120x96 pixels of size 250x250 metres.
#'
#'@format A `matrix` with 120 rows and 96 columns. Values are either 0 (non-urban) or 1 (urban).
#' @source EEA (2011). Corine land cover 2000 raster data. Technical Report, downloadable at
#' http://www.eea.europa.eu/data-and-maps/ data/corine-land-cover-2000-raster-1

"data_bologna"
###########################################
