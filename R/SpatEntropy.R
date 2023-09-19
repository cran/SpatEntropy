#'SpatEntropy: a package for computing spatial entropy measures.
#'
#'The heterogeneity of spatial data presenting a finite number of categories
#'can be measured via computation of spatial entropy. Functions are available
#'for the computation of the main entropy and spatial entropy measures in the
#'literature. They include the traditional version of Shannon's entropy,
#'Batty's spatial entropy, O'Neill's entropy, Li and Reynolds' contagion index,
#'Karlstrom and Ceccato's entropy, Leibovici's entropy, Parresol and Edwards'
#'entropy and Altieri's entropy.
#'The package is able to work with lattice and point data. A step-by-step guide
#'for new users can be found in the first referenced article.
#'
#'References:
#'
#'ALTIERI L., COCCHI D., ROLI G. (2021). Spatial entropy for biodiversity and environmental data: The
#'R-package SpatEntropy. Environmental Modelling and Software
#'
#'ALTIERI L., COCCHI D., ROLI G. (2019). Advances in spatial entropy measures.
#'Stochastic Environmental Research and Risk Assessment
#'
#'ALTIERI L., COCCHI D., ROLI G. (2019). Measuring heterogeneity in urban expansion via spatial entropy.
#'Environmetrics, 30(2), e2548
#'
#'ALTIERI L., COCCHI D., ROLI G. (2018). A new approach to spatial entropy measures.
#'Environmental and Ecological Statistics, 25(1), 95-110
#'
#'Altieri, L., D. Cocchi, and G. Roli (2017). The use of spatial information in entropy measures.
#' arXiv:1703.06001
#'
#'Batty, M. (1974). Spatial entropy. Geographical Analysis 6, 1-31.
#'
#'Batty, M. (1976). Entropy in spatial aggregation. Geographical Analysis 8, 1-21.
#'
#'EEA (2011). Corine land cover 2000 raster data. Technical Report, downloadable at
#'http://www.eea.europa.eu/data-and-maps/ data/corine-land-cover-2000-raster-1.
#'
#'Karlstrom, A. and V. Ceccato (2002). A new information theoretical measure of global and local
#'spatial association. The Review of Regional Research 22, 13-40.
#'
#'Leibovici, D. (2009). Defining spatial entropy from multivariate distributions of co-occurrences.
#' Berlin, Springer: In K. S. Hornsby et al. (eds.): 9th International Conference on Spatial Information
#' Theory 2009, Lecture Notes in Computer Science 5756, 392-404.
#'
#'Li, H. and J. Reynolds (1993). A new contagion index to quantify spatial patterns of landscapes.
#' Landscape Ecology 8(3), 155-162.
#'
#'O'Neill, R., J. Krummel, R. Gardner, G. Sugihara, B. Jackson, D. DeAngelis, B. Milne, M. Turner,
#' B. Zygmunt, S. Christensen, V. Dale, and R. Graham (1988). Indices of landscape pattern. Landscape
#' Ecology 1(3), 153-162.
#'
#'Parresol, B. and L. Edwards (2014). An entropy-based contagion index and its sampling properties
#' for landscape analysis. Entropy 16(4), 1842-1859.
#'
#'Shannon, C. (1948). A mathematical theory of communication. Bell Dyditem Technical Journal 27,
#' 379-423, 623-656.
#'
#'@docType package
#'@name SpatEntropy
#'@import spatstat
#'@aliases SpatEntropy-package

NULL
