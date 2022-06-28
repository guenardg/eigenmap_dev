## **************************************************************************
##
##    (c) 2018-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Package eigenmap description**
##
##    This file is part of eigenmap
##
##    eigenmap is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    eigenmap is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with eigenmap. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' \packageTitle{eigenmap}
#' 
#' @description \packageDescription{eigenmap}
#' The eigenfunctions are obtained in three steps: 1) a distance matrix is
#' calculated from the locations of samples in space (or the sampling
#' organisation through time). 2) From that distance matrix, a matrix of Moran
#' spatial weights is obtained; this is the same matrix as used to calculate
#' Moran's autocorrelation index, hence the name. And 3) the spatial weight
#' matrix is eigenvalue-decomposed after centring the rows and columns of the
#' spatial weight matrix.
#' 
#' @docType package
#' 
#' @name eigenmap-package
#' 
#' @details Function \code{\link{eigenmap}} calculates spatial eigenvector maps
#' following the approach outlined in Dray et al. (2006). It returns a
#' \link{eigenmap-class} object. The package also features methods to print
#' (\code{\link{print.eigenmap}}) and plot (\code{\link{plot.eigenmap}}) these
#' objects. Function \code{\link{eigenmap.score}} can be used to make
#' predictions for spatial models built from the eigenfunctions of
#' \code{\link{eigenmap}} using distances between one or more target locations
#' and the sampled locations for which the spatial eigenvector map was built.
#' 
## The package also features an exemplary data set...
#' 
#' The DESCRIPTION file:
#' \packageDESCRIPTION{eigenmap}
#' \packageIndices{eigenmap}
#' 
#' @author \packageAuthor{eigenmap}
#' Maintainer: \packageMaintainer{eigenmap}
#' 
#' @references
#' Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbor
#' matrices (PCNM). Ecol. Modelling 196: 483-493
#' 
#' @seealso
#' Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English edition.
#' Elsevier Science B.V., Amsterdam, The Neatherlands.
#' 
NULL
#' 
