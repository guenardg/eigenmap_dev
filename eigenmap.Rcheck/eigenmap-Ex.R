pkgname <- "eigenmap"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "eigenmap-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('eigenmap')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Euclid")
### * Euclid

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Euclid
### Title: Calculation of the Euclidean Distance
### Aliases: Euclid

### ** Examples

### A set of reference points:
x <- cbind(c(1,4,5,2,8,4), c(3,6,7,1,3,2))
dimnames(x) <- list(LETTERS[1:6], c("x", "y"))

## The pairwise Euclidean distances among the reference points: 
d1 <- Euclid(x)
d1

## That result is the same as that obtained from function dist:
d2 <- dist(x, method = "euclidean")
all(d1 == d2)

## A second set of points:
y <- cbind(c(3,5,7), c(3,6,8))
dimnames(y) <- list(LETTERS[7:9], c("x", "y"))

## The distances between the points in y (rows) and x (columns):
Euclid(x, y)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Euclid", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("eigenmap")
### * eigenmap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: eigenmap
### Title: Spatial Eigenvector Maps
### Aliases: eigenmap eigenmap.score

### ** Examples

### A unevenly sampled surface.

data(mite)

## Example using the principal coordinates of the square root of the
## (Euclidean) distances:
map <- eigenmap(x = as.matrix(mite.geo), weighting = wf.sqrd)
map
## plot(map)

## Example using the radial basis functions (RBF):
map <- eigenmap(x = as.matrix(mite.geo), weighting = wf.RBF)
map
## plot(map)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("eigenmap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("geodesics")
### * geodesics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: geodesics
### Title: Calculation of Geodesic Distances
### Aliases: geodesics

### ** Examples

### First example: locations spread throughout the world

coords <- cbind(c(43,22,9,12,-40,72,-86,-22),
                c(-135,22,0,1,-45,12,27,-139))

res_hav <- geodesics(coords)  ## Default: the haversine formula
res_hav

res_vif <- geodesics(coords, method = "Vincenty")
res_vif

attr(res_vif,"niter") ## The numbers of iterations
res_vif-res_hav       ## Absolute difference
200*(res_vif-res_hav)/(res_vif+res_hav) ## Large relative difference

### Second example: locations nearer from one another

coords <- cbind(c(45.01,44.82,45.23,44.74),
                c(72.03,72.34,71.89,72.45))

res_hav <- geodesics(coords)
res_vif <- geodesics(coords, method = "Vincenty")
res_vif-res_hav       ## Absolute difference
200*(res_vif-res_hav)/(res_vif+res_hav) ## Relative difference are smaller




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("geodesics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mite")
### * mite

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mite
### Title: The Oribatid Mite Data Set
### Aliases: mite mite.species mite.env mite.geo
### Keywords: mite

### ** Examples

data(mite)
summary(mite.species)
summary(mite.env)
summary(mite.geo)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mite", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("weighting-functions")
### * weighting-functions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: weighting-functions
### Title: Weighting Functions for Spatial Eigenvector Map
### Aliases: weighting-functions wf.sqrd wf.RBF wf.PCNM wf.binary wf.Drayf1
###   wf.Drayf2 wf.Drayf3

### ** Examples

locations <- c(1,2,4,7,10,14,17,21)
D <- dist(locations)
wf.sqrd(D)
wf.RBF(D, wpar = 0.1)
wf.binary(D, c(0,5))
wf.PCNM(D, c(0,5))
wf.Drayf1(D, c(0,5))
wf.Drayf2(D, c(0,5), 0.5)
wf.Drayf3(D, c(0,5), 0.5)

emap <- eigenmap(D, locations, wf.Drayf2, c(0,5), 0.5)
emap

emap <- eigenmap(D, locations, wf.Drayf3, c(0,5), 0.25)
emap

emap <- eigenmap(D, locations, wf.RBF, wpar = 0.1)
emap




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("weighting-functions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
