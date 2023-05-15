\name{pycno-package}
\alias{pycno-package}
\docType{package}
\title{
Pycnophylactic Interpolation
}
\description{
Given a SpatialPolygonsDataFrame and a set of populations for each polygon,
compute a population density estimate based on Tobler's pycnophylactic interpolation
algorithm. The result is a SpatialGridDataFrame.
}
\details{
\tabular{ll}{
Package: \tab pycno\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-04-15\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
For use in conjunction with the sp package,  computes pycnophylactic surfaces given a
SpatialPolygonsDataFrame and a population for each poplygon.  A pycnophylactic surface 
is smooth,  but populations allocated to each pixel sum up to the initial polygon counts, 
when summed over the polygons contained in each pixel.
}
\author{
Chris Brunsdon

Maintainer: Francesca Bitonti \email{francesca.bitonti@phd.unict.it}
}
\references{
Tobler, W.R. (1979) \emph{Smooth Pycnophylactic Interpolation for Geographical Regions}. Journal of the American Statistical Association, v74(367) pp. 519-530.
}

\keyword{ package }

\examples{
# Read in data for North Carolina as a SpatialPolygonsDataFrame
nc.sids <- readShapeSpatial(system.file("shapes/sids.shp", package="maptools")[1], 
  IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))

# Compute the pycnophylactic surface for 1974 births as a SpatialGridDataFrame
# Note probably shouldn't really base grid cells on Lat/Long coordinates
# This example just serves to illustrate the use of the function
births74 <- pycno(nc.sids,nc.sids$BIR74,0.05, converge=1)

# Draw it
image(births74)

# Overlay North Carolina county boundaries for reference
plot(nc.sids,add=TRUE)
}
