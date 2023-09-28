\name{Estimation of Population from Pycnophylatic Interpolation}
\alias{estimate.pycno}

\title{
Estimate populations for a set of zones based on pycnophylactic interpolation
}
\description{
Given a \code{\link{SpatialGridDataFrame}} of population estimates and a set of polygons,
compute a population  estimate based on Tobler's pycnophylactic interpolation
algorithm for each zone. The result is a vector.
}
\usage{
estimate.pycno(sgdf, spdf)
}

\arguments{
  \item{sgdf}{
A \code{\link{SpatialGridDataFrame}} containing the ouput of a pycnophylatic interpolation,
such as those produced by \code{\link{pycno}}.
}
  \item{spdf}{
A \code{\link{SpatialPolygonsDataFrame}}, giving the polygons for which estimates are wanted.
}
}

\details{
Takes the estimate of population density for each pixel,  checks which polygon each pixel is in, 
and aggregates them. Accuracy depends on the scale of pixels in the initial interpolation.
}
\value{
A vector in which each each pixel set at the estimated population aggregation to each zone in \code{spdf}.
}
\references{
Tobler, W.R. (1979) \emph{Smooth Pycnophylactic Interpolation for Geographical Regions}. Journal of the American Statistical Association, v74(367) pp. 519-530.
}
\author{
Chris Brunsdon
}
\note{
Pycnophylatic interpolation has the property that the sum of the estimated values associated with all of the pixels in any polygon equals the supplied population for that polygon.  A further property is that all pixel values are greater than or equal to zero.  The method is generally used to obtain pixel-based population estimates when total populations for a set of irregular polygons (eg. counties) are known.
}

\seealso{
\code{\link{pycno}}
}

\examples{
library(sp)
# Read in data for North Carolina as a SpatialPolygonsDataFrame
#nc.sids <- readShapeSpatial(system.file("shapes/sids.shp", package="maptools")[1], 
#  IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))

nc.sids <- as(sf::st_read(system.file("shape/nc.shp", package="sf")), "Spatial")
row.names(nc.sids) <- as.character(nc.sids$FIPSNO)

# Compute the pycnophylactic surface for 1974 births as a SpatialGridDataFrame
# Note probably shouldn't really base grid cells on Lat/Long coordinates
# This example just serves to illustrate the use of the functions
births74 <- pycno(nc.sids,nc.sids$BIR74,0.05,converge=1)

# Create a new 'blocky' set of zones
#blocks <- gUnionCascaded(nc.sids,1*(coordinates(nc.sids)[,2] > 36) + 
#  2*(coordinates(nc.sids)[,1] > -80))

crds <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(sf::st_as_sf(nc.sids)),
 of_largest_polygon = TRUE))
block_ID <- 1*(crds[,2] > 36) + 2*(crds[,1] > -80)
temp <- sf::st_as_sf(nc.sids)
temp$block_ID <- block_ID
blocks <- as(aggregate(temp, by=list(temp$block_ID), head, n=1), "Spatial")

# Plot the blocky zones
plot(blocks)
# Aggregate data to them
estimates <- estimate.pycno(births74,blocks)
# Write the estimates on to the map
text(coordinates(blocks),as.character(estimates))
}
\keyword{ smoothing }
\keyword{ spatial }
