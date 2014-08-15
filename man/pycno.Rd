\name{Pycnophylactic Interpolation}
\alias{pycno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Perform Tobler's pycnophylactic algorithm.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a \code{\link{SpatialPolygonsDataFrame}} and a set of populations for each polygon,
compute a population density estimate based on Tobler's pycnophylactic interpolation
algorithm. The result is a \code{\link{SpatialGridDataFrame}}.
}
\usage{
pycno(x, pops, celldim, r = 0.2, converge = 3, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A \code{\link{SpatialPolygonsDataFrame}} representing a set of zones for which total populations are known.
}
  \item{pops}{
%%     ~~Describe \code{pops} here~~
A list of populations,  one for each zone represented in \code{x}.
}
  \item{celldim}{
%%     ~~Describe \code{celldim} here~~
The size of the cell sides for the output SpatialGridDataFrame,  
or a SpatialGrid or SpatialGridDataFrame with the grid topology of the desired output. 
}
  \item{r}{
A relaxation parameter for the iterative step in the pycnophylactic algorithm.  Prevents
over-compensation in the smoothing step.  In practice the default value works well.
%%     ~~Describe \code{r} here~~
}
  \item{converge}{
A convergence parameter,  informing the decision on when iterative improvements on the smooth surface have converged sufficiently - see details.
%%     ~~Describe \code{converge} here~~
}
  \item{verbose}{
If set to \code{TRUE} (the default value) makes the function report the maximum change in any pixel value for each iterative step.  If set to \code{FALSE} the function is silent.
%%     ~~Describe \code{converge} here~~
}
}

\details{
This method uses an iterative approach,  and for each iteration notes the maximum change in a pixel.  When this value falls below a certain level (actually 10^(-converge) times the largest initial pixel value) the iterations stop and the current surface is reported.  
}
\value{
A \code{\link{SpatialGridDataFrame}} with each pixel set at the estimated density computed with the Pycnophylactic algorithm,  these values being stored in the data item \code{dens}. The extent of this object is set at by the bounding box of the input \code{\link{SpatialPolygonsDataFrame}}.  Pixels outside of all of the polygons are set to NA.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# Read in data for North Carolina as a SpatialPolygonsDataFrame
nc.sids <- readShapeSpatial(system.file("shapes/sids.shp", package="maptools")[1], 
   IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))

# Compute the pycnophylactic surface for 1974 births as a SpatialGridDataFrame
# Note probably shouldn't really base grid cells on Lat/Long coordinates
# This example just serves to illustrate the use of the function
# It is suggested to use a hihger value for 'converge' - this value justs speeds 
# things up for the example.
births74 <- pycno(nc.sids,nc.sids$BIR74,0.05,converge=1)

# Draw it
image(births74)

# Overlay North Carolina county boundaries for reference
plot(nc.sids,add=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smoothing }
\keyword{ spatial }% __ONLY ONE__ keyword per line
