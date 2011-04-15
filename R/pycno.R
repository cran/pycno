#                          Pycnophylactic interpolation
# =============================================================================
#
# The key function is called pycno 
# it takes the form
#    result <- pycno(shapefile,popns,cellsize)
# Where shapefile is a polygon shape object
#       popns is a list of populations for each polygon in the shapefile
#       cellsize is the size of the grid cell for the output
#       result is a grid object with the pycnophylactic interpolated surface
#    
# You can visualize the result with something like
#
#       image(result)
#
# ==============================================================================

#require(stats)
#require(GISTools)


  
.pycno.core <- function(zones,pops,r=0.2,converge=3,verbose=TRUE) {
 smooth2D <- function(x) {
     mval <- mean(x)
     s1d <- function(s) unclass(filter(s,c(0.5,0,0.5)))
     pad <- rbind(mval,cbind(mval,x,mval),mval)
     pad <- (t(apply(pad,1,s1d)) + apply(pad,2,s1d))/2
     return(pad[2:(nrow(x)+1),2:(ncol(x)+1)])}
   
   
 correct2Dm <- function(x,zones,pops) {
     zone.list <- sort(unique(array(zones)))
     for (item in zone.list) {
       zone.set <- (zones == item)
       correct <- pops[item]/sum(x[zone.set])
       x[zone.set] <- x[zone.set]*correct }
     return(x)}
     
 correct2Da <- function(x,zones,pops) {
     zone.list <- sort(unique(array(zones)))
     for (item in zone.list) {
       zone.set <- (zones == item)
       correct <- (pops[item] - sum(x[zone.set]))/sum(zone.set)
       x[zone.set] <- x[zone.set] + correct }
     return(x)}
   
 populate <- function(zones,pops) {
     x <- zones*0
     zone.list <- sort(unique(array(zones)))
     for (item in zone.list) {
       zone.set <- (zones == item)
       x[zone.set] <- pops[item]/sum(zone.set)}
     return(x)}
   
 x <- populate(zones,pops)
 stopper <- max(x)
 stopper <- stopper*10^(-converge)
 repeat {
   old.x <- x
   sm <- smooth2D(x)
   x <- x*r + (1-r)*sm
   x <- correct2Da(x,zones,pops)
   x[x<0] <- 0
   x <- correct2Dm(x,zones,pops)
   if (verbose) {
     cat(sprintf("Maximum Change: %12.5f - will stop at %12.5f\n", max(abs(old.x - x)),stopper))}
   if (max(abs(old.x - x)) < stopper) break }
 return(x)} 
 
.poly2grid <- function(x,celldim) {
  bbx <- slot(x,'bbox')
  offset <- bbx[,1]
  extent <- bbx[,2] - offset
  shape <- ceiling(extent / celldim)
  sg <- SpatialGrid(GridTopology(offset,c(celldim,celldim),shape))
  res <- overlay(SpatialPoints(coordinates(sg)),x)
  spdf <- SpatialPixelsDataFrame(coordinates(sg),data.frame(zone=res))
  return(as(spdf,"SpatialGridDataFrame"))}

.grid.matrix <- function(gr,index=1) {
  gr.dim <- slot(getGridTopology(gr),"cells.dim")
  res <- gr[[index]]
  dim(res) <- gr.dim
  attr(res,'na') <- is.na(res)
  res[is.na(res)] <- max(res,na.rm=T) + 1
  return(res)}

 
.matrix2grid <- function(gr,mat) {
  if (!is.null(attr(mat,'na'))) mat[attr(mat,'na')] <- NA
  spdf <- SpatialPixelsDataFrame(coordinates(gr),data.frame(dens=array(mat)))
  return(as(spdf,"SpatialGridDataFrame"))}

pycno <- function(x,pops,celldim,r=0.2,converge=3,verbose=TRUE) {
  gr <- .poly2grid(x,celldim)
  gm <- .grid.matrix(gr)
  pops2 <- c(pops,0)
  pm <- .pycno.core(gm,pops2,r=r,converge=converge,verbose=verbose)
  result <- .matrix2grid(gr,pm)
  return(result) }
  
