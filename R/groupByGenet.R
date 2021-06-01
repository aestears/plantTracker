#' Assign genetIDs to grouped ramet polygons
#'
#' This function assigns a unique 'genetID' to polygons if they are within a
#' user-defined distance from one another.
#'
#' If polygons are 'grouped,' they are given the same number in the 'genetID'
#' column. This assignment is made using network analysis to group together
#' polygons that are closest to one another. In the context of
#' \pkg{PlantTracker}, this function was designed to group ramets together into
#' one genet, which is a genetically distinct individual. This function was
#' designed for use within the \code{\link{assign}} function (and then the
#' \code{\link{trackSpp}} function), but can be used independently.
#'
#' @param sf An sf object that contains polygons to be grouped. Typically should
#' include data only for one species, one quadrat, and one year. For intended
#' use, this dataset should be of the format described in
#' \code{\link{grasslandData}}. This function will run if 'sf' contains
#' only 'geometry' data, but it is *strongly* recommended that other columns
#' specified in \code{\link{grasslandData}} are included.
#' @param buffGenet A numeric argument indicating the maximum distance
#' between individual ramets (polygons) that the function will group together as
#' one genet. Note that this argument is in the same units that are used in the
#' 'sf' argument. For example, if buffGenet = 0.005 and we use the
#' \code{\link{grasslandData}} (in which measurements are in meters), then
#' polygons that are 0.005 m (0.5 cm) apart or less will be grouped as a genet.
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return A numeric vector of unique genetIDs that is as long as the number of
#' rows in 'sf.' Each \code{i} in the output is the genetID for the \code{i}th
#' row in 'sf'.
#'
#' @examples
#'
<<<<<<< HEAD
#' @importFrom igraph clusters graph_from_adjacency_matrix
#' @import Matrix
#' @import sf
=======
#' dat <- grasslandData[grasslandData$Site=="CO" &
#' grasslandData$Quad == "ungz_5a" &
#' grasslandData$Species == "Bouteloua gracilis" &
#' grasslandData$Year == 1997,]
#'
#' groupByGenet(sf = dat, buffGenet = 0.005)
>>>>>>> f42d0628dc7bfb9579e6c64712754982addb65cb
#'
#' @seealso [assign()] and [trackSpp()], \pkg{PlantTracker} functions that apply
#' this function across multiple species, quadrats, and years.
#'
#' @import sf
#' @import Matrix
#' @importFrom igraph clusters graph_from_adjacency_matrix
#' @importFrom methods as
#' @export
#'
groupByGenet <-  function(sf, buffGenet,...){
  buffDat <- sf::st_buffer(sf, dist = buffGenet) ## buffers the data for the
  # focal year by the buff argument
  overlaps = sf::st_intersects(buffDat,buffDat) ## identifies which polygons
  # overlap with each other

  i <- rep(1:length(overlaps), lengths(overlaps))
  j <- factor(unlist(overlaps))
  tab <- Matrix::sparseMatrix(i = i, j = as.integer(j), x = TRUE,
                             dimnames = list(NULL, levels(j)))

  connects <- Matrix::tcrossprod(tab, boolArith = TRUE)
  group <- igraph::clusters(igraph::graph_from_adjacency_matrix(
    methods::as(connects, "lsCMatrix")))$membership

  tapply(overlaps, group, function(x) sort(unique(unlist(x))))
  rametIDs <- tapply(1:length(overlaps), group, toString)

  groupID <- rep(NA, nrow(sf))

  for(i in 1:nrow(sf)){
    for(j in 1:length(rametIDs)){
      vec <- as.numeric(#str_trim(
        unlist(strsplit(rametIDs[j], ",")))
      #)
      if(i %in% vec){
        groupID[i] <- j

      } else{
        next
      }
    }
  }
 return(groupID)
}

