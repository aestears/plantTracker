#' This function will assign grouping variables for each feature of an sf
#' object. Grouping is based on buffers of a user-defined distance overlapping.
#'
#' @param sf your sf object
#' @param buffGenet the distance you want the buffer to be (put 2.5 to group
#' polygons that are 5cm from each other)
#'
#' @return
#' @export
#'
#' @examples##
#'
#'  Usage:
#' @import igraph, Matrix, sf
library(igraph)
library(Matrix)
library(sf)



groupByGenet <-  function(sf, buffGenet){
  buffDat <- sf::st_buffer(sf, dist = buffGenet) ## buffers the data for the
  # focal year by the buff argument
  overlaps = sf::st_intersects(buffDat,buffDat) ## identifies which polygons
  # overlap with each other

  i <- rep(1:length(overlaps), lengths(overlaps))
  j <- factor(unlist(overlaps))
  tab <- Matrix::sparseMatrix(i = i, j = as.integer(j), x = TRUE,
                             dimnames = list(NULL, levels(j)))

  connects <- Matrix::tcrossprod(tab, boolArith = TRUE)
  group <- igraph::clusters(graph_from_adjacency_matrix(
    as(connects, "lsCMatrix")))$membership

  tapply(overlaps, group, function(x) sort(unique(unlist(x))))
  trackIDs <- tapply(1:length(overlaps), group, toString)

  groupID <- rep(NA, nrow(sf))

  for(i in 1:nrow(sf)){
    for(j in 1:length(trackIDs)){
      vec <- as.numeric(#str_trim(
        unlist(strsplit(trackIDs[j], ",")))
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



