#' Group polygon observations together into a 'genet' based on proximity
#'
#' @description This function assigns a unique 'genetID' to polygons if they are
#' within a user-defined distance from one another. Each ramet is still
#' represented by a single row of data, but all ramets of the same genet have
#' the same 'genetID'.
#'
#' @details If polygons are 'grouped,' they are given the same number in the
#' 'genetID' column. This assignment is made using network analysis to group
#' together polygons that are closest to one another. In the context of
#' \pkg{plantTracker}, this function was designed to group ramets (vegetative
#' clones) together into one genet (genetically distinct individual). This
#' function was designed for use within the \code{\link{assign}} function (and
#' then the \code{\link{trackSpp}} function), but can be used independently.
#'
#' @param dat An sf object that contains polygons to be grouped. Typically
#' should include data only for one species, one quadrat, and one year. For
#' intended use, this dataset should be of the format described in
#' \code{\link{grasslandData}}. This function will run if 'dat' contains
#' only 'geometry' data, but it is *strongly* recommended that other columns
#' specified in \code{\link{grasslandData}} are included.
#' @param buffGenet A numeric argument indicating half of the maximum distance
#' between individual ramets (polygons) that the function will group together as
#' one genet. Note that this argument is in the same units that are used in the
#' 'dat' argument. For example, if buffGenet = 0.005 and we use the
#' \code{\link{grasslandData}} (in which measurements are in meters), then
#' polygons that are 0.01 m (1 cm) apart or less will be grouped as a genet.
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return A numeric vector of unique genetIDs that is as long as the number of
#' rows in 'dat.' Each element \code{i} in the output is the genetID for the
#' \code{i}th row in 'dat'.
#'
#' @examples
#' dat <- grasslandData[grasslandData$Site == "AZ" &
#' grasslandData$Quad == "SG2" &
#' grasslandData$Species == "Bouteloua rothrockii" &
#' grasslandData$Year == 1922,]
#'
#' groupByGenet(dat = dat, buffGenet = 0.01)
#' groupByGenet(dat = dat, buffGenet = 0.1)
#'
#' @seealso [assign()] and [trackSpp()], \pkg{plantTracker} functions that apply
#' this function across multiple species, quadrats, and years.
#'
#' @import Matrix
#' @importFrom igraph clusters graph_from_adjacency_matrix
#' @importFrom methods as
#' @export
#'
groupByGenet <-  function(dat, buffGenet,...){
  ## argument checks ---------------------------------------------------
  ## check the 'dat' argument
  ## is the 'dat' argument in the correct format? (is it an 'sf' object of type
  # 'POLYGON' or 'MULTIPOLYGON'?)
  if (sum(!st_is(dat, c("POLYGON", "MULTIPOLYGON"))) == 0) {
    ## is the name of the sf column containing data called 'geometry'?
    if (sum(names(dat) == "geometry") == 1) {
      ## does the 'geometry' column contain sf data?
      if (sum(!st_is(dat$geometry, c("POLYGON", "MULTIPOLYGON"))) != 0) {
        stop("The 'dat' data.frame must contain the spatial data in the column
               called 'geometry'")
      }
    } else {
      stop("The 'dat' data.frame must have its 'sf' data in a column called
             'geometry.'")
    }
  } else { ## if 'dat' is not in the correct sf format
    stop("'dat' is not in correct sf format.
         sfc must be POLYGON or MULTIPOLYGON")
  }

  ## if there is a column for 'species,' it must contain data for only one
  # species
  if (sum(names(dat) %in% "Species") == 1) {
    if (length(unique(dat$Species)) > 1  ## must be data for only one species
    ) {
      stop("The 'Species' column in 'dat' must contain only one species name,
           because 'dat' must be subset by species before being used in this
           function.")
    }
  } else { ## if the dat d.f. does not have a column called "Species"
    warning("It is recommended that the 'dat' data.frame contains a column
            called 'Species' that contains character vectors of species names.")
  }

  ## if there is a column for 'quad,' it must contain data for only one quadrat
  if (sum(names(dat) %in% "Quad") == 1) {
    if (length(unique(dat$Quad)) > 1  ## must be data for only one quadrat
    ) {
      stop("The 'Quad' column in 'dat' must contain only one quadrat name,
           because 'dat' must be subset by quadrat before being used in this
           function.")
    }
  } else { ## if the dat d.f. does not have a column called "Quad"
    warning("It is recommended that the 'dat' data.frame contains a column
            called 'Quad' that contains character vectors of quadrat names.")
  }

  ## if there is a column for 'Year,' it must contain data for only one year
  if (sum(names(dat) %in% "Year") == 1) {
    if (length(unique(dat$Year)) > 1  ## must be data for only one species
    ) {
      stop("The 'Year' column in 'dat' must contain only one year value, because
           'dat' must be subset by year before being used in this function.")
    }
  } else { ## if the dat d.f. does not have a column called "Year"
    warning("It is recommended that the 'dat' data.frame contains a column
            called 'Year' that contains character vectors of species names.")
  }

  ## check the 'buffGenet' argument
  if(!is.numeric(buffGenet) | ## buffGenet must be numeric
     buffGenet > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buffGenet must not
     # be larger than the dimensions of the quadrat
     buffGenet < 0 ## buffGenet must be greater than or equal to zero
  ) {
    stop("'buffGenet' argument must be numeric and cannot exceed the maximum
         dimensions of the quadrat")
  }

  ## work --------------------------------------------------------------
  buffDat <- sf::st_buffer(dat, dist = buffGenet) ## buffers the data for the
  # focal year by the buff argument
  overlaps = sf::st_intersects(buffDat,buffDat) ## identifies which polygons
  # overlap with each other

  i <- rep(1:length(overlaps), lengths(overlaps))
  j <- factor(unlist(overlaps))
  tab <- Matrix::sparseMatrix(i = i, j = as.integer(j), x = TRUE,
                              dimnames = list(NULL, levels(j)))

  connects <- Matrix::tcrossprod(tab, boolArith = TRUE)
  group <- igraph::clusters(igraph::graph_from_adjacency_matrix(
    methods::as(connects, "lMatrix")))$membership

  tapply(overlaps, group, function(x) sort(unique(unlist(x))))
  rametIDs <- tapply(1:length(overlaps), group, toString)

  groupID <- rep(NA, nrow(dat))

  for(i in 1:nrow(dat)){
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
  ## output -----------------------------------------------------------
  return(groupID)
}

