#' Calculates basal area for each species in a quadrat
#'
#' @description This function calculates the total basal area for each species
#' present in a quadrat for each year of data collection. Both the absolute
#' basal area (in the same units of area as the input data.frame) is returned,
#' as well as the percentage of the total occupied basal area in the quadrat
#' that is occupied by a given species ("percent total basal area").
#'
#' @param dat An sf data.frame in which each row represents a unique polygon
#' (either a genet or a ramet) in a unique site/quadrat/year combination. A
#' data.frame returned by \code{\link{trackSpp}} can be put directly into this
#' function. However, it is not necessary for 'dat' to have demographic data or
#' unique identifiers (i.e. 'trackIDs') assigned. 'dat' must have columns that
#' contain a unique identification for each research site (default name is
#' "Site"), species name (default name is "Species"), quadrat identifier
#' (default name is "Quad"), year of data collection (default name is "Year"),
#' and an s.f 'geometry' column that contains a polygon or multipolygon data
#' type for each individual observation.
#' @param species An optional character string argument. Indicates
#' the name of the column in 'dat' that contains species name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Species" (default value is 'Species').
#' @param quad An optional character string argument. Indicates
#' the name of the column in 'dat' that contains quadrat name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Quad" (default is 'Quad').
#' @param site An optional character string argument. Indicates
#' the name of the column in 'dat' that contains site name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Site" (default value is 'Site').
#' @param year An optional character string argument. Indicates
#' the name of the column in 'dat' that contains data for year of sampling. It
#' is unnecessary to include a value for this argument if the column name is
#' "Year" (default is 'Year').
#' @param geometry An optional character string argument. Indicates
#' the name of the column in 'dat' that contains sf geometry data. It is
#' unnecessary to include a value for this argument if the column name is
#' "geometry" (default is 'geometry').
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @import sf
#' @return [getBasalAreas()] returns a data.frame with the columns "Site",
#' "Quad", "Year" and "Species". No two rows will have the same values for all
#' four of these columns. The column "absolute_basalArea" has the area of the
#' quadrat that is occupied by a species in a given unique site/quadrat/year
#' combination. This is in the same units as the area for area in 'dat'.
#' "quad_basalArea" gives the combined basal area of all plants in a given
#' site/quadrat/year. "percent_basalArea" gives the percentage of occupied basal
#' area within a quadrat that is occupied by each species in a given
#' site/quadrat/year. For example, species A has a basal area of 22 cm^2 in
#' quadrat 1 in 2005 ("absolute_basalArea = 22). In 2005, there are 50 cm^2 of
#' quadrat 1 that are occupied by plants ("quad_basalArea" = 55). 44% of the
#' occupied basal area in quadrat 1 in 2005 is occupied by species A
#' ("percent_basalArea" = 44).
#'
#' @export
#'
#' @examples
#' dat <- grasslandData[grasslandData$Site == "CO" &
#'  grasslandData$Year %in% c(1998:2002),]
#' names(dat)[1] <- "speciesName"
#' inv <- grasslandInventory[unique(dat$Quad)]
#' outDat <- trackSpp(dat = dat,
#'  inv = inv,
#'  dorm = 1,
#'  buff = .05,
#'  buffGenet = 0.005,
#'  clonal = data.frame("Species" = unique(dat$speciesName),
#'  "clonal" = c(1,0)),
#'  species = "speciesName",
#'  aggregateByGenet = TRUE
#'  )
#' getBasalAreas(dat = outDat,
#' species = "speciesName")
getBasalAreas <- function(dat,
                          species = "Species",
                          quad = "Quad",
                          site = "Site",
                          year = "Year",
                          geometry = "geometry",
                          ...) {
  # argument checks ---------------------------------------------------------
  ## check the 'dat' data.frame and the column names (change if needed)
  newNames <- list("species" = species, "site" = site, "quad" = quad,
                   "year" = year, geometry = "geometry")

  ## check that each arg. is a character vector
  if (sum(sapply(newNames, is.character)) != 5) { ## if there is one or more
    # elements of the newNames list that is not a character vector
    ## find which elements are not character vectors
    badArgs <- paste("'",names( which(sapply(newNames, is.character) == FALSE)),
                     "'", collapse = ", and ")

    stop(paste0("The argument(s) ", badArgs, " must each contain a single
    character string that gives the name(s) of the column(s) in 'dat' that
    contain the data for ", badArgs))

  } else { ## if each of the elements of 'newNames' is a character vector
    ## make sure that each of the elements of newNames is present as a column
    # name in 'dat'
    if (sum(unlist(newNames) %in% names(dat)) != 5) { ## if the column names of
      # 'dat' do NOT match the values provided in 'newNames'
      badBadArgs <- paste("'",names(newNames)[which(!unlist(newNames) %in%
                                                      names(dat))],"'",
                          collapse = ", and ")
      stop(paste0("The argument(s) ", badBadArgs, " contain values that are not
      column names in 'dat'. These arguments must be character vectors that give
      the name(s) of the column(s) in 'dat' that contain the data for ",
      badBadArgs, ". Check for spelling errors, or make sure that you have
      included values for these arguments that give the name of the columns in
      'dat' that contain these data types." ))
    }
  }

  ## re-assign the names of dat to the default column names
  ## make a vector of 'user names'
  usrNames <- unlist(newNames)
  ## make a vector of 'default names'
  defaultNames <- c( "Species", "Site", "Quad", "Year", "geometry")

  ## replace the user-provided names in 'dat' with the default names
  names(dat)[match(usrNames, names(dat))] <- defaultNames

  ## proceed with remaining checks
  ## check the 'dat' argument (with default names)

  ## does 'dat' contain sf data?
  if (sum(class(dat) %in% "sf") > 0) { ## if there IS sf data
    if (sum(!sf::st_is(x = dat, type = c("POLYGON", "MULTIPOLYGON"))) != 0) {
      ## if the sf data are not polygon or multipolygon ...
      stop("The sf data in 'dat' must be of either class POLYGON
           or MULTIPOLYGON")
    }
  } else { ## error if 'dat' does not have sf data
    stop("The 'dat' data.frame must contain sf spatial data for each
         observation.")
  }

  ## check the Species column
  if (is.null(dat$Species) == FALSE) {
    if (sum(is.na(dat$Species)) != 0 | ## cannot have 'NA' values for species
        !is.character(dat$Species)  ## must be a character vector
    ) {
      stop("The 'Species' column must be a character column with no 'NA's.")
    }
  } else {
    stop("The 'dat' data.frame must contain values in the column labeled
         'Species'."
    )
  }

  ## check the 'Year' column
  if (is.null(dat$Year) == FALSE) { ## does the 'Year' column exist?
    if (sum(is.na(dat$Year)) != 0 | ## cannot have 'NA' values for Year
        !is.integer(dat$Year) ## must be an integer vector
    ) {
      stop("The 'Year' column must be an integer column with no 'NA's.")
    }
  } else {
    stop("The 'dat' data.frame must contain values in the column labeled
          'Year'.")
  }

  ## check the 'Quad' column
  if (is.null(dat$Quad) == FALSE) { ## does the 'Quad' column exist?
    if (sum(is.na(dat$Quad)) != 0 | ## cannot have 'NA' values for Quad
        !is.character(dat$Quad) ## must be a character vector
    ) {
      stop("The 'Quad' column must be an character column with no 'NA's.")
    }
  } else {
    stop("The 'dat' data.frame must contain values in the column labeled
         'Quad'.")
  }

  ## check the 'geometry' column
  if (is.null(dat$geometry) == TRUE) { ## does the 'geometry' column exist?
    stop("The 'dat' data.frame must contain values in the column
         labeled 'geometry'")
  }

  # work --------------------------------------------------------------------

  ## drop all of the columns except those that you'll need for this fxn
  dat <- dat[,c("Site", "Quad", "Species", "Year", "geometry")]

  ## sum the area in each year/species/quad/site combo
  datArea <- aggregate(x = sf::st_area(dat), by = list(
    Year = dat$Year,
    Species = dat$Species,
    Quad = dat$Quad,
    Site = dat$Site
  ),
  FUN = sum)

  ## rename the 'x' column to the correct value
  names(datArea)[which(names(datArea) == "x")] <- "absolute_basalArea"

  ## reorder the names of columns
  datArea <- datArea[,c("Site", "Quad", "Species", "Year",
                        "absolute_basalArea")]

  ## calculate the percent of total basal area in that year
  # (the basal area of species A / basal area of quadrat occupied by any plants)
  ## get the total plant area by site/quad/year
  quadBasalArea <- aggregate(x = datArea$absolute_basalArea,
                             by = list(
                               Year = datArea$Year,
                               Quad = datArea$Quad,
                               Site = datArea$Site
                             ),
                             FUN = sum)

  names(quadBasalArea)[which(names(quadBasalArea) == 'x')] <- "quad_basalArea"

  datArea<- merge(x = datArea, y = quadBasalArea,
                  by = c("Site","Quad","Year"))

  ## calculate percentBasalArea
  datArea$percent_basalArea <- (datArea$absolute_basalArea /
                                  datArea$quad_basalArea)*100

  ## revert the names of the output data.frame to the names that the user input
  ## re-name the appropriate columns in the output data.frame with the
  # user-provided names of 'dat'
  ## from above, user-provided names are stored in 'usrNames'
  defaultNames <- defaultNames[1:4]
  ## reset the names for the columns that we changed to 'default' values
  names(datArea)[match(defaultNames, names(datArea))] <-
    usrNames[1:4]

  # output ------------------------------------------------------------------
  return(datArea)
}
