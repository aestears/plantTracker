#' Title
#'
#' @param dat
#' @param species
#' @param quad
#' @param site
#' @param year
#' @param geometry
#'
#' @return
#' @export
#'
#' @examples
getBasalAreas <- function(dat,
                     species = "Species",
                     quad = "Quad",
                     site = "Site",
                     year = "Year",
                     geometry = "geometry") {
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

    stop(paste0("The argument(s) ", badArgs, " must each contain a single character
string that gives the name(s) of the column(s) in 'dat' that contain the data
for ", badArgs))

  } else { ## if each of the elements of 'newNames' is a character vector
    ## make sure that each of the elements of newNames is present as a column
    # name in 'dat'
    if (sum(unlist(newNames) %in% names(dat)) != 5) { ## if the column names of
      # 'dat' do NOT match the values provided in 'newNames'
      badBadArgs <- paste("'",names(newNames)[which(!unlist(newNames) %in%
                                                      names(dat))],"'",
                          collapse = ", and ")
      stop(paste0("The argument(s) ", badBadArgs, " contain values that are not column
names in 'dat'. These arguments must be character vectors that give the name(s)
of the column(s) in 'dat' that contain the data for ", badBadArgs, ". Check for
spelling errors, or make sure that you have included values for these arguments that give the name of the columns in 'dat' that contain these data types." ))
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
    if (sum(!st_is(x = dat, type = c("POLYGON", "MULTIPOLYGON"))) != 0) {
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
  datArea<- aggregate(x = st_area(dat), by = list(
    Year = dat$Year,
    Species = dat$Species,
    Quad = dat$Quad,
    Site = dat$Site
  ),
  FUN = sum)

  ## rename the 'x' column to the correct value
  names(datArea)[which(names(datArea) == "x")] <- "absolute_basal_area"

  ## reorder the names of columns
  datArea <- datArea[,c("Site", "Quad", "Species", "Year", "absolute_basal_area")]

  ## calculate the percent of total basal area in that year
  # (the basal area of species A / basal area of quadrat occupied by any plants)
  ## get the total plant area by site/quad/year
  quadBasalArea <- aggregate(x = datArea$absolute_basal_area,
            by = list(
              Year = datArea$Year,
              Quad = datArea$Quad,
              Site = datArea$Site
            ),
            FUN = sum)

  names(quadBasalArea)[which(names(quadBasalArea) == 'x')] <- "percent_total_basal_area"

  datArea<- merge(x = datArea, y = quadBasalArea,
                  by = c("Site","Quad","Year"))

  ## calculate percentBasalArea
  datArea$percentBasalArea <- (datArea$absolute_basal_area /
                                      datArea$quad_basal_area)*100

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
