#' Calculates the number of recruits of each species per year in each quadrat
#'
#' @param dat
#' @param species
#' @param quad
#' @param site
#' @param year
#' @param trackID
#' @param recruit
#' @param byGenet
#'
#' @return
#' @export
#'
#' @examples
getRecruits <- function(dat,
                        species = "Species",
                        quad = "Quad",
                        site = "Site",
                        year = "Year",
                        trackID = "trackID",
                        recruit = "recruit",
                        byGenet = TRUE
) {
  # argument checks ---------------------------------------------------------
  ## check the 'dat' data.frame and the column names (change if needed)
  newNames <- list("species" = species, "site" = site, "quad" = quad,
                   "year" = year, "trackID" = trackID)

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
  defaultNames <- c("Species", "Site", "Quad", "Year",
                    "trackID")

  ## replace the user-provided names in 'dat' with the default names
  names(dat)[match(usrNames, names(dat))] <- defaultNames

  ## proceed with remaining checks
  ## check the 'dat' argument (with default names)

  ## does 'dat' contain sf data?
  if (sum(class(dat) %in% "sf") > 0) { ## if there IS sf data, then drop it
    dat <- st_drop_geometry(dat)
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

  ## check the 'trackID' column
  if (is.null(dat$trackID) == FALSE) {## does the 'trackID' column exist?
    if (sum(is.na(dat$trackID)) != 0) {
      stop("The 'trackID' column must not have any 'NA's.")
    }
  } else {
    stop("The 'dat' data.frame must contain values in the column labeled
         'trackID'.")
  }

  ## check the 'recruit' column
  if (is.null(dat$recruit) == FALSE) {## does the 'recruit' column exist?
    if (nrow(dat[is.na(dat$recruit) |
                 dat$recruit == 1 |
                 dat$recruit == 0,]) != nrow(dat)
    ) { ## the values must be either NA, 1, or 0
      stop("The 'recruit' can only have values of NA, 1, or 0")
    }
  } else {
    stop("The 'dat' data.frame must contain values in the column labeled
         'recruit'.")
  }

  ## check the 'byGenet' argument
  if (is.null(byGenet) == FALSE) { ## if the 'byGenet' arg exists
    if (is.logical(byGenet) == FALSE) { ## if the arg isn't logical
      stop("The 'byGenet' argument must be either TRUE or FALSE. TRUE means that
      you want to count recruits so that each genet counts as one recruit.
      FALSE means that you want to count recruits so that each ramet counts
      as one recruit.")
    }
  } else { ## if the 'byGenet' arg doesn't exist
    stop("The function call must include a 'byGenet' argument that is either
         TRUE or FALSE.")
  }

  # work --------------------------------------------------------------------
  ## subset 'dat' to only get the individuals that are recruits
  datRecs <- dat[dat$recruit == 1 & is.na(dat$recruit) == FALSE,]

  ## drop all of the columns except those that you'll need for this fxn
  names(datRecs)
  datRecs <- datRecs[,
                     c("Site", "Quad", "Species", "trackID", "Year", "recruit")]

  if (byGenet == TRUE) {
    if (nrow(unique(datRecs[,c("Year","trackID")])) != nrow(datRecs) ) { ## if
      # the number of year/trackID combos is not the same as the number of rows
      # in 'datRecs', then 'datRecs' needs to be aggregated by genet
      ## make sure that each genet has only one row in each year
      datRecs<- aggregate(x = datRecs[,c('recruit')],
                          by = list("Site" = datRecs$Site,
                                    "Quad" = datRecs$Quad,
                                    "Species" = datRecs$Species,
                                    "trackID" = datRecs$trackID,
                                    "Year" = datRecs$Year),
                          FUN = sum)
      ## the output will make some values >1, so change them back to 1
      names(datRecs)[6] <- "recruit"
      datRecs[datRecs$recruits > 1,"recruit"] <- 1
    }
  }

  ## count the number of recruits in each year/species/quad/site combo
  datRecruits <- aggregate(x = datRecs[,c("recruit")], by = list(
    Year = datRecs$Year,
    Species = datRecs$Species,
    Quad = datRecs$Quad,
    Site = datRecs$Site
  ),
  FUN = length)

  names(datRecruits)[which(names(datRecruits) == "x")] <- "recruits"

  ## reorder the names of columns
  datRecruits <- datRecruits[,c("Site", "Quad", "Species", "Year", "recruits")]


  ## revert the names of the output data.frame to the names that the user input
  ## re-name the appropriate columns in the output data.frame with the
  # user-provided names of 'dat'
  ## from above, user-provided names are stored in 'usrNames'
  defaultNames <- defaultNames[1:4]
  ## reset the names for the columns that we changed to 'default' values
  names(datRecruits)[match(defaultNames, names(datRecruits))] <-
    usrNames[1:4]

  # output ------------------------------------------------------------------
  return(datRecruits)
}
