#' aggregateByGenet
#' @description This function aggregates a data.frame by unique trackIDs so that
#' each row represents a genet (unique genetic individual) in a given year.
#'
#' @details  This function is a [plantTracker]-specific, user-friendly wrapper
#' for the 'aggregate' function. This function was designed for use within the
#' \code{\link{trackSpp}} function, but can be also called independently by the
#' user. The function is designed to take a data.frame of the same format that
#' is returned by [trackSpp()], but it can be used to aggregate any dataset by
#' genet (or some other unique identifier), as long as there is a column that
#' identifies each distinct individual.
#'
#' @param dat An sf data.frame. Typically this is a data.frame that has been
#' returned from the [trackSpp()] function, although 'dat' can be any sf
#' data.frame of plant demographic data in which genet has a unique value
#' associated with it (what we call here a 'trackID'). If 'dat' already only has
#' one row for each unique trackID in each unique year (i.e. there are no
#' vegetative individuals--no ramets), then the output of [aggregateByGenet()]
#' will have the same number of rows as 'dat'. If 'dat' (in some cases) has
#' multiple rows (each of which is a ramet) for the same trackID in the same
#' year, then the output of [aggregateByGenet()] will have fewer rows than the
#' input 'dat' data.frame. 'dat' MUST have columns called 'basalArea_genet',
#' 'age', 'recruit', 'survives_tplus1', 'nearEdge', and 'size_tplus1', although
#' they can be populated with NAs. The [trackSpp()] function adds these columns,
#' so if you have made no changes to the data.frame that was returned from
#' [trackSpp()], then your data.frame should have these columns already.
#' @param site A character string giving the name of the column in 'dat' that
#' contains the values for site. The default is "Site".
#' @param quad A character string giving the name of the column in 'dat' that
#' contains the values for quadrat. The default is "Quad".
#' @param species A character string giving the name of the column in 'dat' that
#' contains the values for species. The default is "Species".
#' @param year A character string giving the name of the column in 'dat' that
#' contains the values for year. The default is "Year".
#' @param trackID A character string giving the name of the column in 'dat' that
#' contains the values for the unique identifier, or trackID. The default
#' is "trackID".
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return An sf data.frame that has been aggregated by genet so that each row
#' represents one unique genetic individual in one year. This data.frame has
#' columns containing data for site, quadrat, species, year, and trackID, as
#' well as columns called 'basalArea_genet', 'age', 'recruit',
#' 'survives_tplus1', 'nearEdge', and 'size_tplus1'. This data.frame will most
#' likely be shorter than the input 'dat', because a genet that was previously
#' broken into multiple rows representing each ramet in each year will now have
#' only one row in each year.
#'
#' @seealso This function is called inside the [trackSpp()] function.
#'
#' @export
#'
#' @examples
#'
aggregateByGenet <-  function(dat,
                              site = "Site",
                              quad = "Quad",
                              species = "Species",
                              year = "Year",
                              trackID = "trackID",
                              ...) {

# argument checking -------------------------------------------------------

#dat ## an sf d.f that contains data returned from the 'trackSpp' function.
  ## check that 'dat' is an s.f data.frame
  if (sum(!st_is(dat, type = c("POLYGON", "MULTIPOLYGON"))) != 0) {
stop("The 'dat' argument must be an sf data.frame with only 'POLYGON' or
'MULTIPOLYGON' geometries.")
  } else {
    ## check that the 'dat' d.f. has the required columns
    if (sum(c('basalArea_genet', 'age', 'recruit','survives_tplus1',
              'size_tplus1', 'nearEdge') %in% names(dat)) != 6) {
stop("The 'dat' argument must contain columns that are called 'basalArea_genet',
'age', 'recruit', 'survives_tplus1', 'nearEdge', and 'size_tplus1'.")
    }
  }

  ## check that the exact required columns in 'dat' ('basalArea_genet','age',
  # 'recruit', 'survives_tplus1', 'nearEdge' and 'size_tplus1') have the correct
  # data types
  ## 'basalArea_genet'
  if (is.numeric(dat$basalArea_genet)==FALSE |
      sum(!dat$basalArea_genet > 0) != 0 ) {
stop("The 'basalArea_genet' column in 'dat' must be a numeric vector with
         only positive values.")
  }
  ## 'age'
  if (is.numeric(dat$age) == FALSE | ## age must be numeric
      sum(dat$age < 0, na.rm = TRUE) != 0 ## age values must be positive
      ) {
stop("The 'age' column in 'dat' must be a numeric vector with positive values.")
  }
  ## 'recruit'
  if (is.numeric(dat$recruit) == FALSE | ## recruit must be numeric
      sum((dat$recruit == 0 | dat$recruit == 1 | is.na(dat$recruit))) !=
      nrow(dat)
      ## recruit values must be 1/0/NA
  ) {
stop("The 'recruit' column in 'dat' must be a numeric vector with values of only
1, 0, or NA.")
  }
  ## 'survives_tplus1'
  if (is.numeric(dat$survives_tplus1) == FALSE | ## survives_tplus1 must be numeric
      sum((dat$survives_tplus1 == 0 | dat$survives_tplus1 == 1 | is.na(dat$survives_tplus1))) !=
      nrow(dat)
      ## survives_tplus1 values must be 1/0/NA
  ) {
stop("The 'survives_tplus1' column in 'dat' must be a numeric vector with values
of only 1, 0, or NA.")
  }
  ## 'size_tplus1'
  if (is.numeric(dat$size_tplus1)==FALSE |
      sum(!dat$size_tplus1 > 0, na.rm = TRUE) != 0
      ) {
stop("The 'size_tplus1' column in 'dat' must be a numeric vector with only
positive values (or NA).")
  }
  ## 'nearEdge'
  if (is.logical(dat$nearEdge) == FALSE) {
stop("The 'nearEdge' column in 'dat' must be a logical vector.")
  }
#site ## the name of the column in 'dat' that contains the values for site

#quad ## the name of the column in 'dat' that contains the data for quadrat

#species ## the name of the column in 'dat' that contains the data for species

#year ## the name of the column in 'dat' that contains the data for Year

#trackID ## the name of the column in 'dat' that contains the data for trackID
  ## check these column name arguments
  ## put column name args. into a list for basic checks
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
spelling errors, or make sure that you have included values for these arguments
that give the name of the columns in 'dat' that contain these data types." ))
    }
  }

# work --------------------------------------------------------------------

## make sure the names of 'dat' are what the function expects
  ## re-assign the names of dat to the default column names
  ## make a vector of 'user names'
  usrNames <- unlist(newNames)
  ## make a vector of 'default names'
  defaultNames <- c("Species", "Site", "Quad", "Year",
                    "trackID")


  ## replace the user-provided names in 'dat' with the default names
  names(dat)[match(usrNames, names(dat))] <- defaultNames

  ## aggregate the 'dat' argument by trackID
  datOut <- aggregate(x = dat[,c('basalArea_genet', 'age', 'recruit',
                               'survives_tplus1', 'size_tplus1', 'nearEdge')],
            by = list("Site" = dat$Site,
                      "Quad" = dat$Quad,
                      "Species" = dat$Species,
                      "trackID" = dat$trackID,
                      "Year" = dat$Year),
            do_union = TRUE,
            FUN = mean)

  ## fix the 'nearEdge' mean issue--is averaged to a numeric value, not logical
  # b/c some of the ramets might be w/in the buffer from the quadrat edge while
  # others aren't.
  datOut[datOut$nearEdge > 0, 'nearEdge'] <- 1
  datOut$nearEdge <- as.logical(datOut$nearEdge)

  ## change the column name sback to what were present in 'dat'
  ## reset the names for the columns that we changed to 'default' values
  names(datOut)[match(defaultNames, names(datOut))] <-
    usrNames
# output ------------------------------------------------------------------
return(datOut)
  }
