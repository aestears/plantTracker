#' checkDat
#' @description A function that both checks and prepares a data.frame for use in
#' the \code{\link{trackSpp}} function.
#'
#' @details This function is used internally in \code{\link{trackSpp}} and other
#' plantTracker functions to check the 'dat' and 'inv' arguments and ensure that
#' they are in the correct format with the correct column names. [checkDat()]
#' can also be used independently to check that a data.frame is in the correct
#' format and contains the correct data to be used in the [trackSpp()] or other
#' plantTracker functions.
#'
#' @param dat An sf data.frame of the same format as
#' \code{\link{grasslandData}}. It must have columns that contain...
#' * a unique identification for each research site in character format
#' with no NAs (the default column name is "Site")
#' * species name in character format with no NAs (the default column
#' name is "Species")
#' * unique quadrat identifier in character format with no NAs (the default
#' column name is "Quad")
#' *  year of data collection in integer format with no NAs (the
#' default column name is "Year")
#' * an s.f 'geometry' column that contains a polygon or multipolygon data type
#' for each individual observation (the default column name is "geometry")
#' @param inv The name of each element of the list is a
#' quadrat name in 'dat', and the contents of that list element is a numeric
#' vector of all of the years in which that quadrat (or other unique spatial
#' area) was sampled. Make sure this is the years the quadrat was actually
#' sampled, not just the years that have data in the 'dat' argument! This
#' argument allows the function to differentiate between years when the quadrat
#' wasn't sampled and years when there just weren't any individuals of a species
#' present in that quadrat.
#' @param species An optional character string argument. Indicates
#' the name of the column in 'dat' that contains species name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Species" (default value is 'Species').
#' @param site An optional character string argument. Indicates
#' the name of the column in 'dat' that contains site name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Site" (default value is 'Site').
#' @param quad An optional character string argument. Indicates
#' the name of the column in 'dat' that contains quadrat name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Quad" (default is 'Quad').
#' @param year An optional character string argument. Indicates
#' the name of the column in 'dat' that contains data for year of sampling. It
#' is unnecessary to include a value for this argument if the column name is
#' "Year" (default is 'Year').
#' @param geometry An optional character string argument. Indicates
#' the name of the column in 'dat' that contains sf geometry data. It is
#' unnecessary to include a value for this argument if the column name is
#' "geometry" (default is 'geometry').
#' @param reformatDat A TRUE/FALSE argument. If 'FALSE', which is the default
#' value, then [checkDat()] prints a message that says the 'dat' and 'inv'
#' datasets are ready for use in the [trackSpp()] function. This message
#' includes any column name arguments that must be included in the [trackSpp()]
#' function call for these datasets.
#' If 'TRUE', [checkDat()] returns a list with three elements: a version of
#' 'dat' that has been checked for errors and is ready to be input directly into
#' the  \code{\link{trackSpp}} or  \code{\link{assign}} functions, a version of
#' 'inv' that is checked and ready for input into [trackSpp()] or [assign()],
#' and a character vector called 'userColNames' that has all of the user-defined
#' column names that were used in the original version of 'dat'.
#'
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return If the 'reformatDat' argument is FALSE (the default), the return of
#' this function is a message that says the 'dat' and 'inv' datasets are ready
#' for use in the [trackSpp()] function. This message includes any column name
#' arguments that must be included in the [trackSpp()] function call for these
#' datasets.
#' If the 'reformatDat' argument is TRUE, then [checkDat()] returns a list that
#' contains versions of 'dat' and 'inv' that are ready to go directly into the
#' [trackSpp()] or [assign()] functions. The list has the following elements:
#'
#'\item{dat}{An sf data.frame that has five columns called 'Species', 'Site',
#''Quad', 'Year', and 'geometry'. Any additional columns that were present in
#'the version of 'dat' input into the [checkDat()] function will also be
#'included in this version, although the column names will be appended with
#'"_USER".}
#'\item{inv}{A named list. The name of each element of the list is a quadrat
#' name in 'dat', and the contents of that list element is a numeric vector of
#' all of the years in which that quadrat (or other unique spatial area) was
#' sampled. For each list element, the vector of years is in sequential order
#' from oldest to most recent.}
#'\item{userColNames}{A named character vector. This vector contains the column
#'names provided by the user for all of the required data columns in the
#'original version of 'dat'. The name of each vector element indicates the type
#'of data that is contained in the column with the name in that vector element.
#'For example, if the 'dat' data.frame input into [checkDat()] has the names
#''Species','location' 'quadrat', 'Year', and 'geometry', then this list element
#'will be a character vector of these values with the names 'species', 'site',
#''quad', 'year', and 'geometry'. }
#'
#' @export
#'
#' @examples
#' checkDat(dat = grasslandData,
#' inv = grasslandInventory)
#'
#' checkDat(dat = grasslandData,
#' inv = grasslandInventory,
#' reformatDat = TRUE)
#'
#' @import sf
#'

## check that the dat and inv arguments contain the appropriate values/formats?
## is also reformatting the column names
checkDat <- function (dat, inv = NULL,
                      species = "Species",
                      site = "Site",
                      quad = "Quad",
                      year = "Year",
                      geometry = "geometry",
                      reformatDat = FALSE,
                      ...) {



  # work (all arg checking) ---------------------------------------------------
  ## check the datNames argument AND convert the names of the 'dat' argument to
  # be consistent with what this function expects

  ## put column name args. into a list for basic checks
  newNames <- list("species" = species, "site" = site, "quad" = quad,
                   "year" = year, "geometry" = geometry)
  ## check that each arg. is a character vector
  if (sum(sapply(newNames, is.character)) != 5) { ## if there is one or more
    # elements of the newNames list that is not a character vector
    ## find which elements are not character vectors
    badArgs <- paste("'",names( which(sapply(newNames, is.character) == FALSE)),
                     "'", collapse = ", and ")

    stop(paste0("The checkDat() argument(s) ", badArgs, " must each contain a
    single character string that gives the name(s) of the column(s) in 'dat'
    that contain the data for ", badArgs))

  } else { ## if each of the elements of 'newNames' is a character vector
    ## make sure that each of the elements of newNames is present as a column
    # name in 'dat'
    if (sum(unlist(newNames) %in% names(dat)) != 5) { ## if the column names of
      # 'dat' do NOT match the values provided in 'newNames'
      badBadArgs <- paste("'",names(newNames)[which(!unlist(newNames) %in%
                                                      names(dat))],"'",
                          collapse = ", and ")
      stop(paste0("The checkDat argument(s) ", badBadArgs, " contain values that
      are not column names in 'dat'. These arguments must be character vectors
      that give the name(s) of the column(s) in 'dat' that contain the data
                  for ",
      badBadArgs, ". Check for spelling errors, or make sure that you have
      included values for these arguments that give the name of the columns in
      'dat' that contain these data types." ))
    }
  }

  ## re-assign the names of dat to the default column names
  ## make a vector of 'user names'
  usrNames <- unlist(newNames)
  ## make a vector of 'default names'
  defaultNames <- c("Species", "Site", "Quad", "Year",
                    "geometry")

  ## add on a suffix to the column names that are not the column names that
  # must be specified in this function call
  otherNames <- names(dat)[!names(dat) %in% unlist(newNames)]

  ## change the names of the 'extra' columns so they don't get overwritten
  names(dat)[names(dat) %in% otherNames] <- paste0(
    names(dat)[!names(dat) %in% unlist(newNames)],
    "_USER"
  )

  ## replace the user-provided names in 'dat' with the default names
  names(dat)[match(usrNames, names(dat))] <- defaultNames

  ## proceed with remaining checks
  ## check the 'dat' argument (with default names)

  ## is the 'dat' argument in the correct format? (is it an 'sf' object of type
  # 'POLYGON' or 'MULTIPOLYGON'?)
  if (sum(st_is(dat, c("POLYGON", "MULTIPOLYGON"))) != nrow(dat)) {
    stop("'dat' is not in correct sf format.
         sfc must be POLYGON or MULTIPOLYGON")
  }
  ## does the 'dat' argument contain 'units'? If so, return an error, since
  # units make everything more complicated
  ## does the 'dat' argument contain any 'invalid' geometries?
  if (sum(!sf::st_is_valid(dat)) > 0) {
    inv_row <- which(sf::st_is_valid(dat)==FALSE)
      stop("'dat' contains an invalid geometry in row ", paste(inv_row, collapse = ", "), ". This issue must be resolved before 'dat' can be used in any plantTracker functions.")
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
        (!is.integer(dat$Year) &
        !is.numeric(dat$Year))## must be an integer vector
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

  if (is.null(inv) == FALSE) { ## if there is an argument for 'inv' (the default
    # is NULL, so if it is not NULL, then we need to check it)
    ## check the 'inv' argument
    ## is inv a list?
    if (is.list(inv) == TRUE &  ## 'inv' argument must be a list
        sum(!sapply(inv, is.numeric)) == 0 & ## each element of 'inv' must be a
        # numeric vector
        sum(!sapply(inv, length) > 0) == 0 ## each element must contain at least
        #one value
    ) {
      ## get unique quadrats in 'dat'
      datQuads <- unique(dat[,"Quad"]$Quad)
      if (sum(!datQuads %in% names(inv))==0) { ## do the quadrats in 'dat' have
        #corresponding inventory data in 'inv'? If yes:
        ## check that the years for each quadrat in 'dat' match the years for
        # each quadrat in 'inv'
        ## get the quadrat/year combos present in 'dat'
        datQuadYear <- unique(st_drop_geometry(dat[,c("Quad", "Year")]))
        datQuadYear<- paste0(datQuadYear$Quad,":",datQuadYear$Year)
        ## get the years present for each quad in 'inv'
        invYears <- unlist(inv)
        names(invYears) <- NULL
        ## make a data.frame that contains the quads included in 'inv' and the
        # years present in 'inv' for each quadrat
        invQuadYear <- data.frame(
          "Quad" = rep(names(inv), times = sapply(inv,length)),
          "Year" = invYears)
        invQuadYear <- paste0(invQuadYear$Quad, ":", invQuadYear$Year)
        ## compare the 'datQuadYear' vector and the 'invQuadYear' to make sure
        # that the years in 'dat' are represented in 'inv'
        if (sum(!datQuadYear %in% invQuadYear) != 0) {
          ## find those quads that have a mismatch between 'dat' and 'inv' (the
          # 'inv' data for this quad does not contain some/all years that are
          # present in 'dat' for this quad)
          misMatchQuads <- paste0("'",unique(sapply(strsplit(datQuadYear[
            which(!datQuadYear %in% invQuadYear)], ":"), function(x) x[1])),"'",
            collapse = ", and ")
          stop(paste0("Mismatch between years in 'dat' and years in 'inv' for
                      quadrat(s) ", misMatchQuads, ". The mismatch is for the
                      following quadrat/year combinations: ",
                      paste0(datQuadYear[which(!datQuadYear %in% invQuadYear)],
                             collapse = ", "), " . Either 'inv' does not contain
                             all the years in which these quadrats were
                      measured, or the years in 'dat' for these observations are
                      incorrect."))
        }
      } else { ## there is NOT data in 'inv' that corresponds to every quadrat
        # in 'dat'
        quadMissing <- paste0("'",datQuads[!datQuads %in% names(inv)],"'",
                              collapse = ", and ")
        stop(paste0("The 'inv' argument does not contain sampling year data for
quadrat(s) ", quadMissing, ", which have data in the 'dat' argument. The 'inv'
list must contain element(s) for each quadrat, and each must contain an integer
vector of years in which that quadrat was sampled."))
      }
    } else {
      stop("The 'inv' argument must be a list, and each element of that list
      must be a numeric vector with at least one value.")
    }
  }

## check to see if there are any 'duplicates' in the data (i.e. geometries that
                                             # are exact copies of one another)

  datDups <- dat
  datDups$centroidX <- sapply(suppressWarnings(sf::st_centroid(datDups)$geometry), function(x) x[1])
  datDups$centroidY <- sapply(suppressWarnings(sf::st_centroid(datDups)$geometry), function(x) x[2])
  datDups <- sf::st_drop_geometry(datDups)
  if (sum(duplicated(datDups)) > 0) {
    dup <- datDups[duplicated(datDups),]
    stop(paste0(
      "It seems that this row in 'dat' has the exact same values as another row(s): ",
                paste0(rownames(dup))))
  }

 ## check that the area of the 'geometry' column is >0 for each obs.
 # if there are areas that are 0, then give a 'warning'
  if (is.na(sf::st_crs(dat)) == FALSE) {
    if (sum(units::drop_units(sf::st_area(dat)) == 0) > 0) {
      badRows <- paste0(which(units::drop_units(sf::st_area(dat)) == 0),
                        collapse = ", ")
      warning(paste0("There are some observations in 'dat' that have an area of 0.
                   \n You should double-check the geometry for these rows to \n
                   make sure it's correct! Rows with 0 area: ",
                     badRows))
    }
  } else {
    if (sum(sf::st_area(dat) == 0) > 0) {
      badRows <- paste0(which(sf::st_area(dat) == 0),
                        collapse = ", ")
      warning(paste0("There are some observations in 'dat' that have an area of 0.
                   \n You should double-check the geometry for these rows to \n
                   make sure it's correct! Rows with 0 area: ",
                     badRows))
    }
  }


  # output ------------------------------------------------------------------
  ## prepare the output

  if (is.null(inv) == FALSE) { ## if there IS an argument for 'inv'
    if (reformatDat == TRUE) {
      ## return the 'dat' argument, with column names that are appropriate for
      # use directly in 'assign' or 'trackSpp'
      datReturn <- dat
      invReturn <- inv
      nameReturn <- usrNames
      ## return the reformatted data
      return(list(dat = datReturn,
                  inv = invReturn,
                  userColNames = nameReturn))

    } else if (reformatDat == FALSE) { ## if the user does not want the function
      # to return an output that is ready to go into the trackSpp function, but
      # just wants to know that their dataset is in the correct format
      ## make datReturn and invReturn empty so that the function has no output
      nameReturn <- NULL
      datReturn <- NULL
      invReturn <- NULL

      ## determine if there are any user-specified column names (that differ
      # from the defaults)
      if (sum(!usrNames %in% defaultNames) == 0) { ## if there are NO
        # differences in column names between the input 'dat' d.f. and the
        # defualt required names
        message("The data you put into the 'checkDat()' function for the 'dat' and
        'inv' arguments are ready to be used in the 'trackSpp()' function! You
        do not need to include any values for the 'species', 'site', 'quad',
        'year', and 'geometry' arguments in 'trackSpp()")

      } else if (sum(!usrNames %in% defaultNames) > 0) { ## if there ARE
        # differences in column names between the input 'dat' d.f. and the
        # defualt required names
        ## get the usrNames that are different than the default Names
        neededArgs <- newNames[!usrNames %in% defaultNames]
        message(paste0("The data you put into the 'checkDat()' function for the
        'dat' and 'inv' arguments are ready to be used in the 'trackSpp()'
        function! However, make sure that you include the character value(s): "
                   , paste0("'", unlist(neededArgs), "'", collapse = ", and "),
                   " in the corresponding ",paste0("'", names(neededArgs), "'",
                                                   collapse = ", and "),
                   " arguments"))
      }
    } else {
      stop("The 'reformatDat' argument must be logical (i.e. TRUE or FALSE).")
    }
  } else if (is.null(inv) == TRUE) { ## if there is not an input for 'inv'
    if (reformatDat == TRUE) {
      ## return the 'dat' argument, with column names that are appropriate for
      # use directly in 'assign' or 'trackSpp'
      datReturn <- dat
      nameReturn <- usrNames
      ## return the reformatted data
      return(list(dat = datReturn,
                  userColNames = nameReturn))

    } else if (reformatDat == FALSE) { ## if the user does not want the function
      # to return an output that is ready to go into the trackSpp function, but
      # just wants to know that their dataset is in the correct format
      ## make datReturn and invReturn empty so that the function has no output
      nameReturn <- NULL
      datReturn <- NULL

      ## determine if there are any user-specified column names (that differ
      # from the defaults)
      if (sum(!usrNames %in% defaultNames) == 0) { ## if there are NO
        # differences in column names between the input 'dat' d.f. and the
        # defualt required names
        message(
          paste0("The data you put into the 'checkDat()' function for the",
        "'dat' and 'inv'arguments are ready to be used in the 'trackSpp()'",
        "function! You do not need to include any values for the 'species',",
        "'site', 'quad', 'year', and 'geometry' arguments in 'trackSpp()"))
      } else if (sum(!usrNames %in% defaultNames) > 0) { ## if there ARE
        # differences in column names between the input 'dat' d.f. and the
        # defualt required names
        ## get the usrNames that are different than the default Names
        neededArgs <- newNames[!usrNames %in% defaultNames]
        message(
          paste0("The data you put into the 'checkDat()' function for the",
        "'dat' and 'inv' arguments are ready to be used in the 'trackSpp()'",
        "function! However, make sure that you include the character value(s): "
        , paste0("'", unlist(neededArgs), "'", collapse = ", and "), " in the ",
        "corresponding ",paste0("'", names(neededArgs), "'", collapse =
                                  ", and ")," arguments"))
      }
    } else {
      stop("The 'reformatDat' argument must be logical (i.e. TRUE or FALSE).")
    }
  }
}


# testing -----------------------------------------------------------------
#
# dat <- grasslandData
# names(dat)[1] <- "species"
# names(dat)[8] <- "quadrat"
#
# inv <- grasslandInventory
#
# datNames =  c(
#   "Species = species",
#   "Site = ",
#   "Quad = quadrat",
#   "Year = Year",
#   "geometry = geometry")
#
# checkDat(dat, inv,  quad = "location", year = "YeAr", reformatDat = TRUE)
