#' checkDat
#'
#' @return
#' @export
#'
#' @examples

## check that the dat and inv arguments contain the appropriate values/formats?
## is also reformatting the column names
checkDat <- function (dat, inv,
                      species = "Species",
                      site = "Site",
                      quad = "Quad",
                      year = "Year",
                      geometry = "geometry",
  reformatDat = FALSE,
  ...) {

# arguments ---------------------------------------------------------------
# dat ## an sf d.f that contains all of the raw digitized map data
  # (in grasslandData format) for as many sites and quads as you'd like.

#inv ## a list of the sampling years for each quadrat included in dat (in the
# same format as grasslandInventory).

# species/site/quad/year/geometry ## each arg. is a a character vector
# indicating the name in dat for the corresponding required columns
# (Species, Site, Quad, Year, sp_code_6, geometry)

# reformatDat ## a T/F argument. If 'T', then this function to output a list of
# data that is ready to go into the trackSpp function. If 'F', this function
# will output a message that indicates whether your dataset is ready for use in
# the trackSpp function.

# work --------------------------------------------------------------------
  ## check the datNames argument AND convert the names of the 'dat' argument to
  # be consistent with what this function expects

  ## put column name args. into a list for basic checks
  newNames <- list("species" = species, "site" = site, "quad" = quad,
                   "year" = year, "geometry" = geometry)
  ## check that each arg. is a character vector
  if (sum(sapply(newNames, is.character)) != 5) { ## if there is one or more
    # elements of the newNames list that is not a character vector
    ## find which elements are not character vectors
    badArgs <- paste(names(which(sapply(newNames, is.character) == FALSE)),
                     collapse = ", and ")

    stop(paste0("The argument(s) ", badArgs, " must each contain a single
    character string that gives the name(s) of the column(s) in 'dat' that
                contain the data for ", badArgs))

  } else { ## if each of the elements of 'newNames' is a character vector
    ## make sure that each of the elements of newNames is present as a column
    # name in 'dat'
    if (sum(unlist(newNames) %in% names(dat)) != 5) { ## if the column names of
      # 'dat' do NOT match the values provided in 'newNames'
      badBadArgs <- paste(names(newNames)[which(!unlist(newNames) %in%
                                                  names(dat))],
                          collapse = ", and ")
      stop(paste0("The argument(s) ", badBadArgs, " contain values that are not column names in 'dat'. These arguments must be character vectors that give the name(s) of the column(s) in 'dat' that contain the data for ", badBadArgs, ". Check for spelling errors." ))
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
  names(dat)[which(names(dat) %in% usrNames)] <- defaultNames


  ## proceed with remaining checks
  ## check the 'dat' argument (with default names)

  ## is the 'dat' argument in the correct format? (is it an 'sf' object of type
  # 'POLYGON' or 'MULTIPOLYGON'?)
  if(sum(st_is(dat, c("POLYGON", "MULTIPOLYGON"))) != nrow(dat)) {
    stop("'dat' is not in correct sf format.
         sfc must be POLYGON or MULTIPOLYGON")
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
         'Species'.")
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
      ## check that the years for each quadrat in 'dat' match the years for each
      # quadrat in 'inv'
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
        misMatchQuads <- paste0(unique(sapply(strsplit(datQuadYear[
          which(!datQuadYear %in% invQuadYear)], ":"), function(x) x[1])),
          collapse = ", and ")
        stop(paste0("Mismatch between years in 'dat' and years in 'inv' for quadrat(s) ", misMatchQuads, ". The mismatch is for the following quadrat/year combinations: ",
                    paste0(datQuadYear[which(!datQuadYear %in% invQuadYear)],
                           collapse = ", " ), " . Either 'inv' does not contain all the years in which these quadrats were measured, or the years in 'dat' for these observations are incorrect."))
      }
    } else { ## there is NOT data in 'inv' that corresponds to every quadrat
      # in 'dat'
      quadMissing <- datQuads[!datQuads %in% names(inv)]
      stop(paste0("The 'inv' argument does not contain sampling year data for quadrat ", quadMissing, ", which is present in the 'dat' argument. The 'inv' list must contain an element that contains the sampling year data for quadrat ", quadMissing))
    }
  } else {
    stop("The 'inv' argument must be a list, and each element of that list must
         be a numeric vector with at least one value.")
  }


# output ------------------------------------------------------------------
  ## prepare the output
  if (reformatDat == TRUE) {
    ## return the 'dat' argument, with column names that are appropriate for use
    # directly in 'assign' or 'trackSpp'
    datReturn <- dat
    invReturn <- inv
    nameReturn <- usrNames
  } else if (reformatDat == FALSE) { ## if the user does not want the function
    # to return an output that is ready to go into the trackSpp function, but
    # just wants to know that their dataset is in the correct format
    ## make datReturn and invReturn empty so that the function has no output
    nameReturn <- NULL
    datReturn <- NULL
    invReturn <- NULL

    ## determine if there are any user-specified column names (that differ
    # from the defaults)
    if (sum(!usrNames %in% defaultNames) == 0) { ## if there are NO differences
      # in column names between the input 'dat' d.f. and the defualt required
      # names
      print("The data you put into the 'checkDat()' function for the 'dat' and
      'inv'arguments are ready to be used in the 'trackSpp()' function! You do
      not need to include any values for the 'species', 'site', 'quad', 'year',
      and 'geometry' arguments in 'trackSpp()")

    } else if (sum(!usrNames %in% defaultNames) > 0) { ## if there ARE
      # differences in column names between the input 'dat' d.f. and the defualt
      # required names
      ## get the usrNames that are different than the default Names
      neededArgs <- newNames[!usrNames %in% defaultNames]
      print(paste0("The data you put into the 'checkDat()' function for the 'dat' and 'inv' arguments are ready to be used in the 'trackSpp()' function! However, make sure that you include the character value(s): ", paste0("'", unlist(neededArgs), "'", collapse = ", and "), " in the corresponding ",paste0("'", names(neededArgs), "'", collapse = ", and ")," arguments"))
    }
  } else {
    stop("The 'reformatDat' argument must be logical (i.e. TRUE or FALSE).")
  }

  return(list(dat = datReturn,
              inv = invReturn,
              nameReturn = usrNames))

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
#   "sp_code_6 = sp_code_6",
#   "geometry = geometry")
#
# checkDat(dat, inv,  quad = "location", year = "YeAr", reformatDat = TRUE)
