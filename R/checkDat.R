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
    ###AES can put them all into a list to do basic checks (i.e. that they are are all character vectors) -- then check them against dat.
  reformatDat = FALSE, ##AES change the name of this? 'reformatted data'?
  ...) {

# arguments ---------------------------------------------------------------
# dat ## an sf d.f that contains all of the raw digitized map data
  # (in grasslandData format) for as many sites and quads as you'd like.

#inv ## a list of the sampling years for each quadrat included in dat (in the
# same format as grasslandInventory).

# species/site/quad/year/geometry ## each arg. is a a character vector
# indicating the name in dat for the corresponding required columns
# (Species, Site, Quad, Year, sp_code_6, geometry)

# reformatDat ## a T/F argument that indicates whether this function outputs a
# data.frame that is in the exact format (including column names) required by
# trackSpp, assign, and groupByGenet functions. If 'TRUE', then 'dat' and 'inv'
# are returned in this format. Default is FALSE


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
      stop(paste0("The argument(s) ", badBadArgs, " contain values that are not
      column names in 'dat'. These arguments must be character vectors that give
                  the name(s) of the column(s) in 'dat' that contain the data
                  for ", badBadArgs, ". Check for spelling errors." ))
    }
  }

  ## proceed with remaining checks
  ## check the 'dat' argument (with default names)

  ## re-assign the names of dat to the default column names
  ## assign an arbitrary index number to each row so we can re-join later
  dat$nameIndex <- c(1:nrow(dat))
  ## remove the 'extra' columns and store to rejoin with 'dat' later
  datStore <- dat[, !names(dat) %in% unlist(newNames)]
  dat <- dat[, names(dat) %in% c(unlist(newNames), "nameIndex") ]



  names(dat)[which(names(dat) %in% userDatNames)] <- defaultDatNames

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

  ## check the 'sq_code_6' column
  if (is.null(dat$sp_code_6) == FALSE) { ## does the 'sp_code_6' column exist?
    if (sum(is.na(dat$sp_code_6)) != 0 | ## cannot have 'NA' values for
        # sp_code_6
        !is.character(dat$sp_code_6)  ## must be a character vector
    ) {
      stop("The 'sp_code_6' column must be an character column with no 'NA's.")
    }
  } else {
    stop("The 'dat' data.frame must contain values in the column labeled
          'sp_code_6'.")
  }

  ## are there the same number of unique values for 'Species' and 'sp_code_6'?
  if (length(unique(dat$Species)) != length(unique(dat$sp_code_6))) {
    stop("The number of unique values in the 'Species' column and the 'sp_code_6' column do not match. Every species listed in 'dat' should have a six-letter code in the 'sp_code_6' column. Check spelling, and check that multiple species do not have the same six-letter code.")
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
        stop(paste0("Mismatch between years in 'dat' and years in 'inv' for quadrat(s) ", misMatchQuads, ". The mismatch is for the following quadrat/year combinations: ", paste0(datQuadYear[which(!datQuadYear %in% invQuadYear)], collapse = ", " ), " . Either 'inv' does not contain all the years in which these quadrats were measured, or the years in 'dat' for these observations are incorrect."))
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
  if (trackerFormat == TRUE) {
    ## return the 'dat' argument, with column names that are appropriate for use
    # directly in 'assign' or 'trackSpp'
    datReturn <- dat
    invReturn <- inv
  } else if (trackerFormat == FALSE) { ## if user does want an output, but
    # wants names of 'dat' to be the same as the input (names they provided)
    ## re-name the 'dat' input data.frame with the user-defined arguments
    names(dat)[which(names(dat) %in% defaultDatNames)] <- userDatNames
    datReturn <- dat
    if (printGoAhead == TRUE){
      print("The format of 'dat' and 'inv' are ready to be used in the 'trackSpp' or 'assign' functions, although you'll need to include the same 'datNames' argument that you included in this function.")
    }
  } else {
    stop("The 'trackerFormat' argument must be logical (i.e. TRUE or FALSE).")
  }

  return(list(dat = datReturn,
              inv = invReturn,
              userDatNames = userDatNames))

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
#   "Site = Site",
#   "Quad = quadrat",
#   "Year = Year",
#   "sp_code_6 = sp_code_6",
#   "geometry = geometry")
#
# checkDat(dat, inv, datNames)
