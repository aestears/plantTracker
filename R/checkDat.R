#' checkDat
#'
#' @return
#' @export
#'
#' @examples

## check that the dat and inv arguments contain the appropriate values/formats?
checkDat <- function (dat, inv, datNames =  c(
  "Species = Species",
  "Site = Site",
  "Quad = Quad",
  "Year = Year",
  "sp_code_6 = sp_code_6",
  "geometry = geometry"),
  trackerFormat = FALSE,
  inheritFromTrackSpp = FALSE,
  printGoAhead = TRUE,
  ...) {

# arguments ---------------------------------------------------------------
# dat ## an sf d.f that contains all of the raw digitized map data
  # (in grasslandData format) for as many sites and quads as you'd like.

#inv ## a list of the sampling years for each quadrat included in dat (in the
# same format as grasslandInventory).

# datNames ## a character vector indicating the names in dat for each of the
# required columns (Species, Site, Quad, Year, sp_code_6, geometry)

# trackerFormat ## a T/F argument that indicates whether this function outputs a
# data.frame that is in the exact format (including column names) required by
# trackSpp, assign, and groupByGenet functions. If 'TRUE', then 'dat' and 'inv'
# are returned in this format. Otherwise, dat and inv are not returned from this
# function (there is no output as long as there are no errors). Default is FALSE

# inheritFromTrackSpp ## a T/F argument that indicates whether 'dat' is being inherited rom 'trackSpp' function. If 'TRUE', then the error-checks of
# 'checkDat' are not conducted, since this function was already run in
# 'trackSpp' as part of the 'dat' argument-checking process. Default == FALSE.

# printGoAhead ## a T/F argument that indicates whether you want the function to
# tell you if the 'dat' and 'inv' arguments are ready to go into the 'trackSpp'
# and 'assign' function. Default is 'TRUE'.

# work --------------------------------------------------------------------
if (inheritFromTrackSpp == FALSE) {
  ## check the datNames argument AND convert the names of the 'dat' argument to
  # be consistent with what this function expects
  if(is.character(datNames) == TRUE & ## datNames must be a character arg.
     is.vector(datNames) == TRUE & ## must be a vector
     length(datNames) == 6 ## must have 6 values, one for each required column
  ){ ## if the datNames argument is a character vector of length 6,
    # proceed with the following testing
    if (sum(!sapply(regmatches(datNames,gregexpr(pattern = "=", datNames)),
                    length) %in% c(1,1,1,1,1,1)) == 0 ## must have an '=' in
        # each value of datNames, but only one '='
    ) {
      ## proceed with the following re-assignment of column names in dat
      ## separate the default from 'new' values
      datNamesTemp <- strsplit(datNames, "=")
      ## remove any spaces from default and new values
      datNamesTemp <- sapply(datNamesTemp, FUN = function(x) gsub(pattern = " ",
                                                                  x, replacement = ""))
      ## check that the datNamesTemp 'default' values contain the required names
      if (sum(!datNamesTemp[1,] %in% c("Species", "Site", "Quad", "Year",
                                       "sp_code_6", "geometry")) != 0) {
        stop("The first characters in each of the elements of the 'datNames'
             argument must be exactly--including case-- 'Species', 'Site',
             'Quad', 'Year', 'sp_code_6', and 'geometry'.")
      } else { ## if the 'default' values of datNames are correct...
        ## get a vector of 'default' names
        defaultDatNames <- datNamesTemp[1,]
        ## get a vector of 'new' names
        userDatNames <- datNamesTemp[2,]
      }
    } else {
      stop("The 'datNames' arg must have a single '=' in each value")
    }
  } else { ## if the datNames argument is NOT a character vector
    stop("The 'datNames' arg, if specified, must be a character vector and
         contain values for each of the required columns in dat ('Species =',
         'Site =', 'Quad =', 'sp_code_6 =', 'geometry =')")
  }

  ## check the 'dat' argument (with default names)
  ## are the user-defined column names from namDat arg. actually present in the
  # user-defined 'dat' argument?
  if(sum(!userDatNames %in% names(dat)[which(names(dat) %in%
                                             userDatNames)]) > 0) {
    ## get the name(s) of the user columns which are missing from the 'dat' arg.
    missingName <- userDatNames[which(userDatNames %in%
                                        names(dat)[which(names(dat) %in% userDatNames)]==0)]
    missingNameString <- paste(missingName, collapse = ", and ")
    ## stop the function and give an error
    stop(paste0("The column names ", missingNameString, " were provided in the 'datNames' argument of this function, but there are no columns in 'dat' called ", missingNameString, ". Either a column is missing in 'dat', or the column name in 'datNames' has been misspelled."))
  }

  ## re-assign the names of dat to the default
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
} else { ## if the 'inheritFromTrackSpp' arg. is TRUE
  ## make sure that the trackerFormat arg. is set to 'TRUE', since we want an
  # output dataset in the 'correct' format because we want an output that will
  # then be used in the 'assign' function
  trackerFormat <- TRUE
}

# output ------------------------------------------------------------------
  ## prepare the output
  if (trackerFormat == TRUE) {
    datReturn <- dat
  } else if (trackerFormat == FALSE) { ## if user does want an output, and wants names of 'dat' to be the same as the input
    ## re-name the 'dat' input data.frame with the user-defined arguments
    names(dat)[which(names(dat) %in% defaultDatNames)] <- userDatNames
    datReturn <- NULL
    if (printGoAhead == TRUE){
      print("The format of 'dat' and 'inv' are ready to be used in the 'trackSpp' or 'assign' functions, although you'll need to include the same 'datNames' argument that you included in this function.")
    }
  } else {
    stop("The 'trackerFormat' argument must be logical (i.e. TRUE or FALSE).")
  }

  return(datReturn)
}
