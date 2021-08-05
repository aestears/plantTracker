#' Create maps of quadrats over time, color-coded either by species or
#' by trackID
#'
#' @param dat
#' @param type
#' @param species
#' @param site
#' @param quad
#' @param year
#' @param geometry
#' @param trackID
#'
#' @return
#' @export
#'
#' @examples

drawQuadMap <- function (dat, type = "Species",
                         species = "Species",
                         site = "Site",
                         quad = "Quad",
                         year = "Year",
                         geometry = "geometry",
                         trackID = "trackID") {

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

  ## check the 'type' argument
  if (is.character(type) == TRUE) {
    if (type != "Species" &
        type != "trackID") {
      stop("'type' must have either the value 'Species' or 'trackID'")
    }
  } else {
  stop("'type' must be a single character argument with either the value
       'Species' or 'trackID'.")
  }

  ## check the 'trackID' argument (only if 'type' is 'trackID')
  if (type == "trackID") {
    if (is.character(trackID) == "TRUE") { ## make sure trackID is a chracter
      ## make sure that the value for 'trackID' is a column name in 'dat'
      if (sum(names(dat) %in% trackID) != 1) {
        stop("'trackID' must be the name of a column in 'dat'.")
      }
    } else {
      stop("'trackID' must be a single character argument that gives the name of
         the column in 'dat' that contains the unique identifying information
         for each genet in each year. ")
    }
    ## change the name of the column in 'dat' to the correct value
    names(dat)[which(names(dat) == trackID)] <- "trackID"
  }

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

  ## get the boundary box dimensions
  quadBox <- st_bbox(dat)

  ## get the number of years in the dat
  numYears <- length(unique(dat$Year)) + 1
  ## get the number of rows we need
  numRows <- ceiling(numYears/4)

  ## change the margins of the plots
  par(mar = c(1,1,1,1))
  ## change the size of the plotting area
  par(mfrow = c(4, numRows))

  ## plot the quadrat maps
  ## (for 'Species' method)
  if (type == "Species") {
    for (i in sort(unique(dat$Year))) {
      datTemp <- dat[dat$Year == i, ]
      plot(
        datTemp$geometry,
        col = as.factor(datTemp$Species),
        main = paste(unique(datTemp$Year)),
        xlim = c(quadBox$xmin, quadBox$xmax),
        ylim = c(quadBox$ymin, quadBox$ymax)
      )
      lines(
        x = c(quadBox$xmin, quadBox$xmin),
        y = c(quadBox$ymin, quadBox$ymax)
      )
      lines(
        x = c(quadBox$xmax, quadBox$xmax),
        y = c(quadBox$ymin, quadBox$ymax)
      )
      lines(
        x = c(quadBox$xmin, quadBox$xmax),
        y = c(quadBox$ymin, quadBox$ymin)
      )
      lines(
        x = c(quadBox$xmin, quadBox$xmax),
        y = c(quadBox$ymax, quadBox$ymax)
      )
    }
    ## make a 'plot' that is a legend
    plot(
      x = NULL,
      y = NULL,
      xlim = c(quadBox$xmin, quadBox$xmax),
      ylim = c(quadBox$ymin, quadBox$ymax),
      axes = FALSE
    )
    legend(
      x = quadBox$xmin,
      y = quadBox$ymax,
      legend = unique(datTemp$Species),
      fill = as.factor(unique(datTemp$Species)),
      bty = "n",
      cex = 1,
      x.intersp = .2
    )
  } else if (type == "trackID") {
    for (i in sort(unique(dat$Year))) {
      datTemp <- dat[dat$Year == i, ]
      plot(
        datTemp$geometry,
        col = as.factor(datTemp$trackID),
        main = paste(unique(datTemp$Year)),
        xlim = c(quadBox$xmin, quadBox$xmax),
        ylim = c(quadBox$ymin, quadBox$ymax)
      )
      lines(
        x = c(quadBox$xmin, quadBox$xmin),
        y = c(quadBox$ymin, quadBox$ymax)
      )
      lines(
        x = c(quadBox$xmax, quadBox$xmax),
        y = c(quadBox$ymin, quadBox$ymax)
      )
      lines(
        x = c(quadBox$xmin, quadBox$xmax),
        y = c(quadBox$ymin, quadBox$ymin)
      )
      lines(
        x = c(quadBox$xmin, quadBox$xmax),
        y = c(quadBox$ymax, quadBox$ymax)
      )
    }
    ## make a 'plot' that is a legend
    plot(
      x = NULL,
      y = NULL,
      xlim = c(quadBox$xmin, quadBox$xmax),
      ylim = c(quadBox$ymin, quadBox$ymax),
      axes = FALSE
    )
    legend(
      x = quadBox$xmin,
      y = quadBox$ymax,
      legend = unique(datTemp$trackID),
      fill = as.factor(unique(datTemp$trackID)),
      bty = "n",
      cex = 1,
      x.intersp = .2
    )
  }
}
# test --------------------------------------------------------------------

# dat <- dat[dat$Site == "CO" &
#            dat$Quad == "unun_11",]
