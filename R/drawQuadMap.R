#' Create maps of a quadrat over time
#'
#' @description This function creates maps of a quadrat over time, color-coded
#' by either species or by genet (trackID).
#'
#' @param dat An sf data.frame in which each row represents a unique polygon
#' (either a genet or a ramet) in a unique site/quadrat/year combination. It is
#' recommended that you only input data for one quadrat at a time. A data.frame
#' returned by \code{\link{trackSpp}} can be put into this function after being
#' subset by quadrat. 'dat' must have columns that contain...
#' * a unique identification for each research site in character format
#' with no NAs (the default column name is "Site")
#' * species name in character format with no NAs (the default column
#' name is "Species")
#' * unique quadrat identifier in character format with no NAs (the default
#' column name is "Quad")
#' *  year of data collection in integer format with no NAs (the
#' default column name is "Year")
#' * a unique identifier for each genet in character format with no NAs (the
#' default column name is "trackID")
#' * an s.f 'geometry' column that contains a polygon or multipolygon data type
#' for each individual observation (the default column name is "geometry")
#' @param type A character argument indicating how the plots returned by
#' `drawQuadMap()` will be color coded. If `type = "bySpecies"`, then
#' observations are color-coded by species. If `type = "bytrackID"`, then
#' observations are color-coded by trackID. The default value is "bySpecies".
#' @param addBuffer A logical argument indicating whether `drawQuadMap()`
#' will draw a small buffer around each polygon in the returned maps to make the
#' observations more visible. This is particularly useful for observations that
#' were originally mapped as points, which are hard to see when plotted in their
#' original dimensions. The buffer distance is 1/20th of the quadrat width. The
#' default value is `FALSE`.
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
#' @param trackID An optional character string argument. Indicates the name of
#' the column in 'dat' that contains unique identifiers for each genet. It is
#' unnecessary to include a value for this argument if the column name is
#' "trackID" (default is 'trackID')
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @import sf
#' @importFrom grDevices hcl.colors
#' @importFrom graphics legend lines par
#'
#' @return This function returns a multipanel plot where each panel shows a map
#' of the quadrat in a unique year. Panels are arranged in chronological order,
#' and plots are color-coded either by species or trackID (unique
#' genet identifier).
#'
#' #' @seealso [trackSpp()], which can be used to assign trackIDs
#' to observations.
#'
#' @export
#' @examples
#'  dat <- grasslandData[grasslandData$Site == c("AZ") &
#'  grasslandData$Quad == "SG2" &
#'  grasslandData$Year %in% c(1922:1925),]
#' names(dat)[1] <- "speciesName"
#' inv <- grasslandInventory[unique(dat$Quad)]
#' outDat <- trackSpp(dat = dat,
#'  inv = inv,
#'  dorm = 1,
#'  buff = .05,
#'  buffGenet = 0.005,
#'  clonal = data.frame("Species" = unique(dat$speciesName),
#'  "clonal" = c(TRUE,FALSE)),
#'  species = "speciesName",
#'  aggByGenet = TRUE
#'  )
#' drawQuadMap(dat = outDat,
#' type = "bySpecies",
#' addBuffer = FALSE,
#' species = "speciesName"
#' )
drawQuadMap <- function (dat,
                         type = "bySpecies",
                         addBuffer = FALSE,
                         species = "Species",
                         site = "Site",
                         quad = "Quad",
                         year = "Year",
                         geometry = "geometry",
                         trackID = "trackID",
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
                  badBadArgs, ". Check for spelling errors, or make sure that
                  you have included values for these arguments that give the
                  name of the columns in 'dat' that contain these data types."))
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
    if (type != "bySpecies" &
        type != "bytrackID") {
      stop("'type' must have either the value 'bySpecies' or 'bytrackID'")
    }
  } else {
    stop("'type' must be a single character argument with either the value
       'bySpecies' or 'bytrackID'.")
  }

  ## check the 'trackID' argument (only if 'type' is 'trackID')
  if (type == "bytrackID") {
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

  ## check the 'addBuffer' argument
  if (is.logical(addBuffer) == FALSE) {
    stop("The 'addBuffer' must be a logical argument. If 'addBuffer' = TRUE,
         then a small buffer will be drawn around individuals to make them
         easier to see. If 'addBuffer' = FALSE (the default), then a buffer will
         not be drawn.")
  }

  ## see how many quadrats have data in 'dat'
  if (length(paste0(unique(dat$Site), "_", unique(dat$Quad))) > 1) {
    stop("This function expects data from only one quadrat, and it looks like
         you've included data for multiple. Make sure to subset your 'dat'
         data.frame by quadrat, and use the 'drawQuadMap()' function
         independently on each quadrat's data.")
  }

  # work --------------------------------------------------------------------

  ## get the boundary box dimensions
  quadBox <- sf::st_bbox(dat)

  ## get the number of years in the dat
  numYears <- length(unique(dat$Year)) + 1
  ## get the number of rows we need
  numRows <- ceiling(numYears/4)

  ## save the previous margin and mfrow parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar)) # code line i + 1

  ## change the margins of the plots
  graphics::par(mar = c(1,1,1,1))
  ## change the size of the plotting area
  graphics::par(mfrow = c(numRows, 4))

  ## plot the quadrat maps
  ## (for 'Species' method)
  if (type == "bySpecies") {
    ## add a column to dat (that will later be removed) for color plotting
    ## make a data.frame with unique spp and colour data
    colorDat <- data.frame(Species = unique(dat$Species),
                           ## get color codes for full alpha colors
                           col = grDevices::hcl.colors(n = length(unique(dat$Species)),
                                            palette = "viridis"),
                           ## get color odes for alpha = .5 colors
                           colBuff = grDevices::hcl.colors(n = length(unique(dat$Species)),
                                                palette = "viridis",
                                                alpha = .5))
    ## add to the 'dat' data.frame
    dat$col <-
      colorDat$col[match(x = dat$Species, table = unique(dat$Species))]
    dat$colBuff <-
      colorDat$colBuff[match(x = dat$Species, table = unique(dat$Species))]
    if (addBuffer == FALSE) {
      for (i in sort(unique(dat$Year))) {
        datTemp <- dat[dat$Year == i, ]
        plot(
          datTemp$geometry,
          col = datTemp$col,
          main = paste(unique(datTemp$Year)),
          xlim = c(quadBox$xmin, quadBox$xmax),
          ylim = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmin),
          y = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmax, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymin)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymax, quadBox$ymax)
        )
      }
    } else if (addBuffer == TRUE) {
      ## calculate the buffer (based on the units of the data) -- 1/20th of the
      # quadrat width
      buffDist <- (quadBox$xmax)/20
      ## buffer the 'dat' data
      datBuff <- sf::st_buffer(x = dat, dist = buffDist)
      ## make the plots
      for (i in sort(unique(dat$Year))) {
        datTemp <- dat[dat$Year == i, ]
        datBuffTemp <- datBuff[datBuff$Year == i,]
        ## plot buffered data
        plot(datBuffTemp$geometry,
             col = datBuffTemp$colBuff,
             main = paste(unique(datBuffTemp$Year)),
             xlim = c(quadBox$xmin, quadBox$xmax),
             ylim = c(quadBox$ymin, quadBox$ymax)
        )
        plot(
          datTemp$geometry,
          col = datTemp$col,
          xlim = c(quadBox$xmin, quadBox$xmax),
          ylim = c(quadBox$ymin, quadBox$ymax),
          add = TRUE
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmin),
          y = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmax, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymin)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymax, quadBox$ymax)
        )
      }
    }
    ## make a 'plot' that is a legend
    plot(
      x = NULL,
      y = NULL,
      xlim = c(quadBox$xmin, quadBox$xmax),
      ylim = c(quadBox$ymin, quadBox$ymax),
      axes = FALSE
    )
    graphics::legend(
      x = quadBox$xmin,
      y = quadBox$ymax,
      legend = unique(dat$Species),
      fill = unique(dat$col),
      bty = "n",
      cex = 1,
      x.intersp = .2
    )
  } else if (type == "bytrackID") {
    ## add a column to dat (that will later be removed) for color plotting
    ## make a data.frame with unique trackID and colour data
    colorDat <- data.frame(trackID = unique(dat$trackID),
                           ## get color codes for full alpha colors
                           col = grDevices::hcl.colors(n = length(unique(dat$trackID)),
                                            palette = "viridis"),
                           ## get color odes for alpha = .5 colors
                           colBuff = grDevices::hcl.colors(n = length(unique(dat$trackID)),
                                                palette = "viridis",
                                                alpha = .5))
    ## add to the 'dat' data.frame
    dat$col <-
      colorDat$col[match(x = dat$trackID, table = unique(dat$trackID))]
    dat$colBuff <-
      colorDat$colBuff[match(x = dat$trackID, table = unique(dat$trackID))]
    ## make the plots
    if (addBuffer == TRUE) {## calculate the buffer (based on the units of the
      # data) -- 1/20th of the quadrat width
      buffDist <- (quadBox$xmax)/20
      ## buffer the 'dat' data
      datBuff <- sf::st_buffer(x = dat, dist = buffDist)
      ## make the plots
      for (i in sort(unique(dat$Year))) {
        datTemp <- dat[dat$Year == i, ]
        datBuffTemp <- datBuff[datBuff$Year == i,]
        plot(
          datBuffTemp$geometry,
          col = datBuffTemp$colBuff,
          main = paste(unique(datTemp$Year)),
          xlim = c(quadBox$xmin, quadBox$xmax),
          ylim = c(quadBox$ymin, quadBox$ymax)
        )
        plot(
          datTemp$geometry,
          col = datTemp$col,
          xlim = c(quadBox$xmin, quadBox$xmax),
          ylim = c(quadBox$ymin, quadBox$ymax),
          add = TRUE
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmin),
          y = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmax, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymin)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymax, quadBox$ymax)
        )
      }
    } else if (addBuffer == FALSE) {
      for (i in sort(unique(dat$Year))) {
        datTemp <- dat[dat$Year == i, ]
        plot(
          datTemp$geometry,
          col = unique(datTemp$col),
          main = paste(unique(datTemp$Year)),
          xlim = c(quadBox$xmin, quadBox$xmax),
          ylim = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmin),
          y = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmax, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymax)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymin)
        )
        graphics::lines(
          x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymax, quadBox$ymax)
        )
      }
    }

    ## make a 'plot' that is a legend
    plot(
      x = NULL,
      y = NULL,
      xlim = c(quadBox$xmin, quadBox$xmax),
      ylim = c(quadBox$ymin, quadBox$ymax),
      axes = FALSE
    )
    graphics::legend(
      x = quadBox$xmin,
      y = quadBox$ymax,
      legend = unique(dat$trackID),
      fill = unique(dat$col),
      bty = "n",
      cex = 1,
      x.intersp = .2
    )
  }
}
