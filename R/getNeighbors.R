#' getNeighbors
#' @description
#'
#'
#' @param dat
#' @param buff
#' @param method
#' @param compType
#' @param focal
#' @param trackID
#' @param species
#' @param quad
#' @param year
#' @param site
#' @param geometry
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @details
#'
#' @return
#' @export
#'
#' @examples
#'
#'@import sf
getNeighbors <- function (dat, buff, method,
                          compType = 'allSpp',
                          trackID = 'trackID',
                          species = "Species",
                          quad = "Quad",
                          year = "Year",
                          site = "Site",
                          geometry = "geometry",
                          ...
                          ) {
  #dat ## a df that has been assigned demographic data by the 'trackSpp'
  # function. If compType == 'allSpp', 'dat' must include ALL species present in
  # the quadrat
  #buff # a numeric value that has the desired buff around each individual
  # for which neighborhood area/count will be calculated.
  #method # must equal either 'count' or 'area'. If 'count', then the number of
  # genets within the specified buffer will be counted. If 'area', then the area
  # within the specified buffer that is occupied will be calculated.
  ###AES question: for count method--should I use ramets or genets??
  #compType ## an argument that indicates whether you want to calculate only
  # conspecific local neighborhood (only look at individuals of the same
  # species -- 'oneSpp'), or heterospecific local neighborhood (look at
  # individuals of all species -- 'allSpp'). Default is 'allSpp'.
  #trackID ## the name of the column that contains the 'trackID'
  # (unique genet identifier) data

# argument checks ---------------------------------------------------------
## use the checkDat function to check the 'dat' function
checkedDat <- checkDat(dat = dat, species = species, site = site, quad = quad,
                       year = year, geometry = geometry, reformatDat = TRUE)

dat <- checkedDat$dat

## check/change trackID arg.
trackIDuserName <- paste0(trackID,"_USER")
trackIDdefaultName <- "trackID"
## change the name of the 'trackID' column to have the default column name
names(dat)[names(dat) == trackIDuserName] <- trackIDdefaultName

## check other args.
#buff
if (missing(buff)) {
stop("The 'buff' argument must have a value. This value idicates the buff
of the 'local neighborhood' around the focal individual, and must be a single
numeric value that is in the same units as the spatial attributes  of 'dat'.")
} else if (!is.numeric(buff) |  ## buff must be a  numeric integer
  buff < 0 | ## buff must be greater than or equal to 0
  length(buff)!=1 | ## buff must be a vector of length 1
  sum(buff > st_bbox(dat)[3:4]) > 0 ) { ##  must not be larger than the
  # largest values of the boundary box of 'dat' (an approximate way of checking
  # whether 'buff' and 'dat' have the same units)
stop("'buff' must be a single numeric value that is greater than zero, and in
the same units as the spatial attributes of 'dat' (i.e. in cm if the area of
individuals in 'dat' is recorded in cm^2")
  }

#method
if (missing(method)) {
stop("The 'method' argument must have a value. It must be a character vector
with either the value 'area' (meaning you want to calculate the total basal area
within the specified buffer of the focal individual that is occupied by other
individuals) or 'count' (meaning you want to know how many genets are within the
specified buffer of the focal individual")
} else if (is.character(method) & length(method) == 1) { ## must be a character
  # of length one
  ## make sure that the 'method' argument is lowercase
  method <- tolower(method)
  if (method != 'area' & method != 'count') {
   stop("The 'method' argument must have a value of 'area' or 'count'.")
  }
} else {
    stop("'method' must be a character vector of length one.")
  }

#compType
if (is.character(compType) & length(compType) == 1) { ## must be a
  # character of length one
  if (compType != 'allSpp' & compType != 'oneSpp') {
    stop("The 'compType' argument must have a value of 'allSpp' or 'oneSpp'.")
  }
} else {
  stop("'compType' must be a character vector of length one.")
}


# work --------------------------------------------------------------------
## assign a unique index to each row to simplify the looping process
dat$index <- 1:nrow(dat)


## subset the 'dat' d.f so that it only has the required columns
# (and dstore the rest)
## data to 'store'
datStore <- st_drop_geometry(dat[,!names(dat) %in% c("Species", "Site", "Quad", "Year",
                                    "trackID", "geometry")])
## trimmed 'dat' to use in the function
dat <- dat[,names(dat) %in% c("Species", "Site", "Quad", "Year",
                              "trackID", "geometry", "index")]


## make an empty column in 'dat' to contain the output neighborhood data
dat$neighbors <- NA

## put a buffer around each of the trackIDs in the entire data.frame
datBuffTemp <- st_buffer(x = dat, dist = buff)

## subtract the focal individuals from the buffered dataset
tempBuffGeometry <- mapply(FUN = st_difference, x = st_geometry(datBuffTemp),
                           y = st_geometry(dat))

## replace datBuffTemp geometry with the 'new' geometry ('tempBuffGeometry')
datBuff <- st_set_geometry(x = datBuffTemp, value = st_as_sfc(tempBuffGeometry))

## trim the buffers by the boundaries of the quadrat
datBuff <- st_set_geometry(x = datBuff,
                      value = st_intersection(st_as_sfc(st_bbox(dat)), datBuff))

if (method == 'count') {
  for (i in unique(dat$Site)) { ## loop through each site
    for (j in unique(dat[dat$Site== i ,"Quad"]$Quad)) { ## loop through each
      # quadrat
      for (k in unique(dat[dat$Site == i & dat$Quad == j, "Year"]$Year)) {
        ## loop through each year
        if (compType == 'oneSpp') { ## calculating conspecific neighborhood
          # (need to subset by species)
          for (l in unique(dat[dat$Site == i & dat$Quad == j & dat$Year == k,
                               "Species"]$Species)) {
            ## get the data for this site/quad/year/species combo
            datOneSp <- dat[dat$Site == i & dat$Quad == j &
                              dat$Year == k & dat$Species == l,]
            ## get the buffered data for this site/quad/year/species combo
            datOneBuff <- datBuff[datBuff$Site == i & datBuff$Quad == j &
                                    datBuff$Year == k & datBuff$Species == l,]

            ## calculate the neighborhood for each genet
            ## get a matrix of which polygons overlap each other
            overlapM <- st_intersects(datOneBuff, datOneSp, sparse= FALSE)
            ## make the diagonal of the matrix FALSE, because a genet can't
            # overlap with itself
            diag(overlapM) <- FALSE
            ## make a list such that the list element is the row index in
            # datOneBuff of the focal indvidual, and the values in each element
            # are the row indices of the polygons in datOneSp that overlap with
            # the focal individual
            overlapList <- apply(overlapM, MARGIN = 1, FUN = function(x)
              c(which(x==TRUE)))
            if (length(overlapList) > 0) {
              ## get the number of genets that overlap w/ the focal buffer
              datOneSp$neighbors <-  unlist(lapply(overlapList, length))
            } else {
              datOneSp$neighbors <- 0
            }
            ## put the neighbor counts into the 'dat' data.frame
            dat[match(datOneSp$index, dat$index),]$neighbors <-
              datOneSp$neighbors
          }
        } else if (compType == 'allSpp') { ## calculating heterospecific
          # neighborhood (don't need to subset by species)
          ## get the data for this site/quad/year combo
          datSpp<- dat[dat$Site == i & dat$Quad == j &
                         dat$Year == k ,]
          ## get the buffered data for this site/quad/year combo
          datSppBuff <- datBuff[datBuff$Site == i & datBuff$Quad == j &
                              datBuff$Year == k ,]
          ## calculate the number of neighbors for each genet
          ## get a matrix of which polygons overlap each other
          overlapM <- st_intersects(datSppBuff, datSpp, sparse= FALSE)
          ## make the diagonal of the matrix FALSE, because a genet can't
          # overlap with itself
          diag(overlapM) <- FALSE
          ## make a list such that the list element is the row index in
          # datSppBuff of the focal indvidual, and the values in each element
          # are the row indices of the polygons in datSpp that overlap with
          # the focal individual
          overlapList <- apply(overlapM, MARGIN = 1, FUN = function(x)
            c(which(x==TRUE)))

          ## put the neighbor counts into the 'datSpp' data.frame
          if (length(overlapList) > 0) {
            ## get the number of genets that overlap w/ the focal buffer
            datSpp$neighbors <-  unlist(lapply(overlapList, length))
          } else {
            datSpp$neighbors <- 0
          }

          ## put the neighbor counts into the 'dat' data.frame
          dat[match(datSpp$index, dat$index),]$neighbors <-
            datSpp$neighbors
          }
      }
    }
  }
} else if (method == "area") {
  ## get the overlapping polygon areas
  tempAreas <- suppressWarnings(
    st_intersection(x = datBuff, y = dat))
  if (compType == 'oneSpp') {
    ## match sites quads and years between buffered and raw data
      tempAreas2 <- tempAreas[tempAreas$Site == tempAreas$Site.1 &
                               tempAreas$Quad == tempAreas$Quad.1 &
                               tempAreas$Year == tempAreas$Year.1 &
                               tempAreas$Species == tempAreas$Species.1,]

      ## now aggregate by focal genet (column name = 'trackID')
      temp3 <- aggregate(tempAreas2$geometry,
                         by = list("Site" = tempAreas2$Site,
                                   "Quad" = tempAreas2$Quad,
                                   "Species" = tempAreas2$Species,
                                   "trackID" = tempAreas2$trackID,
                                   "Year" = tempAreas2$Year,
                                   "index" = tempAreas2$index),
                         FUN = function (x) sum(st_area(x)))

      ## proportionalize the area--divide the area of the 'neighbors' by the
      # area of the buffer
      ## make sure the do d.fs are in the same order (sort by index col.)
      datBuff <- datBuff[order(datBuff$index),]
      temp3 <- temp3[order(temp3$index),]
      ## calculate area of neighbors as a proportion of buffered area
      proportionalArea <- temp3$geometry/st_area(datBuff)
      ## put the 'area' data in the correct rows in 'dat'
      ## join by 'index' column
      dat[match(temp3$index, dat$index),]$neighbors <-
        proportionalArea

  } else if (compType == 'allSpp') {
    ## don't need to subset by species! (so use dat and datBuff)

    ## match sites quads and years between buffered and raw data
    tempAreas2 <- tempAreas[tempAreas$Site == tempAreas$Site.1 &
                             tempAreas$Quad == tempAreas$Quad.1 &
                             tempAreas$Year == tempAreas$Year.1,]

    ## now aggregate by focal genet (column name = 'trackID')
    temp3 <- aggregate(tempAreas2$geometry,
                       by = list("Site" = tempAreas2$Site,
                                 "Quad" = tempAreas2$Quad,
                                 "Species" = tempAreas2$Species,
                                 "trackID" = tempAreas2$trackID,
                                 "Year" = tempAreas2$Year,
                                 "index" = tempAreas2$index),
                       FUN = function (x) sum(st_area(x)))

    ## proportionalize the area--divide the area of the 'neighbors' by the
    # area of the buffer
    ## make sure the do d.fs are in the same order (sort by index col.)
    datBuff <- datBuff[order(datBuff$index),]
    temp3 <- temp3[order(temp3$index),]
    ## calculate area of neighbors as a proportion of buffered area
    proportionalArea <- temp3$geometry/st_area(datBuff)
    ## put the 'area' data in the correct rows in 'dat'
    ## join by 'index' column
    dat[match(temp3$index, dat$index),]$neighbors <-
      proportionalArea
  }
}

outputDat <-  dat

# output ------------------------------------------------------------------
## prepare the output
## change names of 'dat' back to the user-input values
## get the user-provided column names
userColNames <- c(checkedDat$userColNames[c(1:4)],
                  trackIDuserName,
                  checkedDat$userColNames[c(5)])

## rejoin the trackSppOut d.f with the 'extra' data stored in 'datStore'
outputDat <- merge(outputDat, datStore, by = "index")
## remove the 'index' value
outputDat <- outputDat[,names(outputDat) != "index"]

## re-name the appropriate columns in 'dat' data.frame with the
# user-provided names of 'dat'
## from above, user-provided names are stored in 'userColNames'
## make a vector of default column names
defaultNames <- c("Species", "Site", "Quad", "Year", "trackID", "geometry")

## reset the names for the columns that we changed to 'default' values
names(outputDat)[match(defaultNames, names(outputDat))] <- userColNames

## remove the '_USER' from the 'extra' column names
names(outputDat) <- gsub(names(outputDat),
                           pattern = "_USER", replacement = "")

## return
return(outputDat)
} ## end of 'getNeighbors()'


# testing -----------------------------------------------------------------
#
# dat <- grasslandData[grasslandData$Site == "CO" &
#                     grasslandData$Quad == "ungz_5a",]
# datIDs<- trackSpp(dat = dat, inv = grasslandInventory, dorm = 1, buff = 0.05,
#                   buffGenet = 0.005,
#                   clonal = data.frame("Species" = c("Bouteloua gracilis",
#                                                     "Agropyron smithii",
#                                                     "Sphaeralcea coccinea"),
#                   "clonal" = c(1,1,0)))
#
# names(datIDs)[c(3,4)] <- c("speciesName", "uniqueID")
#
# dataTest <- getNeighbors(dat = datIDs, buff = .15, method = "count",
#              compType = 'allSpp',
#              focal = 'genet',
#              species = "speciesName",
#              trackID = "uniqueID")
#
#
#
# plot(st_buffer(dataTest[dataTest$uniqueID == "AGRSMI_1997_13" &
#                           dataTest$Quad == "ungz_5a" & dataTest$Year == 1997
#                         ,]$geometry, .15), col = "pink")
# plot(dat[dat$Quad == "ungz_5a" & dat$Year == 1997,]$geometry, add = TRUE)
# plot(datIDs[datIDs$Quad == "ungz_5a" & datIDs$Year == 1997,]$geometry,
#      add = TRUE, col = as.factor(datIDs$uniqueID))
# plot(dataTest[dataTest$uniqueID == "AGRSMI_1997_13" ,]$geometry, col = "red", add = TRUE)
# labels <- datIDs[datIDs$Quad == "ungz_5a" & datIDs$Year == 1997,]
# st_centroid(labels$geometry)
# text(x = sapply(st_centroid(labels$geometry), FUN = function(x) x[1]),
#      y = sapply(st_centroid(labels$geometry), FUN = function(x) x[2]),
#      labels = labels$uniqueID)



