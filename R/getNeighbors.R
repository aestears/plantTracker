#' getNeighbors
#' @description
#'
#'
#' @param dat
#' @param buffer
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
getNeighbors <- function (dat, buffer, method,
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
  #buffer # a numeric value that has the desired buffer around each individual
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
#buffer
if (missing(buffer)) {
stop("The 'buffer' argument must have a value. This value idicates the buffer
of the 'local neighborhood' around the focal individual, and must be a single
numeric value that is in the same units as the spatial attributes  of 'dat'.")
} else if (!is.numeric(buffer) |  ## buffer must be a  numeric integer
  buffer < 0 | ## buffer must be greater than or equal to 0
  length(buffer)!=1 | ## buffer must be a vector of length 1
  sum(buffer > st_bbox(dat)[3:4]) > 0 ) { ##  must not be larger than the
  # largest values of the boundary box of 'dat' (an approximate way of checking
  # whether 'buffer' and 'dat' have the same units)
stop("'buffer' must be a single numeric value that is greater than zero, and in
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

###AES buffer entire d.f outside the for loops

for (i in unique(dat$Site)) { ## loop through each site
  for (j in unique(dat[dat$Site== i ,"Quad"]$Quad)) { ## loop through each
    # quadrat
    for (k in unique(dat[dat$Site == i & dat$Quad == j, "Year"]$Year)) { ## loop
      # through each year
      if (compType == 'oneSpp') { ## calculating conspecific neighborhood (need to
        # subset by species)
        for (l in unique(dat[dat$Site == i & dat$Quad == j & dat$Year == k,
                             "Species"]$Species)) {
          ## get the data for this site/quad/year/species combo
          datOneSp <- dat[dat$Site == i & dat$Quad == j &
                         dat$Year == k & dat$Species == l,]
             ## calculate the neighborhood for each genet
              ## loop through each unique trackID
              for (m in unique(datOneSp$trackID)) {
                ## make a buffer around the focal individual of the specified
                # 'buffer'
                rametBuff <- suppressWarnings(st_buffer(
                  x = datOneSp[datOneSp$trackID == m,],
                  dist = buffer))
                ## make into one large polygon
                rametBuff <- st_union(rametBuff)
                ## subtract the focal ramet from the buffer
                rametNeigh <- suppressWarnings(
                  st_difference(x = rametBuff,
                          y = st_union(datOneSp[datOneSp$trackID == m,])))
                ## compare 'rametNeigh' (without the focal ramet) to the rest of
                # the ramets in this quad (of this spp.) (make sure to remove
                # the overlap with the focal individual)
                if ( method == 'count') {
                  ## get the number of genets that overlap w/ the focal buffer
                  overlaps <- datOneSp[st_intersects(rametNeigh, datOneSp,
                                                   sparse = FALSE),]
                  ## remove focal individuals
                  overlaps <- overlaps[!overlaps$trackID %in%
                                         datOneSp[datOneSp$trackID == m,]$trackID,]

                  ## get the number of genets that are w/in the buffer
                  count <- length(unique(overlaps$trackID))
                  ## put this value into 'dat' in the 'neighbors' column
                  datOneSp[datOneSp$trackID == m,"neighbors"] <- count
                } else if (method == 'area') {
                  ## get the overlapping polygons
                  overlapArea <- suppressWarnings(
                    st_intersection(rametNeigh, datOneSp)[st_is(st_intersection(
                      rametNeigh, datOneSp), type = c("MULTIPOLYGON", "POLYGON"))
                      ,])
                  ## get the area of the ramets that are w/in the buffer
                  area <- sum(st_area(overlaps))
                  ## put this value into 'dat' in the 'neighbors' column for
                  # each ramet that belongs to the focal genet
                  datOneSp[datOneSp$trackID == m,"neighbors"] <- area
                }
              }
          ## append the 'datOneSp' d.f to the output d.f
          if (exists("outputDat") == FALSE) {
            outputDat <- datOneSp
          } else if (nrow(outputDat) > 0) {
            outputDat <- rbind(outputDat, datOneSp)
          }
        } ## end of loop to get data by species
      } else if (compType == 'allSpp') { ## calculating heterospecific
        # neighborhood (don't need to subset by species)
        ## get the data for this site/quad/year combo
        datSpp<- dat[dat$Site == i & dat$Quad == j &
                          dat$Year == k ,]
        ## calculate the neighborhood for each genet
          ## loop through each unique trackID
          for (m in unique(datSpp$trackID)) {
            ## make a buffer around the focal individual of the specified
            # 'buffer'
            rametBuff <- suppressWarnings(st_buffer(
              x = datSpp[datSpp$trackID == m,],
              dist = buffer))
            ## make into one large polygon
            rametBuff <- st_union(rametBuff)
            ## subtract the focal ramet from the buffer
            rametNeigh <- suppressWarnings(
              st_difference(x = rametBuff,
                            y = st_union(datSpp[datSpp$trackID == m,])))
            ## compare 'rametNeigh' (without the focal ramet) to the rest of
            # the ramets in this quad (of this spp.) (make sure to remove
            # the overlap with the focal individual)
###AES make the 'method' if/else into an internal function
            if ( method == 'count') {
              ## get the number of genets that overlap w/ the focal buffer
              overlaps <- datSpp[st_intersects(rametNeigh, datSpp,
                                               sparse = FALSE),]
              ## remove focal individuals
              overlaps <- overlaps[!overlaps$trackID %in%
                                     datSpp[datSpp$trackID == m,]$trackID,]

              ## get the number of genets that are w/in the buffer
              count <- length(unique(overlaps$trackID))
              ## put this value into 'dat' in the 'neighbors' column
              datSpp[datSpp$trackID == m, "neighbors"] <- count
            } else if (method == 'area') {
              ## get the overlapping polygons
              overlapArea <- suppressWarnings(
                st_intersection(rametNeigh, datSpp)[st_is(st_intersection(
                  rametNeigh, datSpp), type = c("MULTIPOLYGON", "POLYGON"))
                  ,])
              ## get the area of the ramets that are w/in the buffer
              area <- sum(st_area(overlapArea))
              ## put this value into 'dat' in the 'neighbors' column for
              # each ramet that belongs to the focal genet
              datSpp[datSpp$trackID == m, "neighbors"] <- area
            }
          }
        ## append results to the output d.f
        if (exists("outputDat") == FALSE) {
          outputDat <- datSpp
        } else if (nrow(outputDat) > 0) {
          outputDat <- rbind(outputDat, datSpp)
        }
      } ## end of loop to get heterospecific dataset (all species present)
    } ## end of loop to get data by year
    } ## end of loop to get data by quadrat
  } ## end of loop to get data by site

# output ------------------------------------------------------------------
## prepare the output
## change names of 'dat' back to the user-input values
## get the user-provided column names
userColNames <- c(checkedDat$userColNames[c(1:4)],
                  trackIDuserName,
                  checkedDat$userColNames[c(5)])

## rejoin the trackSppOut d.f with the 'extra' data stored in 'datStore'
outputDat <- merge(outputDat, datStore, by = "index")
## remove the 'indexStore' value
outputDat <- outputDat[,names(outputDat) != "index"]

## re-name the appropriate columns in 'dat' data.frame with the
# user-provided names of 'dat'
## from above, user-provided names are stored in 'userColNames'
## make a vector of default column names
defaultNames <- c("Species", "Site", "Quad", "Year",  "geometry", "trackID")

## reset the names for the columns that we changed to 'default' values
names(outputDat)[which(names(outputDat) %in% defaultNames)] <- userColNames

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
# names(datIDs)[c(1,6)] <- c("speciesName", "uniqueID")
#
# dataTest <- getNeighbors(dat = datIDs, buffer = .15, method = "area",
#              compType = 'allSpp',
#              focal = 'genet',
#              species = "speciesName",
#              trackID = "uniqueID")
#
#
# datTestTest <- testOut[testOut$Site == "CO" &
#                          testOut$location == "ungz_5a",]
#
# dataTest <- getNeighbors(dat = datTestTest, buffer = .15, method = "area",
#           compType = 'allSpp',
#           focal = 'genet',
#           species = "Species_Name",
#           quad = "location")
