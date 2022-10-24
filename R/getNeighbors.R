#' Calculates local neighborhood density around each individual in a
#' mapped dataset
#'
#' @description This function calculates the density of individuals around
#' each distinct individual in a mapped dataset. It is intended for use on
#' a dataset that has been returned by the \code{\link{trackSpp}} function, but
#' will work with any sf data.frame where each individual row represents a
#' distinct individual (genet) in a unique year.
#'
#' @details This function draws a buffer around each individual genet of a
#' distance specified by the user. Then it either counts the number of other
#' genets within that buffer, or calculates the proportion of that buffer area
#' that is occupied by other individuals. [getNeighbors] can calculate
#' either interspecific local neighborhood density (between the focal individual
#' and all other individuals within the buffer area) or intraspecific local
#' neighborhood density (between the focal individual and other individuals of
#' the same species within the buffer area).
#'
#' @param dat An sf data.frame. Each row must represent a unique individual
#' organism in a unique year. This argument can be a data.frame that is returned
#' by the \code{\link{trackSpp}} function. It must have columns that contain a
#' unique identification for each research site (default name is "Site"),
#' species name (default name is "Species"), quadrat identifier (default name is
#' "Quad"), year of data collection (default name is "Year"), genet identity
#' (default name is "trackID"), and an s.f 'geometry' column that contains a
#' polygon or multipolygon data type for each individual observation.
#' @param buff A numeric value that is greater than or equal to zero. This
#' indicates the distance (in the same units as the spatial data in 'dat')
#' around each focal individual within which you want to look for competitors.
#' OR buff can be a data.frame with the columns "Species" and "buff". This
#' data.frame must have a row for each unique species present in 'dat', with
#' species name as a character string in the "Species" column, and a numeric
#' value in the 'buff' column you'd like to use for that species.
#' @param method A character string, either 'count' or 'area'. This argument
#' determines which method is used to calculate local neighborhood density.
#' If 'method' = 'count', then the number of other individuals within the buffer
#' will be returned in the 'neighbors' column. If method = 'area', then the
#' proportion of the buffer area that is occupied by other individuals will be
#' returned in the 'neighbors' column. If the data in 'dat' was mapped initially
#' as points, it is best to use 'method' = 'count'. If the data was mapped as
#' polygons that are representative of individual basal area, using
#' 'method' = 'area' is likely a more accurate representation of the crowding
#' the focal individual is experiencing.
#' @param compType A character string, either 'allSpp' or 'oneSpp'.
#' If compType = 'allSpp', then local neighborhood density is calculated
#' considering all individuals around the focal individual, regardless of
#' species (interspecific competition). If compType = 'oneSpp', then local
#' neighborhood density is calculated considering only individuals in the buffer
#' area that are the same species as the focal individual (intraspecific
#'  competition).
#' @param output A character string, either 'summed' or 'bySpecies'. The default
#' is 'summed'. This argument is only important to consider if you are using
#' compType = 'allSpp'. If output = 'summed', then only one count/area value is
#' returned for each individual. This value is the total count or area of all
#' neighbors within the focal species buffer zone, regardless of species. If
#' output = 'bySpecies', there is a count or area value returned for each
#' species present in the buffer zone. For example, there are 15 individuals
#' inside a buffer zone. Five are species A, three are species B, and 7 are
#' species C. If output = 'summed', then the 'neighbors_count' column in the
#' output data.frame will have the single value '15' in the row for this focal
#' individual. However, if output = 'bySpecies', the row for this focal
#' individual in the output data.frame will contain a named list in the
#' 'neighbors_count' column that looks like the one below. If 'method' = 'area'
#' and 'output' = 'bySpecies', a similar list will be returned, but will be in
#' the 'neighbors_area' column and will contain areas rather than counts.
#' ```{r}
#' list("Species A "= 5, "Species B" = 3, "Species C" = 7)
#' ```
#' @param trackID An optional character string argument. Indicates the name of
#' the column in 'dat' that contains a value that uniquely identifies each
#' individual/genet. It is unnecessary to include a value for this argument if
#' the column name is 'trackID' (the default value is 'trackID').
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
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return This function returns a data.frame with the same number of rows as
#' 'dat' and all of the same columns, but with an additional column or columns.
#' If method = 'count', then a column called "neighbors_count" is added, which
#' contains a count of the number of individuals within the buffer area that is
#' occupied by other individuals. If method = 'area', two columns are added. The
#' first is called "nBuff_area", which contains the area of the buffer around
#' the focal individual. The second is called "neighbors_area", which contains
#' the area of the individuals within the buffer zone around the
#' focal individual.
#' @seealso The [trackSpp()] function returns a data.frame that can be input
#' directly into this function. If a data.frame is not aggregated by genet such
#' that each unique genet/year combination is represented by only one row, then
#' passing it through the [aggregateByGenet()] function will return a data.frame
#' that can be input directly into this function.
#'
#' @examples
#' dat <- grasslandData[grasslandData$Site == c("AZ") &
#'  grasslandData$Species %in% c("Bouteloua rothrockii",
#'  "Calliandra eriophylla"),]
#' names(dat)[1] <- "speciesName"
#' inv <- grasslandInventory[unique(dat$Quad)]
#' outDat <- trackSpp(dat = dat,
#'  inv = inv,
#'  dorm = 1,
#'  buff = .05,
#'  buffGenet = 0.005,
#'  clonal = data.frame("Species" = unique(dat$speciesName),
#'  "clonal" = c(TRUE)),
#'  species = "speciesName",
#'  aggByGenet = TRUE
#'  )
#'
#'  finalDat <- getNeighbors(dat = outDat,
#'  buff = .15,
#'  method = 'count',
#'  compType = 'oneSpp',
#'  species = "speciesName")
#'
#'@import sf
#'
#'@export


getNeighbors <- function (dat, buff, method,
                          compType = 'allSpp',
                          output = 'summed',
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
  # intraspecific local neighborhood (only look at individuals of the same
  # species -- 'oneSpp'), or interspecific local neighborhood (look at
  # individuals of all species -- 'allSpp'). Default is 'allSpp'.
  #trackID ## the name of the column that contains the 'trackID'
  # (unique genet identifier) data

   # argument checks ---------------------------------------------------------
  ## use the checkDat function to check the 'dat' function
  checkedDat <- checkDat(dat = dat, species = species, site = site, quad = quad,
                         year = year, geometry = geometry, reformatDat = TRUE)

  dat <- checkedDat$dat

  ## check/change trackID arg.
  ## does dat have a column called 'trackID'?
  if (sum(names(dat) %in% paste0(trackID,"_USER")) > 0) {
    trackIDuserName <- paste0(trackID,"_USER")
    trackIDdefaultName <- "trackID"
    ## change the name of the 'trackID' column to have the default column name
    names(dat)[names(dat) == trackIDuserName] <- trackIDdefaultName
    ## make sure the trackID col is in the correct format
    if (sum(is.na(dat$trackID)) != 0) { ## make sure that there are no NAs in
      # the trackID col.
      stop("The column in 'dat' that contains trackID information cannot have
      any NA values")
    }
  } else {
    stop("The 'dat' argument must have a column that contains a unique
         identifier for each genetic individual (i.e. a 'trackID'). This column
         must have the same name that you specified in the 'trackID' argument in
         this function call.  ")
  }

  ## check other args.
  #buff
  if(missing(buff)) {
    stop("The 'buff' argument must have a value.")
  } else {
    if (is.numeric(buff) & length(buff) == 1) { ## is the value of buff a single
      # numeric value?
      if (buff < 0 | ## buff must be greater than or equal to 0
          buff > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buff must
          # not be larger than the dimensions of the quadrat
          length(buff)!=1) { ## buff must be a vector of length 1
        stop("If 'buff' is not specified for every species, it must be a single
        numeric value that is greater than or equal to 0")
      }
      ## make the 'buff' arg. into a data.frame to make things easier
      buff <- data.frame("Species" = unique(dat$Species), "buff" = buff)
    } else if (is.data.frame(buff)) {
      if (sum(!names(buff) %in% c("Species", "buff")) == 0) {
        if(sum(!unique(dat$Species) %in% buff$Species) > 0 | ## buff must have
           # data for all species
           sum(is.na(dat$buff)) > 0 | ## can't have NA values in buff
           !is.numeric(buff$buff) | ## can't have non-numeric values for
           # buff$buff
           sum(buff$buff < 0) > 0 | ## can't be less than 0
           round(buff$buff) != buff$buff ## must be whole numbers
        ) {
          stop("If the 'buff' argument is specified by species, it must be a
          data.frame that includes a 'Species' column with a row for every
          species in 'dat', and a 'buff' column that contains positive, numeric
          values for each species with no NAs.")
        }
      } else {
        stop("If the 'buff' argument is specifed by species, the column names
        must be 'Species' and 'buff'")
      }
    } else {
      stop("The 'buff' argument must be either a single numeric value that is
      greater than or equal to 0, OR a data.frame that has a 'Species' column
      with values for each species in 'dat', and a 'buff' column with numeric
      values for each species.")
    }
  }

  #method
  if (missing(method)) {
    stop("The 'method' argument must have a value. It must be a character vector
with either the value 'area' (meaning you want to calculate the total basal area
within the specified buffer of the focal individual that is occupied by other
individuals) or 'count' (meaning you want to know how many genets are within the
specified buffer of the focal individual")
  } else if (is.character(method) & length(method) == 1) { ## must be a
    # character of length one
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

  #output
  if (is.character(output) & length(output) == 1) { ## must be a
    # character of length one
    if (output != 'summed' & output != 'bySpecies') {
      stop("The 'output' argument must have a value of 'summed' or 'bySpecies'.")
    }
  } else {
    stop("'output' must be a character vector of length one.")
  }

  ## make sure that the 'dat' argument has been aggregated by genet!
  if (nrow(unique(dat[,c("Year","trackID", "Quad")])) != nrow(dat)) {
    stop("In order to be used in this function, the 'dat' argument must have
    only one row for each unique individual (genet) in each year. You can use
    the 'aggregateByGenet()' function to get a dataset with one row per genet
per year.")
  }

  # work --------------------------------------------------------------------
  ## assign a unique index to each row to simplify the looping process
  dat$index <- 1:nrow(dat)

  ## subset the 'dat' d.f so that it only has the required columns
  # (and dstore the rest)
  ## data to 'store'
  datStore <- sf::st_drop_geometry(dat[,!names(dat) %in% c("Species", "Site",
                                                           "Quad", "Year",
                                                           "trackID",
                                                           "geometry")])
  ## trimmed 'dat' to use in the function
  dat <- dat[,names(dat) %in% c("Species", "Site", "Quad", "Year",
                                "trackID", "geometry", "index")]

  ## put a buffer around each of the trackIDs in the entire data.frame
  dat <- merge(dat, buff, by = "Species")

  datBuffTemp <- sf::st_buffer(x = dat, dist = dat$buff)

  ## subtract the focal individuals from the buffered dataset
  tempBuffGeometry <- list()
  for (i in 1:nrow(dat)) {
   tempBuffGeometry[i] <- suppressWarnings(
     sf::st_difference(x = datBuffTemp[i,],y = dat[i,])$geometry)
  }
  #tempBuffGeometry <- m(FUN = sf::st_difference,
                             #x = datBuffTemp$geometry,
                             #y = dat$geometry)

  ## replace datBuffTemp geometry with the 'new' geometry ('tempBuffGeometry')
  datBuff <- sf::st_set_geometry(x = datBuffTemp,
                                 value = sf::st_as_sfc(tempBuffGeometry))

  ## trim the buffers by the boundaries of the quadrat
  datBuff <- sf::st_set_geometry(x = datBuff,
                                 value = sf::st_intersection(
                                   sf::st_as_sfc(sf::st_bbox(dat)), datBuff))

  if (method == 'count') {
    ## make an empty column in 'dat' to contain the output neighborhood data
    dat$neighbors_count <- NA
    for (i in unique(dat$Site)) { ## loop through each site
      for (j in unique(dat[dat$Site== i ,"Quad"]$Quad)) { ## loop through each
        # quadrat
        for (k in unique(dat[dat$Site == i & dat$Quad == j, "Year"]$Year)) {
          ## loop through each year
          if (compType == 'oneSpp') { ## calculating intraspecific neighborhood
            # (need to subset by species)
            for (l in unique(dat[dat$Site == i & dat$Quad == j & dat$Year == k,
                                 "Species"]$Species)) {
              ## get the data for this site/quad/year/species combo
              datOneSpp <- dat[dat$Site == i & dat$Quad == j &
                                 dat$Year == k & dat$Species == l,]
              ## get the buffered data for this site/quad/year/species combo
              datOneBuff <- datBuff[datBuff$Site == i & datBuff$Quad == j &
                                      datBuff$Year == k & datBuff$Species == l,]

              ## calculate the neighborhood for each genet
              ## get a matrix of which polygons overlap each other
              overlapM <- sf::st_intersects(datOneBuff, datOneSpp,
                                            sparse= FALSE)
              ## make the diagonal of the matrix FALSE, because a genet can't
              # overlap with itself
              diag(overlapM) <- FALSE
              ## make a list such that the list element is the row index in
              # datOneBuff of the focal indvidual, and the values in each
              # element are the row indices of the polygons in datOneSpp that
              # overlap with the focal individual
              overlapList <- apply(overlapM, MARGIN = 1, FUN = function(x)
                c(which(x==TRUE)), simplify = FALSE)
              if (length(overlapList) > 0) {
                ## get the number of genets that overlap w/ the focal buffer
                datOneSpp$neighbors <-  unlist(lapply(overlapList, length))
              } else {
                datOneSpp$neighbors <- 0
              }
              ## put the neighbor counts into the 'dat' data.frame
              dat[match(datOneSpp$index, dat$index),]$neighbors_count <-
                datOneSpp$neighbors
            }

          } else if (compType == 'allSpp') { ## calculating interspecific
            # neighborhood (don't need to subset by species)
            ## get the data for this site/quad/year combo
            datSpp<- dat[dat$Site == i & dat$Quad == j &
                           dat$Year == k ,]
            ## get the buffered data for this site/quad/year combo
            datSppBuff <- datBuff[datBuff$Site == i & datBuff$Quad == j &
                                    datBuff$Year == k ,]
            ## calculate the number of neighbors for each genet
            ## get a matrix of which polygons overlap each other
            overlapM <- sf::st_intersects(datSppBuff, datSpp, sparse= FALSE)
            ## make the diagonal of the matrix FALSE, because a genet can't
            # overlap with itself
            diag(overlapM) <- FALSE
            if (output == 'summed') {
              ## make a list such that the list element is the row index in
              # datSppBuff of the focal indvidual, and the values in each element
              # are the row indices of the polygons in datSpp that overlap with
              # the focal individual
              overlapList <- apply(X = overlapM, MARGIN = 1, FUN = function(x)
                c(which(x==TRUE)), simplify = FALSE)

              ## put the neighbor counts into the 'datSpp' data.frame
              if (length(overlapList) > 0) {
                ## get the number of genets that overlap w/ the focal buffer
                datSpp$neighbors <-  unlist(lapply(overlapList, length))
              } else {
                datSpp$neighbors <- 0
              }

              ## put the neighbor counts into the 'dat' data.frame
              dat[match(datSpp$index, dat$index),]$neighbors_count <-
                datSpp$neighbors

            } else if (output == 'bySpecies') {
              ## in overlapM, the ROW is the focal indvidiual, and the COLUMN is
              # the overlapping individual
              ## make the name of each of the columns correspond to the species
              # identity of the overlappers
              colnames(overlapM) <- datSpp$Species
              ## aggregate by column (? can you do that?) so each cell of a row
              # is the # of individuals that overlap with the focal individual
              # for a given species (given by the colname)
              # make empty matrix to hold by-species overlaps
              overlapSpp <- matrix(NA, nrow = nrow(overlapM),
                                   ncol = length(unique(colnames(overlapM))))
              colnames(overlapSpp) <- unique(colnames(overlapM))

              ## get the number of overlaps for each focal ind. (each row) for
              # each 'm' species
              for (m in unique(colnames(overlapM))) {
                # if there is more than one individual in that matrix (i.e. no
                # overlaps)
                if (nrow(overlapSpp) > 1) {
                  # is there more than one overlapping individual of species
                  # "m"? (i.e. is there a matrix of overlaps?)
                  if (is.matrix(overlapM[,colnames(overlapM) == m])) {
                    overlapSpp[,colnames(overlapSpp) == m] <-
                      rowSums(overlapM[,colnames(overlapM) == m])
                  } else {
                    # OR is there only one overlappign individual of species "m"?
                    # (i.e.is there only a vector of overlaps?)
                    overlapSpp[,colnames(overlapSpp) == m] <-
                      as.numeric(overlapM[,colnames(overlapM) == m])
                  }
                } else {# if there is only one individual in the matrix
                  overlapSpp[,colnames(overlapSpp) == m] <- sum(overlapM)
                }
              }

              overlapSpp_List <- apply(overlapSpp, MARGIN = 1,
                                       FUN = function(x) as.list(x), simplify = FALSE)

              ## put the neighbor counts into the 'datSpp' data.frame
              if (length(overlapSpp_List) > 0) {
                ## get the number of genets that overlap w/ the focal buffer
                datSpp$neighbors <- overlapSpp_List
              } else {
                datSpp$neighbors <- 0
              }

              ## put the neighbor counts into the 'dat' data.frame
              dat[match(datSpp$index, dat$index),]$neighbors_count <-
                datSpp$neighbors
            }
          }
        }
      }
      }
  } else if (method == "area") {
    ## make an empty column in 'dat' to contain the output neighborhood data
    dat$neighbors_area <- NA
    dat$nBuff_area <- NA
    ## get the overlapping polygon areas
    tempAreas <- suppressWarnings(
      sf::st_intersection(x = datBuff, y = dat))
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
                         FUN = function (x) sum(sf::st_area(x)))


      ## make sure the d.fs are in the same order (sort by index col.)
      datBuff <- datBuff[order(datBuff$index),]
      temp3 <- temp3[order(temp3$index),]

      ## put the 'area' data in the correct rows in 'dat'
      ## join by 'index' column
      dat[match(temp3$index, dat$index),]$nBuff_area <-sf::st_area(datBuff)
      dat[match(temp3$index, dat$index),]$neighbors_area <- temp3$geometry

    } else if (compType == 'allSpp') {
      ## don't need to subset by species! (so use dat and datBuff)

      ## match sites quads and years between buffered and raw data
      tempAreas2 <- tempAreas[tempAreas$Site == tempAreas$Site.1 &
                                tempAreas$Quad == tempAreas$Quad.1 &
                                tempAreas$Year == tempAreas$Year.1,]

      if (output == "summed") {
        ## now aggregate by focal genet (column name = 'trackID')
        temp3 <- aggregate(x = tempAreas2$geometry,
                           by = list("Site" = tempAreas2$Site,
                                     "Quad" = tempAreas2$Quad,
                                     "Species" = tempAreas2$Species,
                                     "trackID" = tempAreas2$trackID,
                                     "Year" = tempAreas2$Year,
                                     "index" = tempAreas2$index),
                           FUN = function (x) sum(sf::st_area(x)))
        ## make sure the two d.fs are in the same order (sort by index col.)
        datBuff <- datBuff[order(datBuff$index),]
        temp3 <- temp3[order(temp3$index),]

        ## put the 'area' data in the correct rows in 'dat'
        dat[match(temp3$index, dat$index),]$nBuff_area <-sf::st_area(datBuff)
        dat[match(temp3$index, dat$index),]$neighbors_area <- temp3$geometry
      } else if (output == "bySpecies") {
        ## now aggregate by focal genet (column name = 'trackID') AND by species
        # of the neighbors
        temp3_spp <- stats::aggregate(x = tempAreas2$geometry,
                           by = list("Site" = tempAreas2$Site,
                                     "Quad" = tempAreas2$Quad,
                                     "Species" = tempAreas2$Species,
                                     "trackID" = tempAreas2$trackID,
                                     "Year" = tempAreas2$Year,
                                     "index" = tempAreas2$index,
                                     "Species_neighbor" = tempAreas2$Species.1),
                           FUN = function (x) sum(sf::st_area(x)))

        ## put this aggregated d.f into a list format (one element for each
        # focal ind., then one next level element for each neighbor species)
        for (n in unique(temp3_spp$index)) {
          tmp <- temp3_spp[temp3_spp$index == unique(temp3_spp$index)[n],]
          tmpDF <- tmp[1,1:6]
          tmpL <- list()
          tmpL[[1]] <- as.list(tmp[,8])
          names(tmpL[[1]]) <- tmp[,7]
          tmpDF$neighbors_area <- tmpL
          if(n == unique(temp3_spp$index)[1]) {
            bySppDF <- tmpDF
          } else {
            bySppDF <- rbind(bySppDF, tmpDF)
          }
        }

        ## make sure the two d.fs are in the same order (sort by index col.)
        datBuff <- datBuff[order(datBuff$index),]
        bySppDF <- bySppDF[order(bySppDF$index),]

        ## put the 'area' data in the correct rows in 'dat'
        dat[match(bySppDF$index, dat$index),]$nBuff_area <-sf::st_area(datBuff)
        dat[match(bySppDF$index, dat$index),]$neighbors_area <-
          bySppDF$neighbors_area
      }
    }
  }

  outputDat <-  dat
  ## remove the 'buff' column
  outputDat <- outputDat[,names(outputDat)!= "buff"]

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
# dat <- grasslandData[grasslandData$Site == "AZ" &
#                     grasslandData$Quad == "SG2",]
# datIDs<- trackSpp(dat = dat, inv = grasslandInventory, dorm = 1, buff = 0.05,
#                   buffGenet = 0.005,
#                   clonal = data.frame("Species" = c("Heteropogon contortus",
#                                                     "Bouteloua rothrockii",
#                                                     "Ambrosia artemisiifolia",
#                                                     "Calliandra eriophylla" ),
#                   "clonal" = c(TRUE,TRUE,FALSE, FALSE)))
#
# names(datIDs)[c(3,4)] <- c("speciesName", "uniqueID")
#
# dataTest <- getNeighbors(dat = datIDs, buff = .15, method = "count",
#              compType = 'allSpp', output = 'summed',
#              focal = 'genet',
#              species = "speciesName",
#              trackID = "uniqueID")
#


# plot(sf::st_buffer(dataTest[dataTest$uniqueID == "AGRSMI_1997_13" &
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
#
#  dat <- grasslandData[grasslandData$Site == c("AZ") &
#   grasslandData$Species %in% c("Bouteloua rothrockii",
#   "Calliandra eriophylla"),]
#  names(dat)[1] <- "speciesName"
#
#  inv <- grasslandInventory[unique(dat$Quad)]
#  outDat <- trackSpp(dat = dat,
#   inv = inv,
#   dorm = 1,
#   buff = .05,
#   buffGenet = 0.005,
#   clonal = data.frame("Species" = unique(dat$speciesName),
#   "clonal" = c(TRUE)),
#   species = "speciesName",
#   aggByGenet = TRUE
#   )
#
#   finalDat <- getNeighbors(dat = outDat,
#   buff = .15,
#   method = 'count',
#   compType = 'oneSpp', output = "bySpecies",
#   species = "speciesName")

