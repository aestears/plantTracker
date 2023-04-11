#' Tracks genets through time
#'
#' @description This function tracks individual organisms through time, but only
#' for one species in one quadrat. It is designed for use within the
#'  \code{\link{trackSpp}} function, and is not intended for use on its own.
#'
#' @details see \code{\link{trackSpp}} for details of arguments and usage.
#'
#' @param dat An sf data.frame of the same format as
#' \code{\link{grasslandData}}. 'dat' must contain data for only one species in
#' one quadrat. It must have columns that contain...
#' * a unique identification for each research site in character format
#' with no NAs (the default column name is "Site")
#' * species name in character format with no NAs (the default column
#' name is "Species")
#' * unique quadrat identifier in character format with no NAs (the default
#' column name is "Quad")
#' *  year of data collection in integer format with no NAs (the
#' default column name is "Year")
#' * an s.f 'geometry' column that contains a= polygon or multipolygon data type
#' for each individual observation (the default column name is "geometry")
#'
#' This function will add columns called "basalArea_ramet", "trackID",
#' "age", "size_tplus1", "recruit" and "survives_tplus1", so 'dat' should not
#' contain columns with these names.
#' @param inv An integer vector that contains all of years in which this quadrat
#' (or other unique spatial area) was sampled. Years must be
#' ordered sequentially.
#' @param dorm A numeric vector of length 1, indicating the number of years this
#' species is allowed to go dormant, i.e. be absent from the map but be
#' considered the same individual when it reappears. This must be an integer
#' greater than or equal to 0.
#' @param buff A numeric vector of length 1 that is greater than or equal to
#' zero, indicating how far (in the same units as spatial values in 'dat') a
#' polygon can move from year \code{t} to year \code{t+1} and still be
#' considered the same individual.
#' @param buffGenet A numeric vector of length 1 that is greater than or equal
#' to zero, indicating how close (in the same units as spatial values in 'dat')
#' polygons must be to one another in the same year to be grouped as a genet
#' (if 'clonal' argument = TRUE). This argument is passed to
#' the \code{\link{groupByGenet}} function, which is used inside the
#' \code{\link{assign}} function.
#' @param clonal A logical argument of length 1, indicating whether this
#' species is allowed to be clonal or not (i.e. if multiple polygons (ramets)
#' can be grouped into one individual (genet)). If clonal = TRUE, the species is
#' allowed to be clonal, and if clonal = FALSE, the species is not allowed to
#' be clonal.
#' @param flagSuspects A logical argument of length 1, indicating whether
#' observations that are 'suspect' will be flagged. The default is
#' `flagSuspects = FALSE`. If `flagSuspects = TRUE`, then a column called
#' 'Suspect' is added to the output data.frame. Any suspect observations get a
#' 'TRUE' in the 'Suspect' column, while non-suspect observations receive a
#' 'FALSE'. There are two ways that an observation can be classified as
#' 'suspect'. First, if two consecutive observations have the same trackID, but
#' the basal area of the observation in year \code{t+1} is less that a certain
#' percentage (defined by the `shrink` arg.) of the basal area of the
#' observation in year t, it is possible that the observation in year \code{t+1}
#' is a new recruit and not the same individual. The second way an observation
#' can be classified as 'suspect' is if it is very small before going dormant.
#' It is unlikely that a very small individual will survive dormancy, so it is
#' possible that the function has mistakenly given a survival value of '1' to
#' this individual. A 'very small individual' is any observation with an area
#' below a certain percentile (specified by 'dormSize') of the size distribution
#' for this species, which is generated using all of the size data for this
#' species in 'dat'.
#' @param shrink A single numeric value. This value is only used when
#' `flagSuspects = TRUE`. When two consecutive observations have the same
#' trackID, and the ratio of size_t+1 to size_t is smaller than the value of
#' `shrink`, the observation in year \code{t} gets a 'TRUE' in the 'Suspect'
#' column. For example, `shrink = 0.2`, and an individual that the tracking
#' function has identified as 'BOUGRA_1992_5' has an area of 9 cm^2 in
#' year \code{t} and an area of 1.35 cm^2 in year \code{t+1}. The ratio of
#' size \code{t+1} to size \code{t} is 1.35/9 = 0.15, which is smaller than the
#' cutoff specified by `shrink`, so the observation of BOUGRA_1992_5' in year
#' \code{t} gets a 'TRUE' in the 'Suspect' column. The default value
#' is `shrink = 0.10`.
#' @param dormSize A single numeric value. This value is only used when
#' `flagSuspects = TRUE` and `dorm` is greater than or equal to 1. An individual
#' is flagged as 'suspect'
#' if it 'goes dormant' and has a size that is less than or equal to the
#' percentile of the size distribution for this species that is designated by
#' `dormSize`. For example, `dormSize = 0.05`, and the focal individual has a
#' basal area of 0.5 cm^2. The 5th percentile of the distribution of size for
#' this species, which is made using the mean and standard deviation of all
#' observations in 'dat' for the species in question, is 0.6 cm^2. This
#' individual does not have any overlaps in the next year (year \code{t+1}), but
#' does have an overlap in year \code{t+2}. However, because the basal area of
#' this observation is smaller than the 5th percentile of size for this species,
#' the observation in year \code{t} will get a 'TRUE' in the 'Suspect' column.
#' It is possible that the tracking function has mistakenly assigned a '1' for
#' survival in year t, because it is unlikely that this individual is large
#' enough to survive dormancy. The default value is `dormSize = .05`.
#' @param ... Other arguments passed on to methods. Not currently used.
#' @param  inheritsFromTrackSpp A logical argument that is only applicable when
#' 'assign()' is used internally in 'trackSpp()'.
#' @param nearEdgeBox An sf data frame indicating the bounding box of the
#' quadrat. This argument is only used if 'assign()' is used internally
#' in 'trackSpp()'.
#'
#' @seealso [trackSpp()], which is a wrapper for the [assign()] function that
#' applies it over many species and quadrats. The [assign()] function uses the
#' [groupByGenet()] function to group ramets into genets
#' (if 'clonal' argument = TRUE).
#'
#' @import sf
#' @importFrom stats aggregate reshape qnorm sd
#' @importFrom utils stack

 assign <- function(dat,
                    inv,
                    dorm,
                    buff,
                    buffGenet,
                    clonal,
                    flagSuspects = FALSE,
                    shrink = 0.1,
                    dormSize = .05,
                    inheritsFromTrackSpp = FALSE,
                    nearEdgeBox = NULL,
                    ...){
  ## argument checking --------------------------------------------------------
   ## don't really do arg. checking, since this function is not exported, and
   # expects to recieve an input directly from the 'trackSpp' function, which
   # has already done comprehensive arg. checking.

  ## internal functions ------------------------------------------------------
  ## FUNCTION FOR AGGREGATING BY GENET (if clonal or not clonal)
  ## this fxn is internal to the 'assign' fxn
  ifClonal <- function(cloneDat, clonal, buffGenet, ...) {
    ## arguments
    ## cloneDat = cloneDataset for one species/quad/year (input from 'assign()'
    # fxn), is a  subset of the 'dat' argument from the 'assign()' fxn
    ## clonal = inherits from the 'clonal' argument in the 'assign()' fxn
    ## buffGenet = inherits from the 'buffGenet' argument in the 'assign()' fxn,
    # is input into the plantTracker::groupByGenet() fxn
    if(clonal == TRUE) {
      cloneDat$genetID <- groupByGenet(cloneDat, buffGenet)
      ## aggregate size by genetID (total size for each ramet)
      tempCloneDat <- stats::aggregate(basalArea_ramet ~ genetID, sum,
                                       data = cloneDat)
      names(tempCloneDat) <- c("tempGenetID", "basalArea_genet")
      ## add aggregated size data to the cloneData df
      cloneDat$basalArea_genet <- tempCloneDat[match(cloneDat$genetID,
                                               tempCloneDat$tempGenetID),
                                         "basalArea_genet"]
    }
    ## assign unique genetIDs for every polygon (if clonal = FALSE)
    else { #if(clonal==FALSE) {
      cloneDat$genetID <- 1:nrow(cloneDat)
      cloneDat$basalArea_genet <- cloneDat$basalArea_ramet
    }

    return(cloneDat)
  }

  ## work -------------------------------------------------------------------
  ## remove the out data.frame
  if(exists("assignOut") == TRUE){
    suppressWarnings(rm("assignOut"))
  }
  ## make sure dat is in the correct sf format
  dat <- sf::st_as_sf(x = dat, sf_column_name = "geometry")
  ## add columns to the 'dat' dataset needed for output from assign()
  dat$trackID <- as.character(NA)
  dat$age <- as.integer(NA)
  dat$size_tplus1 <- as.numeric(NA)
  dat$recruit <- as.integer(NA)
  dat$survives_tplus1 <- as.integer(NA)
  dat$ghost <- as.integer(NA)
  dat$basalArea_genet <- as.numeric(NA)
  dat$genetID <- as.character(NA)
  ## if not already present, create a column called 'basalArea_ramet' in 'dat'
  if (sum(dat$basalArea_ramet) == 0) {
    dat$basalArea_ramet <-  sf::st_area(dat)
  }
  ## if 'flagSuspects' is TRUE, add a column for the flags to go into
  if (flagSuspects == TRUE) {
    dat$Suspect <- FALSE ## default value is 0
    }

  ## get the 6-letter species code for each observation
  ## make a column in the d.f with the 6-letter species code for each row
  ## only if the species column contains species name with a space
  if (sum(grepl(pattern = "[[:space:]]",x = dat$Species)) > 0  ## does the
      # species column contain a space? (is the length of the rows that contain
      # a space greater than 0?)
      ) {
    dat$sp_code_6  <- sapply(strsplit(dat$Species, " "), function(x)
      paste0(substr(toupper(x[1]), 1, 3), ## species name
             substr(toupper(x[2]), 1, 3)) ## genus name
    )
  } else if (sum(grepl(pattern = "_",x = dat$Species)) > 0) { ## if the species
    # name doesn't contain a space, does it have an underscore?
    dat$sp_code_6  <- sapply(strsplit(dat$Species, "_"), function(x)
      paste0(substr(toupper(x[1]), 1, 3), ## species name
             substr(toupper(x[2]), 1, 3)) ## genus name
    )
  } else {  ## if the species column contains some sort of code (i.e.
    # uppercase letters), or something else...
    # put the name alone in the sp_code_6 column
    dat$sp_code_6 <- dat$Species
  }


  ## assign an arbitrary, unique index number to each row in the dataset
  dat$index <- c(1:nrow(dat))

  ## find the first year in the dataset that was actually sampled
  firstDatYear <- min(dat$Year)
  ## find the index of the first year in the quadrat inventory
  firstYearIndex <- which(inv == firstDatYear)

  ## get the dataset for the first year of actually sampling
  tempPreviousYear <- dat[dat$Year == firstDatYear,]
  ## assign genetIDs to the first year of sampling dataset
  tempPreviousYear <- ifClonal(cloneDat = tempPreviousYear, clonal = clonal,
                              buffGenet = buffGenet)
  ## assign a unique trackID to every unique genetID in this first-year dataset
  IDs <- data.frame( "genetID" = sort(unique(tempPreviousYear$genetID)), ## get
                     # a vector of all of the unique genetIDs
                     "trackID" = paste0(unique(dat$sp_code_6), ## get the unique
                                        # 6-letter species code
                                        "_",unique(tempPreviousYear$Year), ## get
                                        # the unique year
                                        "_",c(1:length(unique(
                                          tempPreviousYear$genetID))))) ## get a
  # vector of  unique numbers that is the same length as the genetIDs in this
  # quad/year

  ## get the row index numbers of the rows in 'IDs' that match the genetID of
  # the row in 'tempPreviousYear'
  tempPreviousYear$trackID <- IDs[match(tempPreviousYear$genetID, IDs$genetID),
                                 "trackID"]

  ## give all individuals in year #1 a '0' in the ghost column
  tempPreviousYear$ghost <- 0

  ## check that first year in quadrat inventory and first year of
  # actual data match
  if (min(dat$Year) > inv[1]) { ## if the year of 'tempPreviousYear' is NOT
    # the same as the first year of data in the quadrat inventory (if it is
    # larger), then each individual in tempPreviousYear gets a '1' in the
    # 'recruit' column, and a '0' in the 'age' column
    tempPreviousYear$recruit <- 1
    tempPreviousYear$age <- 0
    }
  if (min(dat$Year) == inv[1]) { ## if the year of 'tempPreviousYear' is the SAME
    # as the first year of the quadrat inventory, then make sure that there is
    # an 'NA' in both the 'recruit' and 'age' columns
    tempPreviousYear$recruit <- NA
    tempPreviousYear$age <- NA
  }
  if (min(dat$Year) < inv[1]) { ## if the year of 'tempPreviousYear' is SMALLER
    # (earlier) than the first year of the quadrat inventory, then return an
    # error
    stop("Quadrat inventory dataset does not match years in the sample dataset")
  }

  ## 'tempPreviousYear' data.frame will get redefined for
  # each iteration of the for-loop below
  tempPrevI_empty <- FALSE
  ##  i = year in inventory
  if (inv[firstYearIndex] < max(inv)) { ## is the first year of data also the
    # last year of sampling?
    for (i in (firstYearIndex+1):length(inv)) {
       if (tempPrevI_empty == TRUE) {
         tempPreviousYear <- dat[dat$Year == inv[i-1], ]
       if (nrow(tempPreviousYear) > 0) {
         tempPreviousYear <- ifClonal(cloneDat = tempPreviousYear,
                                      clonal = clonal, buffGenet = buffGenet)
         tempPreviousYear$trackID <- paste0(tempPreviousYear$sp_code_6,
                                            "_", tempPreviousYear$Year, "_",
                                            tempPreviousYear$genetID)
         tempPrevI_empty <- FALSE
         } else {
         tempPrevI_empty <- TRUE
         next
       }
    }
      ## CHECK IF YEARS ARE CONTINUOUS -- check to see if the sampling years of
      # 'tempPreviousYear' and 'tempCurrentYear' are not far enough apart to
      # exceed the 'dormancy' argument. If dormancy is not exceeded, then go
      # ahead with this loop. If it is, then freshly redefine 'tempPreviousYear'
      # and proceed to the next 'i'
      if (inv[i] - inv[i-1] > (dorm+1)) { ## if the gap between years EXCEEDS
        # the dormancy argument
        ## put the 'tempPreviousYear' data into the 'assignOut' d.f, making sure
        # that the columns are in the correct order, and that the
        # 'survives_tplus1' and 'size_tplus1' columns contain 'NA'
        if (exists("assignOut") == TRUE) {
          ## if this is not the first year, then add demographic data to the
          # output d.f
          assignOut <- rbind(assignOut, tempPreviousYear)
        } else{ ## if the assignOut df is empty
          assignOut <- tempPreviousYear
        }
       ## is there any more data? if not, then end the loop
       if(inv[i] > max(dat$Year)) {
        break
       }
        ## get data from year i and put in 'tempPreviousYear'
        tempPreviousYear <- dat[dat$Year==inv[i],]

        ## determine if the tempPreviousYear data.frame has any data
        if (nrow(tempPreviousYear) < 1) { ## if there is NOT data in year i,
          # then go to the next i--but first indicate that the previous i
         # did not have any data
         tempPrevI_empty <- TRUE
          next
        }
        if (nrow(tempPreviousYear) > 0 ) { ## if there IS data in year i, then
          # group by genets and assign trackIDs, but first check if this is the
          # first year (if so, then don't have to make trackIDs b/c it already
          # has them!)

          ## is there data for genetIDs in the 'tempPreviousYear' data?
          if (sum(is.na(tempPreviousYear$genetID) == TRUE) >=
            1) { ## if there is
            # NOT data for genetID
            ## group by genet
            tempPreviousYear <- ifClonal(cloneDat  = tempPreviousYear,
                                        clonal = clonal, buffGenet = buffGenet)

            ## assign a unique trackID to every unique genetID
            IDs <- data.frame( "genetID" = sort(unique(
              tempPreviousYear$genetID)), ## get a vector of all the genetIDs
              "trackID" = paste0(unique(dat$sp_code_6), ## get the sp. code
                                 "_",unique(tempPreviousYear$Year), ## get the
                                 # unique year
                                 "_",c(1:length(unique(
                                   tempPreviousYear$genetID))))) ## get a
            # vector of unique numbers that is the same length as the genetIDs
            # in this quad/year

            ## add trackIDs to the tempPreviousYear data.frame
            tempPreviousYear$trackID <- IDs[match(tempPreviousYear$genetID,
                                                 IDs$genetID),"trackID"]

          } ## end of 'if' that determines if there is genetID data, and if not,
          # assigns genetID and trackID
        } ## end of 'if' that determines if there is data in year i
        ## end of 'if' that determines what to do if the gap between years
        # exceeds the dormancy argument
    } else { ## if the gap between years (inv[i] - inv[i-1] <= (dorm+1))
        # does NOT exceed the dormancy argument
        ## 'tempPreviousYear' is the sf data.frame of the 'current' year
        ## need to get the sf data.frame of the 'next' year (year 'i')
        tempCurrentYear <-  sf::st_as_sf(dat[dat$Year==inv[i],])

        ## MAKE SURE THERE IS DATA IN YEAR i-1 (tempPreviousYear isn't empty) If
        # not, then go to the next i (only after replacing 'tempPreviousYear'
        # with the tempCurrentYear data.frame)
        if (nrow(tempPreviousYear) < 1) { ## what to do if 'tempPreviousYear'
          # does NOT exist, then overwrite the tempPreviousYear w/ data from the
          # 'next' year, and go to the next i
          ## but first, give them trackIDs and recruit/age data (if
          # appropriate)-- because they go directly into 'tempPreviousYear,' so
          # don't get these data later in the 'orphans' section
          if (nrow(tempCurrentYear) > 0) {
            ## get a d.f that contains the obs. w/ no trackIDs
            tempTrackIDs <- tempCurrentYear[
              is.na(tempCurrentYear$trackID)==TRUE,]
            ## assign them genetIDs first
            tempTrackIDs <- ifClonal(tempTrackIDs,
                                     clonal = clonal, buffGenet = buffGenet)
             ## add genet basal areas to the tempCurrentYear data.frame
             tempCurrentYear[is.na(tempCurrentYear$trackID) == TRUE,
                             "basalArea_genet"] <- tempTrackIDs$basalArea_genet
            ## add trackIDs to the tempCurrentYear data.frame
             tempCurrentYear[is.na(
               tempCurrentYear$trackID)==TRUE,"trackID"] <- paste0(
               tempTrackIDs$sp_code_6, "_", tempTrackIDs$Year, "_",
               tempTrackIDs$genetID)
            ## then need to add 'age' and 'recruit' data (but first check that
            # this isn't the first year after a gap in sampling)
            tempCurrentYear[(tempCurrentYear$Year - inv[i-1]) < 2,
                         c("age")] <- 0
            tempCurrentYear[(tempCurrentYear$Year - inv[i-1]) < 2,
                         c("recruit")] <- 1
          }
          ## if else to determine what to do next, but depending on whether this
          # is the last year of sampling or not
          if (inv[i] == max(inv)) { ## if this IS the last year
            ## put the 'tempCurrentYear' data into the assignOut d.f
            if (exists("assignOut") == TRUE) {
              ## if this is not the first year, then add demographic data to the
              # output d.f
              assignOut <- rbind(assignOut, tempCurrentYear)
            } else{ ## if the assignOut df is empty
              assignOut <- tempCurrentYear
            }
          } else { ## if this is NOT the last year
            ## put tempCurrentYear data into tempPreviousYear, then go to the next i
            tempPreviousYear <- sf::st_as_sf(tempCurrentYear)
            next
          }
          ## end of 'if' of what to do if tempPreviousYear is empty
        } else  { ## what to do if the 'tempPreviousYear'
          # (nrow(tempPreviousYear)>=1) DOES have data
          ## add a buffer to the previous year data
          tempPreviousBuff <- sf::st_buffer(tempPreviousYear, buff)

          ## MAKE SURE THERE IS DATA IN YEAR i (tempCurrentYear)
          if (nrow(tempCurrentYear) < 1) { ## if the tempCurrentYear data does NOT
            # exist, then keep the tempPreviousYear data.frame for the next i
            ## if the gap between year of measurement (year inv[i-1]) and the
            # next i (year inv[i+1]) does not exceed the dormancy argument, then
            # roll those individuals over into 'tempPreviousYear' for the next i
            # (with a '1' in the 'ghost' column)

            ## do this workflow for years when that are not the 'last' year
            if (inv[i] != max(inv)) {
              ## get 'ghosts' (if there are any)
              ghosts <- tempPreviousYear[((inv[i+1]-tempPreviousYear$Year) <=
                                           (dorm + 1)),]
              ## put a '1' in the 'ghost' column for ghosts
              if (nrow(ghosts)>0) {
                ghosts$ghost <- '1'
              }

              ## get 'deadGhosts' (if there are any) (if the gap between year
              # i-1 and year i+1 exceeds the dormancy argument) and make sure
              # columns are in correct order
              deadGhosts <- tempPreviousYear[((inv[i+1]-tempPreviousYear$Year) >
                                               (dorm + 1)),
                                            names(dat)]

              ## rewrite 'tempPreviousYear' so it just has 'ghosts'
              # (not 'deadGhosts')
              tempPreviousYear <- ghosts

              if (nrow(deadGhosts)>0) {
                ## deadGhosts get a '0' in the 'survival' column
                deadGhosts$survives_tplus1 <- 0
                ## add the 'deadGhosts' to the 'assignOut' df
                if (exists("assignOut") == TRUE) {
                  ## if this is not the first year, then add demographic data
                  # to the output d.f
                  assignOut <- rbind(assignOut, deadGhosts)
                } else  { ## if the assignOut df is empty
                  assignOut <- deadGhosts
                }
              }
          } else {  ## what to do if this i is the 'last' year!
              ## get 'ghosts' (if there are any)
              ghosts <- tempPreviousYear[(((inv[i]+1)-tempPreviousYear$Year) <=
                                           (dorm + 1)),]

              ## get 'deadGhosts' (if there are any) (if the gap between year
              # i-1 and year i+1 exceeds the dormancy argument) and make sure
              # columns are in correct order
              deadGhosts <- tempPreviousYear[(((inv[i]+1)-tempPreviousYear$Year) >
                                               (dorm + 1)),
                                            c(names(dat))]

              if (nrow(deadGhosts)>0) {
                ## deadGhosts get a '0' in the 'survival' column
                deadGhosts$survives_tplus1 <- 0
                ## add the 'deadGhosts' to the 'assignOut' df
                if (exists("assignOut") == TRUE) {
                  ## if this is not the first year, then add demographic data
                  # to the output d.f
                  assignOut <- rbind(assignOut, deadGhosts)
                } else  { ## if the assignOut df is empty
                  assignOut <- deadGhosts
                }
              }
              ## put a '1' in the 'ghost' column for ghosts
              if (nrow(ghosts)>0) {
                ghosts$ghost <- '1'
                ghosts$survives_tplus1 <- NA
                ghosts$size_tplus1 <- NA
                ## add the 'deadGhosts' to the 'assignOut' df
                if (exists("assignOut") == TRUE) {
                  ## if this is not the first year, then add demographic data
                  # to the output d.f
                  assignOut <- rbind(assignOut, ghosts)
                } else  { ## if the assignOut df is empty
                  assignOut <- ghosts
                }
              }
            }
            ## go to next i
            next
            ## end of 'else' that has steps if tempCurrentYear is empty
          } else { ## if the tempCurrentYear data DOES exist
            ## FIND OVERLAPPING POLYGONS
            ## get the amount of overlap between each polygon
            overlapArea <- suppressWarnings(sf::st_intersection(tempPreviousBuff,
                                                                tempCurrentYear))

            ## determine if there are overlaps between the previous year and
            # next
            ## if there is NOT overlap, then skip the 'while' loop below, and
            # proceed with the code to assign trackIDs to new 'orphans' and to
            # assign 'ghost' status
            if (nrow(overlapArea) > 0) { ## if there IS overlap
              ## get the trackID names merged with the index value (same values
              # used for names of rows in 'overlaps' matrix)
              overlapArea$parentName <- paste0(overlapArea$trackID, "__",
                                               overlapArea$index)
              ## get the genetID names merged with the index value (same values
              # used for names of columns in 'overlaps' matrix)
              overlapArea$childName <- paste0("childPoly__",
                                              overlapArea$index.1)

              ## calculate the overlap between each parent poly and each child
              overlapArea$overlappingArea <- sf::st_area(overlapArea$geometry)
             # remove units from the overlap matrix, if they exist
              base::units(overlapArea$overlappingArea) <- NULL

              overlapArea <- sf::st_drop_geometry(overlapArea)

              ## get only the necessary data from the overlapArea d.f
              overlaps <- overlapArea[,c("parentName", "childName",
                                         "overlappingArea")]

              ## transform the overlaps dataframe into a 'wide' dataframe, so
              # that every row is a unique parent trackID, and the columns are
              # children (not completely aggregated yet). The value is the
              # overlap between each parent/child pair
              overlapsTemp <- stats::reshape(overlaps,
                                             v.names = "overlappingArea",
                                             idvar = "parentName",
                                             timevar = "childName",
                                             direction = "wide")

              ## remove the 'index' ID from the parent name
              overlapsTemp$parentName <- sapply(strsplit(
                overlapsTemp$parentName, "__"), unlist)[1,]

              ## aggregate by parent genet, so that each genet has only one
              # row of data (if clonal == TRUE)
              if (clonal == TRUE) {
                oldNames <- names(overlapsTemp)
                overlapsTemp <- stats::aggregate(x = overlapsTemp[,2:ncol(overlapsTemp)],
                                          by = list(
                                            "parentName" =
                                              overlapsTemp$parentName),
                                          FUN = sum, na.rm = TRUE)
                # change '0's to 'NA's
                overlapsTemp[which(overlapsTemp==0, arr.ind = TRUE)] <- NA
                names(overlapsTemp) <- oldNames
              }

              ## add parentTrackID as rownames
              rownames(overlapsTemp) <- as.character(overlapsTemp$parentName)
              ## add childGenetID as the column names
              names(overlapsTemp)  <- c("parentName", sapply(strsplit(names(overlapsTemp)[2:ncol(
                overlapsTemp)], "overlappingArea."), unlist)[2,])
              ## remove old "parentTrackID" column
              overlaps <- as.data.frame(overlapsTemp[,2:ncol(overlapsTemp)])
              rownames(overlaps) <- rownames(overlapsTemp)
              names(overlaps) <- names(overlapsTemp)[2:ncol(overlapsTemp)]

              ## determine if clonality is allowed...
              if (clonal == TRUE) { ## if yes, then each parent can have
                # multiple children (but each child can only have one parent)

                ## each child can have only one parent, so if there is only ever
                # one value in each column, than the next step is easy... each
                # column gets the trackID of the 'parent' that it overlaps with
                multParents <- apply(X = overlaps,
                                     MARGIN = 2,
                                     FUN = function(x) sum(is.na(x)==FALSE))
                # does each child only have one parent?
                if (sum(multParents > 1) != 0) {  ## if no (at least one child has
                  # more than one parent)...
                  ## get a vector of individuals with 'ties' (>1 parent)
                  ties <- names(multParents[multParents > 1])
                  if (length(ties) > 1) {
                    for (m in 1:ncol(as.data.frame(overlaps[,ties]))) {
                     # get the maximum value of overlaps
                      winnerVal <- max(overlaps[,ties][,m], na.rm = TRUE)
                      ## get the location of the highest value between the two ties
                      winner <- which(overlaps[,ties][,m] == winnerVal)
                      # if there is only one 'winner':
                      if (length(winner) ==1) {
                        ## set all other values that aren't the highest in that
                        # column to 'NA'
                        overlaps[,ties][,m][overlaps[,ties][,m] != winnerVal] <- NA
                      } else if (length(winner) > 1) {
                        ##  if there is more than 1 winner (i.e. if there are two
                        # parents with the exact same overlap with the child), then
                        # we have to break the tie some other way. Use the distance
                        # between the centroids to do this

                        ## get the name of the problem 'child'
                        badChild_name <- names(overlaps[,ties])[m]
                        ## get the spatial data for that child
                        badChild <- suppressWarnings(
                          sf::st_centroid(
                            tempCurrentYear[tempCurrentYear$index == as.numeric(
                              strsplit(badChild_name, "__")[[1]][2]),]))

                        ## get the names of the problem 'parents'
                        badParents_names <- rownames(
                          overlaps[is.na(overlaps[,ties][,m])==FALSE,])
                        badParents <-
                          tempPreviousYear[tempPreviousYear$trackID
                                           %in% badParents_names,]
                        ## aggregate badParents by trackID so that each genet has
                        # only one row of data
                        badParents <- suppressWarnings(
                          sf::st_centroid(
                            stats::aggregate(x = badParents[,"trackID"],
                                              by = list("trackID" = badParents$trackID),
                                              FUN = nrow,
                                              do_union = TRUE)))

                        ## compare the centroid distances between parents and child
                        dists <- sf::st_distance(badChild, badParents,
                                                 which = 'Euclidean')
                        rownames(dists) <- badChild_name
                        colnames(dists) <- badParents_names

                        smallDist <- colnames(dists)[which(dists == min(dists))]
                        ## replace the non-winning cell in the 'overlaps' matrix
                        # with an 'NA'
                        overlaps[rownames(overlaps) != smallDist, badChild_name] <- NA
                      }
                    } # end of loop going through each tie
                  } else if (length(ties) == 1) {
                    ## get the highest value between the two ties
                    winner <- max(overlaps[,ties], na.rm = TRUE)
                    # if there is only one 'winner':
                    if (length(which(overlaps == winner)) ==1) {
                      ## set all other values that aren't the highest in that
                      # column to 'NA'
                      overlaps[,ties][overlaps[,ties] != winner] <- NA
                    } else if (length(which(overlaps == winner)) > 1) {
                      ##  if there is more than 1 winner (i.e. if there are two
                      # parents with the exact same overlap with the child), then
                      # we have to break the tie some other way. Use the distance
                      # between the centroids to do this

                      ## get the name of the problem 'child'
                      badChild_name <- ties
                      ## get the spatial data for that child
                      badChild <- suppressWarnings(
                        sf::st_centroid(
                          tempCurrentYear[tempCurrentYear$index == as.numeric(
                            strsplit(badChild_name, "__")[[1]][2]),]))

                      ## get the names of the problem 'parents'
                      badParents_names <- rownames(overlaps)[!is.na(
                       overlaps[,ties])]
                      # get the spatial data for those parents
                      badParents <-
                        tempPreviousYear[tempPreviousYear$trackID
                                         %in% badParents_names,]
                      ## aggregate badParents by trackID so that each genet has
                      # only one row of data
                      badParents <- suppressWarnings(
                        sf::st_centroid(
                          stats::aggregate(x = badParents[,"trackID"],
                                            by = list("trackID" = badParents$trackID),
                                            FUN = nrow,
                                            do_union = TRUE)))

                      ## compare the centroid distances between parents and child
                      dists <- sf::st_distance(badChild, badParents,
                                               which = 'Euclidean')
                      rownames(dists) <- badChild_name
                      colnames(dists) <- badParents_names

                      smallDist <- colnames(dists)[which(dists == min(dists))]
                      ## replace the non-winning cell in the 'overlaps' matrix
                      # with an 'NA'
                      overlaps[rownames(overlaps) != smallDist, badChild_name] <- NA
                    }
                  }
                }

                ## after the previous loop, there is now only one parent for each
                # child. Proceed w/ assigning appropriate trackIDs to children

                ## get the numeric genetIDs of the children (for each parent)
                nameDF <- utils::stack(apply(X = t(overlaps), MARGIN = 1,
                                      function(x) names(x[which(is.na(x)==FALSE)])))
                # rename the columns so they make sense
                names(nameDF) <- c("parentTrackID", "childIndex")

                ## get just the index for the 'children'
                nameDF$childIndex <- as.numeric(sapply(
                  strsplit(
                    as.character(nameDF$childIndex), split = "__"), unlist)[2,])

                ## put the trackID from the parent into the tempCurrentYear d.f
                # rows that correspond to the genetID of the child
                tempCurrentYear[match(nameDF$childIndex, tempCurrentYear$index),
                                "trackID"] <- nameDF$parentTrackID
              } else if (clonal == FALSE) { ## if no, then each parent can have
                # only one child and each child can have only one parent

                ## set up a 'while' loop
                 whileOverlaps <- overlaps
                 done <- FALSE
                 counter <- 0

                 while (!done) {

                   ## get the row and column indices of the maximum value
                   maxInds <- which(whileOverlaps ==
                                      max(whileOverlaps,na.rm = TRUE),
                                    arr.ind = TRUE)
                   ## what if there is a tie?? --i.e. there are >=2 overlaps
                   # that have the same value
                   if (nrow(maxInds) > 1) {
                     # get the names of the potential parents
                     maybeParents_names <- rownames(whileOverlaps)[unique(
                       maxInds[,"row"])]
                     # get the names of the potential children
                     maybeChildren_names <- names(whileOverlaps)[unique(maxInds[,"col"])]
                     ## get a matrix of centroid distances
                     # the the spatial data for maybeChidlren
                     maybeChildren <- suppressWarnings(
                       sf::st_centroid(
                         tempCurrentYear[tempCurrentYear$index %in%
                                           as.numeric(
                                             as.vector(
                                               data.frame(
                                                 strsplit(
                                                   x = maybeChildren_names,
                                                   split = "__")
                                                 )[2,])),]))

                     # get the spatial data for maybeParents
                     maybeParents <-
                       tempPreviousYear[tempPreviousYear$trackID
                                        %in% maybeParents_names,]

                     ## aggregate badParents by trackID so that each genet has
                     # only one row of data
                     maybeParents <- suppressWarnings(
                       sf::st_centroid(
                         stats::aggregate(x = maybeParents[,"trackID"],
                                          by = list("trackID" = maybeParents$trackID),
                                          FUN = nrow,
                                          do_union = TRUE)))

                     ## compare the centroid distances between parents and child
                     dists <- sf::st_distance(maybeChildren, maybeParents,
                                              which = 'Euclidean')
                     rownames(dists) <- maybeChildren_names
                     colnames(dists) <- maybeParents_names

                     ## get the indices of the smallest distance pair
                     smallDistInds <- which(dists == min(dists), arr.ind = TRUE)

                  # if the smallDistInds d.f contains data for TWO parents (>1 column number!)
                  if (length(unique(smallDistInds[,"col"])) > 1) {
                    for (n in unique(smallDistInds[,"col"])) {
                      # get the data just for the nth parent
                      tinyDistInds <- smallDistInds[smallDistInds[,"col"]==n,]
                      # get the name of the nth child
                      tinyChild <- rownames(dists)[tinyDistInds["row"]]
                      # get the name of the nth parent
                      tinyParent <- colnames(dists)[tinyDistInds["col"]]
                      # give the nth child the trackID of the nth parent (in tempCurrentYear)
                      tempCurrentYear[tempCurrentYear$index ==
                                        strsplit(x = tinyChild, split = "__")[[1]][2],
                                      "trackID"] <- tinyParent
                      # update the "whileOverlaps" matrix w/ appropriate zeros
                      whileOverlaps[tinyParent, tinyChild] <- 0
                      whileOverlaps[tinyParent, ] <- 0
                      whileOverlaps[, tinyChild] <- 0
                    }
                  } else {
                    # if the smallDistInds d.f contains data for only one parent (same column number)

                     ## get the index of the smallest distance child
                     smallChild <- rownames(dists)[smallDistInds[,"row"]]
                     ## get the trackID of the smallest distance parent
                     smallParent <- colnames(dists)[smallDistInds[,"col"]]

                     ## put the trackID from the parent into the tempCurrentYear
                     # d.f rows that correspond to the row index of the child
                     tempCurrentYear[
                       tempCurrentYear$index==strsplit(
                         x = smallChild, split = "__")[[1]][2],
                                     "trackID"] <- smallParent

                     ## overwrite the 'max' value with an NA, so we can find the
                     # next largest value
                     whileOverlaps[smallParent,smallChild] <- 0

                     ## overwrite all of the other values in the parent row and
                     # the child column with NAs also (since each parent can
                     #only have one child, and each child can only have
                     # one parent)
                     whileOverlaps[smallParent,] <- 0
                     whileOverlaps[,smallChild] <- 0
                    }

                   } else { ## if there is only one maximum overlap--no ties

                     ## get the trackID of the parent
                     maxParent <- rownames(whileOverlaps)[maxInds[1,1]]

                     ## get the numeric index of the 'child'
                     maxChild <- as.numeric(strsplit(
                       colnames(whileOverlaps)[maxInds[1,2]],"__")[[1]][2])

                     ## put the trackID from the parent into the tempCurrentYear
                     # d.f rows that correspond to the genetID of the child
                     tempCurrentYear[tempCurrentYear$index==maxChild,
                                     "trackID"] <- maxParent

                     ## overwrite the 'max' value with an NA, so we can find the
                     # next largest value
                     whileOverlaps[maxInds[1,1],maxInds[1,2]] <- 0

                     ## overwrite all of the other values in the parent row and
                     # the child column with NAs also (since each parent can only
                     # have one child, and each child can only have one parent)
                     whileOverlaps[maxInds[1,1],] <- 0
                     whileOverlaps[,maxInds[1,2]] <- 0
                   }

                   ## update the 'counter'
                   counter <- counter + 1

                   if (counter > 3000) {

                     stop("tracking 'while' loop is running out of control!")

                   } ## end of 'if' checking that the counter isn't too large

                   if (sum(whileOverlaps, na.rm = TRUE)==0) {
                     ## if the sum of the  matrix is empty (is all NAs), then
                     # stop the while loop
                     done <- TRUE
                   } ## end of 'if' that's redefining 'done' if matrix is NAs
                 } ## end of 'while' that's finding the max values in the matrix
              } # end of 'if' that determines if clonal is 'TRUE' or 'FALSE'
            } # end of 'if' that determines if there ARE overlaps

            ## make optional checks that will flag obs. as 'suspect'
            if (flagSuspects == TRUE) {
              ## check: individuals w/ the same trackID-- can't have a decrease
              # in size more than 90% and still be the same individual (i.e.
              # current year must be > 10% of the size of previous year)

              ## make sure the areas are aggregated by genet
              smallPrevious <- stats::aggregate(x = sf::st_drop_geometry(
                tempPreviousYear[, "basalArea_ramet"]),
                by = list("Year_prev" = tempPreviousYear$Year,
                          "trackID" = tempPreviousYear$trackID),
                FUN = sum
              )
              names(smallPrevious)[3] <- "Area"

              smallCurrent <- stats::aggregate(x = sf::st_drop_geometry(
                tempCurrentYear[, "basalArea_ramet"]),
                by = list("Year_curr" = tempCurrentYear$Year,
                          "trackID" = tempCurrentYear$trackID),
                FUN = sum
              )
              names(smallCurrent)[3] <- "Area"

              ## get a list of the trackIDs that are present in both the
              # previous and current years
              shrinkage <- merge(smallPrevious, smallCurrent,
                                 by = "trackID")

              ## get ratio of current year size to previous year size. Is the
              # ratio less than or equal to .1?
              if (nrow(shrinkage) > 0) {
                shrinkers <- shrinkage$trackID[(shrinkage$Area.y /
                                                  shrinkage$Area.x) <= shrink]

                if (length(shrinkers) > 0) { ## if there are any shrinkers...
                  ## flag the observation in the PREVIOUS year
                  tempPreviousYear[tempPreviousYear$trackID %in% shrinkers,
                                  "Suspect"] <- TRUE
                }
              }

              ## check: for individuals that are dormant (and only if dorm = 1),
              # then a really tiny organism can't become a really big organism
              # (i.e. an organism that is really tiny probably can't go dormant)
              if (dorm >= 1 & ## if the dorm argument is > 0...
                  length(unique(round(dat$basalArea_ramet, 5))) > 3 ## if the
                  # sizes are not all the same (i.e. if the species is  measured
                  # as polygons, not points)
                  ) {
                ## if there are data that survive from year 1 to year 2
                if (nrow(shrinkage) > 0) {
                  ## is there a difference of greater than 1 year between any of
                  # the individuals with the same trackID?
                  ## get individuals that have a gap > 1
                  # between current and previous year
                  dormants <- shrinkage[shrinkage$trackID %in%
                                          which((shrinkage$Year_curr -
                                                   shrinkage$Year_prev ) > 1),]
                  if (nrow(dormants) > 0) {
                    ## get individuals that have a 'very small size' in the
                    # previous year (as long as the size isn't exactly the same
                    # from year to year--since they are likely then points and
                    # have a fixed radius)
                    dormants <- dormants[round(dormants$Area.x,8) !=
                                           round(dormants$Area.y,8),]
                    ## get the trackID of individuals in the previous year that
                    # were smaller than the 5th percentile of the distribution
                    # of sizes for this species
                    smallest <- exp(qnorm(p = dormSize,
                                          mean = mean(log(dat$Area)),
                                          sd = sd(log(dat$Area))))
                    tooSmallIDs <- dormants[dormants$Area.x < smallest,
                                            "trackID"]

                    ## flag individuals in the PREVIOUS year
                    # that are likely too small to have survived dormancy
                    tempPreviousYear[tempPreviousYear$trackID %in% tooSmallIDs,
                                    "Suspect"] <- TRUE
                  }
                }
              }
            }

            ## ORPHANS: deal with 'child' polygons that don't have parents
            ## give them new trackIDs
            orphans <- tempCurrentYear[is.na(tempCurrentYear$trackID)==TRUE,]
            ## make sure that there are orphans
            if (nrow(orphans)>0) {
              ## orphans are NOT grouped by genet, since it is very unlikely
              # that new recruits consist of multiple ramets

              ## add the orphan trackIDs to the 'orphan's data.frame
              orphans$trackID <- paste0(orphans$sp_code_6, "_",orphans$Year, "_", 1:nrow(orphans))

              ## populate the 'basalArea_genet' column for the orphans
              orphans$basalArea_genet <- orphans$basalArea_ramet

              ## check that the orphans don't come after a gap in sampling (if
              # they do, then leave NA's for all demographic values, if they
              # don't, then proceed with the following code)
              ## check that this is NOT the first year of sampling (i=1)
              if (inv[i] - inv[i-1] <= 1) { ## if this year does NOT come after
                # a gap in sampling
                ## assign a '1' in the recruits column
                orphans$recruit <- 1
                orphans$age <- 0
              }
            } ## end of if that determines if there are any orphans

            ## CHILDREN
            ## (first have to assign the parents)
            ## define the PARENTS data.frame
            parents <- tempPreviousYear[tempPreviousYear$trackID %in%
                                         tempCurrentYear$trackID,]
            ## define the children data.frame
            children <- tempCurrentYear[tempCurrentYear$trackID %in%
                                       tempPreviousYear$trackID,]

            ## ASSIGN DEMOGRAPHIC DATA TO CHILDREN
            if (nrow(children)>0) {
              ## give children a 0 in the recruit column, since they have a
              # parent (as long as the parent doesn't have an NA--was recruited
              # in a year when we couldn't know when it was recruited)
              if (inv[i] - inv[i-1] <= 1) {
                children$recruit <- parents[match(children$trackID,
                                                  parents$trackID),]$recruit
                ## if the parent is an 'NA', the child also gets an 'NA', but if
                # the parent is a '1', the child gets a '0'. If the parent is a
                # '0', then the child gets a '0'
                if (nrow(children[children$recruit==1 &
                                  is.na(children$recruit) == FALSE,]) > 0) {
                  children[children$recruit==1 &
                             is.na(children$recruit) == FALSE,]$recruit <- 0
                }
              } else {
                children$recruit <- NA
              }
              ## give the children the appropriate age column (only if the
              # parents don't have an 'NA' for age) (parent's age + 1)
              ## get the trackID and age+1 of the parents in the 'tempParents'
              # df
              ## get the two relevant columns and get rid of geometry column
              tempParents <- sf::st_drop_geometry(parents[,c("trackID", "age")])
              ## add 1 to the parent age (get an 'NA' if the parent age is NA)
              tempParents$age <- (tempParents$age + 1)
              names(tempParents) <- c("trackIDtemp", "age")

              ## add the'age'column from the appropriate parents (+1), joined by
              # trackID
              children$age <- tempParents[match(children$trackID,
                                                tempParents$trackIDtemp),]$age

              ## calculate the appropriate 'basalArea_genet' data for each child
              childGenet_area <- stats::aggregate(children$basalArea_ramet,
                                           by = list(
                                             "trackID" = children$trackID), sum)
              names(childGenet_area) <- c("trackID", "basalArea_genet")
              children <- merge(x = children[,names(children) !="basalArea_genet"],
                    y = childGenet_area, by = "trackID")
            }

            ## PARENTS
            ## ASSIGN DEMOGRAPHIC DATA FOR PARENTS
            ## make sure that there is actually a parents data.frame
            if (nrow(parents) > 0) {
              ## PARENTS ('parents' data.frame + 'deadGhosts' data.frame)
              ## give the 'parents' a '1' for survival
              parents[,"survives_tplus1"] <- 1
              ## give the 'parents' a 0 in the 'ghosts' column
              parents[,"ghost"] <- 0
              ## assign size in the current year (From 'children' df) to the
              # appropriate rows in the parents data.frame
              childSizeTemp <- unique(sf::st_drop_geometry(
                children[,c("trackID", "basalArea_genet")]))
              names(childSizeTemp) <- c("trackIDtemp", "size_tplus1")
              ## join size_tplus data to 'parents' data.frame
              parents$size_tplus1 <- childSizeTemp[
                match(parents$trackID,childSizeTemp$trackIDtemp),]$size_tplus1
            }

            ## ONLY ALLOW GHOSTS IF DORMANCY > 1
            if (dorm > 0) {
              ## GHOSTS: parent polygons that don't have 'children' in year i.
              # If they were observed in a year that is more than dorm+1 years
              # prior to year i, then they don't get saved to the current year. If
              # they were observed in a year that is >= to dorm+1 years prior to
              # year i, then they get added to the 'tempCurrentYear' data.frame,
              # which both go into the 'tempPreviousYear' data.frame before the
              # next iteration of the loop
              ## get the ghost individuals
              ghostsTemp <- tempPreviousYear[!(tempPreviousYear$trackID %in%
                                                tempCurrentYear$trackID),]

              ## is 'i' the last year of sampling? If not, then do the
              # following:
              if (inv[i] != max(inv)) {
                ## check that these individuals can be 'ghosts' in the current
                # year (if the gap between the year of their observation and
                # year i+1 is greater than the dormancy argument (+1), then)
                # they are not ghosts, and get a 0 for survival
                ghosts <- ghostsTemp[((inv[i+1] - ghostsTemp$Year) <=
                                        (dorm + 1)),]
                ## put the 'ghosts' that exceed the dormancy argument into their
                # own data.frame
                deadGhosts <- ghostsTemp[((inv[i+1] - ghostsTemp$Year) >
                                            (dorm + 1)),]
                ## if the 'i' year is the last year of sampling:
              } else { #if (inv[i] == max(inv)) {
                ## get the 'ghost' data (have to use a different 'i+1' value if
                # this is the current year)
                ghosts <- ghostsTemp[inv[i]+1 - ghostsTemp$Year <= (dorm+1),]
                ## get the 'deadGhost' data
                deadGhosts <- ghostsTemp[inv[i]+1 - ghostsTemp$Year > (dorm+1),]
              }

              ## ASSIGN DEMOGRAPHIC DATA TO GHOSTS
              ##give these survived ghosts a '1' in the 'ghost' column
              if (nrow(ghosts)>0) {
                ghosts$ghost <- 1
                ghosts$survives_tplus1 <- NA
              }
              ## give the deadGhosts a 0 in survival (if there are any
              # deadGhosts)
              if (nrow(deadGhosts)>0) {
                deadGhosts$survives_tplus1 <- 0
                deadGhosts$ghost <- 0
              }
              ## get the data.frames in a consistent format
              ghosts <- ghosts[,names(children)]
              deadGhosts <- deadGhosts[,names(children)]
            } else { ## if (dorm==0); any trackID that doesn't have a child in
              # year 'i' is 'dead'
              deadGhosts <- tempPreviousYear[!(tempPreviousYear$trackID %in%
                                                tempCurrentYear$trackID),]
              if (nrow(deadGhosts) > 0) { ## make sure there are deadGhosts
                ## give the dead ghosts a '0' for survival
                deadGhosts$survives_tplus1 <- 0
              }
              deadGhosts <- deadGhosts[,names(children)]
              ghosts <- NULL
            }

            ## PREPARE FOR NEXT i
            ## arrange columns of children, orphans, and ghosts in the same
            # order
            orphans <- orphans[, names(children)]
            parents <- parents[,names(children)]

            ## bind children, orphans, and ghosts into one data.frame, that will
            # become the data for the previous year in the next i of the loop
            tempCurrentYear <- rbind(children, orphans, ghosts)

            if (exists("assignOut") == TRUE) {
              ## if this is not the first year, then add demographic data to the
              # output d.f
              assignOut <- rbind(assignOut, parents, deadGhosts)
            } else { #if (exists("assignOut") == FALSE) {
              ## if the assignOut df is empty
              assignOut <- suppressWarnings(rbind(parents, deadGhosts))
            }

            ## if this is the last year of sampling, put the 'tempCurrentYear'
            # data into the 'assignOut' d.f, but only after making sure that
            # they have an 'NA' in the 'survives_tplus1' and
            # 'size_tplus1' columns
            if (inv[i] == max(inv)) {
              assignOut <- rbind(assignOut, tempCurrentYear)
            }

            ## assign the data from the current i (tempCurrentYear) to be the
            # past year data (tempPreviousYear) for the next i
            tempPreviousYear <- tempCurrentYear
          } ## end of 'if' statement that determines if the tempCurrentYear data
          # exists
          } ## end of 'if' that determines if the tempPreviousYear data exists
        } ## end of 'if' statement that determines if gap between inv[i-1] and
      # inv[i] is less than or equal to  dorm+1
      } ## end of loop i
    } else if (inv[firstYearIndex] == max(inv)) {
    ## if there are only observations in the last year of 'inv', then there
    # cannot be any demographic data assigned.
    ## make sure there are 'NA's' in the appropriate columns
    tempPreviousYear[, c('size_tplus1', 'survives_tplus1')] <- NA
    assignOut <- tempPreviousYear
    }

  ## make an empty column
  assignOut$nearEdge <- FALSE
  ## populate the 'nearEdge' column (if isn't inheriting data from trackSpp)
  if (inheritsFromTrackSpp == FALSE) {
    ## make a boundary box that is within the 'buff' argument of the actual quad
    buffEdgeOutside <- sf::st_as_sfc(sf::st_bbox(assignOut))
    buffEdgeInside <- sf::st_as_sfc(sf::st_bbox(assignOut) + c(buff, buff,-buff, -buff))
    buffEdge <-  sf::st_difference(buffEdgeOutside, buffEdgeInside)
    } else if (inheritsFromTrackSpp == TRUE) {
    buffEdge <- nearEdgeBox
  }

  ## find out which polys intersect with the buffered quad
  assignOut[sf::st_intersects(assignOut, buffEdge, sparse = FALSE),
            "nearEdge"] <- TRUE

  ## clean up output assignOut data.frame (remove NAs and unneeded columns)
  assignOut <- assignOut[is.na(assignOut$Species)==FALSE,
                         !(names(assignOut) %in% c("ghost","genetID", "index",
                                                   "sp_code_6"))]
  ## output ---------------------------------------------------------------
return(assignOut)
  }

# testing -----------------------------------------------------------------
# example input data ------------------------------------------------------
# #
# ## prepares the dataset to feed into the 'assign' function (the 'trackSpp'
# # function will do this ahead of time when the user calls it)
# sampleDat <- grasslandData[grasslandData$Site == "AZ"
#                            & grasslandData$Quad == "SG4"
#                            & grasslandData$Species == "Bouteloua rothrockii",]
# # this should be a data.frame
# dat <- sampleDat
# #
# # # get the appropriate grasslandInventory data for the "unun_11" quadrat,
# # # to tell the 'assign' function when the quadrat was sampled
# sampleInv<- grasslandInventory[["SG2"]]
# # this should be an integer vector
# inv <- sampleInv
#
# # remove a year to simulate 'dorm' being exceeded
# dat <- sampleDat[sampleDat$Year != 1924,]
# inv <- inv[c(1:2,4:5)]
#
# testOutput <- assign(dat = dat,
#                      inv = inv,
#                      dorm = 1,
#                      buff = .05,
#                      # buffGenet = .001,
#                      clonal =  FALSE,
#                      flagSuspects = TRUE)

#
# ggplot(testOutput) +
#   geom_sf(aes(fill = as.factor(trackID),colour = Year), alpha = 0.5) +
#   scale_fill_discrete(guide = "none") +
#   theme_classic()
# #
# # ggplot() +
# #   geom_sf(data = st_buffer(testOutput[testOutput$Year==2009,],.05), fill = "purple", alpha = 0.5) +
# #   geom_sf(data = st_buffer(dat[dat$Year==2009,],.05), fill = "green", alpha = 0.5) +
# #   theme_linedraw()
# #
# #
# # ggplot() +
# #   geom_sf(data = dat[dat$Year==1923,], fill = "green", alpha = 0.5) +
# #   geom_sf(data = testOutput[testOutput$Year==1923,], fill = "purple", alpha = 0.5) +
# #   theme_linedraw()
# #
# # tempDat <- dat[dat$Year==1923,]
# # tempOut <- testOutput[testOutput$Year==1923,]
# # mapview(tempDat) + mapview(tempOut, col.regions = "green")
#
# ## find duplicated polygons in testOutput d.f
# #testOutput_cent <- st_centroid(testOutput)
#
# dupValues <- testOutput[duplicated(testOutput$geometry),"geometry"]$geometry
#
# dups <- testOutput[testOutput$geometry %in% dupValues,]
# plot(dups[dups$geometry==dupValues[1],"geometry"])
#
# mapview(dups[dups$geometry==dupValues[1],])
# ## all the duplicates are in 2003, what?? that doesn't make sense
# unique(dups$Year)
#
# #see which obs. aren't in the outPut dataset (in comparison to the raw data)
# testDat <- st_drop_geometry(dat)
# testDat$test <- "old"
# testOutputTest <- st_drop_geometry(testOutput)
# testOutputTest$test <- "new"
#
# testTest <- full_join(testDat,testOutputTest, by = c("Species", "Clone", "Seedling", "Stems", "Basal", "Type", "Site", "Quad", "Year", "sp_code_4", "sp_code_6", "Area"))
#
