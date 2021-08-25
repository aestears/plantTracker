#' Tracks genets through time
#'
#' @description This function tracks individual plants through time, but only of
#'  one species in one quadrat. It is designed for use within the
#'  \code{\link{trackSpp}} function, and is not intended for use on its own.
#'
#' @details see \code{\link{trackSpp}} for details of arguments and usage.
#'
#' @param dat An sf data.frame of the same format as
#' \code{\link{grasslandData}}. 'dat' can contain data for only one species in
#' one quadrat. It must have columns that contain a unique
#' identification for each research site (default name is "Site"), species name
#' (default name is "Species"), quadrat identifier (default name is "Quad"),
#' year of data collection (default name is "Year") and an s.f 'geometry' column
#' that contains a polygon or multipolygon data type for each
#' individual observation.
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
#' polygon can move from year \code{i} to year \code{i}+1 and still be
#' considered the same individual.
#' @param buffGenet A numeric vector of length 1 that is greater than or equal
#' to zero, indicating how close (in the same units as spatial values in 'dat')
#' polygons must be to one another in the same year to be grouped as a genet
#' (if 'clonal' argument = TRUE). This argument is passed to
#' the \code{\link{groupByGenet}} function, which is used inside the
#' \code{\link{assign}} function.
#' @param clonal A logical argument of length 1, indicating whether this
#' species is allowed to be clonal or not (i.e. if multiple polygons (ramets)
#' can be grouped as one individual (genet)). If clonal = TRUE, the species is
#' allowed to be clonal, and if clonal = FALSE, the species is not allowed to
#' be clonal.
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @seealso [trackSpp()], which is a wrapper for the [assign()] function that
#' applies it over many species and quadrats. The [assign()] function uses the
#' [groupByGenet()] function to group ramets into genets
#' (if 'clonal' argument = TRUE).
#'
#' @import sf
#' @importFrom stats aggregate reshape

 assign <- function(dat, inv, dorm , buff , buffGenet, clonal,
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
    # fxn), is a  subset of the 'cloneDat' argument from the 'assign()' fxn
    ## clonal = inherits from the 'clonal' argument in the 'assign()' fxn
    ## buffGenet = inherits from the 'buffGenet' argument in the 'assign()' fxn,
    # is input into the PlantTracker::groupByGenet() fxn
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
    dat$basalArea_ramet <-  st_area(dat)
  }

  ## assign an arbitrary, unique index number to each row in the dataset
  dat$index <- c(1:nrow(dat))

  ## find the first year in the dataset that was actually sampled
  firstDatYear <- min(dat$Year)
  ## find the index of the first year in the quadrat inventory
  firstYearIndex <- which(inv==firstDatYear)

  ## get the dataset for the first year of actually sampling
  tempCurrentYear <- dat[dat$Year==firstDatYear,]
  ## assign genetIDs to the first year of sampling dataset
  tempCurrentYear <- ifClonal(cloneDat = tempCurrentYear, clonal = clonal,
                              buffGenet = buffGenet)
  ## assign a unique trackID to every unique genetID in this first-year dataset
  IDs <- data.frame( "genetID" = sort(unique(tempCurrentYear$genetID)), ## get
                     # a vector of all of the unique genetIDs
                     "trackID" = paste0(unique(dat$sp_code_6), ## get the unique
                                        # 6-letter species code
                                        "_",unique(tempCurrentYear$Year), ## get
                                        # the unique year
                                        "_",c(1:length(unique(
                                          tempCurrentYear$genetID))))) ## get a
  # vector of  unique numbers that is the same length as the genetIDs in this
  # quad/year

  ## get the row index numbers of the rows in 'IDs' that match the genetID of
  # the row in 'tempCurrentYear'
  tempCurrentYear$trackID <- IDs[match(tempCurrentYear$genetID, IDs$genetID),
                                 "trackID"]

  ## give all individuals in year #1 a '0' in the ghost column
  tempCurrentYear$ghost <- 0


  ## check that first year in quadrat inventory and first year of
  # actual data match
  if (min(dat$Year) > inv[1]) { ## if the year of 'tempCurrentYear' is NOT
    # the same as the first year of data in the quadrat inventory (if it is
    # larger), then each individual in tempCurrentYear gets a '1' in the
    # 'recruit' column, and a '0' in the 'age' column
    tempCurrentYear$recruit <- 1
    tempCurrentYear$age <- 0
    }
  if (min(dat$Year) == inv[1]) { ## if the year of 'tempCurrentYear' is the SAME
    # as the first year of the quadrat inventory, then make sure that there is
    # an 'NA' in both the 'recruit' and 'age' columns
    tempCurrentYear$recruit <- NA
    tempCurrentYear$age <- NA
  }
  if (min(dat$Year) < inv[1]) { ## if the year of 'tempCurrentYear' is SMALLER
    # (earlier) than the first year of the quadrat inventory, then return an
    # error
    stop("Quadrat inventory dataset does not match years in the sample dataset")
  }

  ## 'tempCurrentYear' data.frame will get redefined for
  # each iteration of the for-loop below

  ##  i = year in inventory
  if (inv[firstYearIndex] < max(inv)) {
    for (i in (firstYearIndex+1):22) { #length(inv)) {
      ## CHECK IF YEARS ARE CONTINUOUS -- check to see if the sampling years of
      # 'tempCurrentYear' and 'tempNextYear' are not far enough apart to exceed
      # the 'dormancy' argument. If dormancy is not exceeded, then go ahead with
      # this loop. If it is, then freshly redefine 'tempCurrentYear' and proceed
      # to the next 'i'
      if (inv[i] - inv[i-1] > (dorm+1)) { ## if the gap between years EXCEEDS
        # the dormancy argument
        ## put the 'tempCurrentYear' data into the 'assignOut' d.f, making sure
        # that the columns are in the correct order, and that the
        # 'survives_tplus1' and 'size_tplus1' columns contain 'NA'
        if (exists("assignOut") == TRUE) {
          ## if this is not the first year, then add demographic data to the
          # output d.f
          assignOut <- rbind(assignOut, tempCurrentYear)
        } else{ ## if the assignOut df is empty
          assignOut <- tempCurrentYear
        }

        ## get data from year i and put in 'tempCurrentYear'
        tempCurrentYear <- dat[dat$Year==inv[i],]

        ## determine if the tempCurrentYear data.frame has any data
        if (nrow(tempCurrentYear) < 1) { ## if there is NOT data in year i, then
          # go to next i (but only after overwriting the tempCurrentYear
          # data.frame with data from the next i)
          tempCurrentYear <- dat[dat$Year==inv[i+1],]
          next
        }
        if (nrow(tempCurrentYear) > 0 ) { ## if there IS data in year i, then
          # group by genets and assign trackIDs, but first check if this is the
          # first year (if so, then don't have to make trackIDs b/c it already
          # has them!)

          ## is there data for genetIDs in the 'tempCurrentYear' data?
          if(is.null(tempCurrentYear$genetID)==TRUE) { ## if there is NOT data
            # for genetID
            ## group by genet
            tempCurrentYear <- ifClonal(cloneDat  = tempCurrentYear,
                                        clonal = clonal, buffGenet = buffGenet)

            ## assign a unique trackID to every unique genetID
            IDs <- data.frame( "genetID" = sort(unique(
              tempCurrentYear$genetID)), ## get a vector of all the genetIDs
              "trackID" = paste0(unique(dat$sp_code_6), ## get the sp. code
                                 "_",unique(tempCurrentYear$Year), ## get the
                                 # unique year
                                 "_",c(1:length(unique(
                                   tempCurrentYear$genetID))))) ## get a
            # vector of unique numbers that is the same length as the genetIDs
            # in this quad/year

            ## add trackIDs to the tempCurrentYear data.frame
            tempCurrentYear$trackID <- IDs[match(tempCurrentYear$genetID,
                                                 IDs$genetID),"trackID"]

          } ## end of 'if' that determines if there is genetID data, and if not,
          # assigns genetID and trackID
        } ## end of 'if' that determines if there is data in year i
        ## end of 'if' that determines what to do if the gap between years
        # exceeds the dormancy argument
      } else { ## if the gap between years does; (inv[i] - inv[i-1] <=
        # (dorm+1)) NOT exceed the dormancy argument
        ## 'tempCurrentYear' is the sf data.frame of the 'current' year
        ## need to get the sf data.frame of the 'next' year (year 'i')
        tempNextYear <-  sf::st_as_sf(dat[dat$Year==inv[i],])

        if (nrow(tempNextYear)>0) { ## if there is data in the
          # tempNextYear d.f (and clonal arg. is true), then assign genetIDs
          tempNextYear <- ifClonal(cloneDat = tempNextYear,
                                   clonal = clonal,
                                   buffGenet = buffGenet)
        }

        ## MAKE SURE THERE IS DATA IN YEAR i-1 (tempCurrentYear isn't empty) If
        # not, then go to the next i (only after replacing 'tempCurrentYear'
        # with the tempNextYear data.frame)
        if (nrow(tempCurrentYear) < 1) { ## what to do if 'tempCurrentYear' does
          # NOT exist, then overwrite the tempCurrentYear w/ data from the
          # 'next' year, and go to the next i
          ## but first, give them trackIDs and recruit/age data (if
          # appropriate)-- because they go directly into 'tempCurrentYear,' so
          # don't get these data later in the 'orphans' section
          if (nrow(tempNextYear) > 0) {
            ## get a d.f that contains the obs. w/ no trackIDs
            tempTrackIDs <- tempNextYear[is.na(tempNextYear$trackID)==TRUE,]

            ## add trackIDs to the tempCurrentYear data.frame
            tempNextYear[is.na(tempNextYear$trackID)==TRUE,"trackID"] <- paste0(
              tempTrackIDs$sp_code_6, "_", tempTrackIDs$Year, "_",
              tempTrackIDs$genetID)

            ## then need to add 'age' and 'recruit' data (but first check that
            # this isn't he first year after a gap in sampling)
            tempNextYear[(tempNextYear$Year - inv[i-1]) < 2,
                         c("age")] <- 0
            tempNextYear[(tempNextYear$Year - inv[i-1]) < 2,
                         c("recruit")] <- 1
          }
          ## if else to determine what to do next, but depending on whether this
          # is the last year of sampling or not
          if (inv[i]==max(inv)) { ## if this IS the last year
            ## put the 'tempNextYear' data into the assignOut d.f
            if (exists("assignOut") == TRUE) {
              ## if this is not the first year, then add demographic data to the
              # output d.f
              assignOut <- rbind(assignOut, tempNextYear)
            } else{ ## if the assignOut df is empty
              assignOut <- tempNextYear
            }
          } else { ## if this is NOT the last year
            ## put tempNextYear data into tempCurrentYear, then go to the next i
            tempCurrentYear <- sf::st_as_sf(tempNextYear)
            next
          }
          ## end of 'if' of what to do if tempCurrentYear is empty
        } else  { ## what to do if the 'tempCurrentYear';
          # (nrow(tempCurrentYear)>=1)DOES have data
          ## add a buffer to the current year data
          tempCurrentBuff <- sf::st_buffer(tempCurrentYear, buff)

          ## MAKE SURE THERE IS DATA IN YEAR i (tempNextYear i
          if (nrow(tempNextYear) < 1) { ## if the tempNextYear data does NOT
            # exist, then keep the tempCurrentYear data.frame for the next i
            ## if the gap between year of measurement (year inv[i-1]) and the
            # next i (year inv[i+1]) does not exceed the dormancy argument, then
            # roll those individuals over into 'tempCurrentYear' for the next i
            # (with a '1' in the 'ghost' column)

            ## do this workflow for years when that are not the 'last' year
            if (inv[i] != max(inv)) {
              ## get 'ghosts' (if there are any)
              ghosts <- tempCurrentYear[((inv[i+1]-tempCurrentYear$Year) <=
                                           (dorm + 1)),]
              ## put a '1' in the 'ghost' column for ghosts
              if (nrow(ghosts)>0) {
                ghosts$ghost <- '1'
              }

              ## get 'deadGhosts' (if there are any) (if the gap between year
              # i-1 and year i+1 exceeds the dormancy argument) and make sure
              # columns are in correct order
              deadGhosts <- tempCurrentYear[((inv[i+1]-tempCurrentYear$Year) >
                                               (dorm + 1)),
                                            names(dat)]

              ## rewrite 'tempCurrentYear' so it just has 'ghosts'
              # (not 'deadGhosts')
              tempCurrentYear <- ghosts

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
              ghosts <- tempCurrentYear[(((inv[i]+1)-tempCurrentYear$Year) <=
                                           (dorm + 1)),]

              ## get 'deadGhosts' (if there are any) (if the gap between year
              # i-1 and year i+1 exceeds the dormancy argument) and make sure
              # columns are in correct order
              deadGhosts <- tempCurrentYear[(((inv[i]+1)-tempCurrentYear$Year) >
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
            ## end of 'else' that has steps if tempNextYear is empty
          } else { ## if the tempNextYear data DOES exist,
            # ((nrow(tempNextYear)>0)) proceed with the loop
            ## AGGREGATE BY GENET for year i (if clonal = TRUE)
            tempNextYear <- ifClonal(cloneDat = tempNextYear, clonal = clonal,
                                     buffGenet = buffGenet)

            ## FIND OVERLAPPING POLYGONS
            ## trying to get the amount of overlap between each polygon
            overlapArea <- suppressWarnings(sf::st_intersection(tempCurrentBuff,
                                                                tempNextYear))

            ## determine if there are overlaps between the current year and
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
              overlapArea$childName <- paste0("genet__",
                                              overlapArea$genetID.1,
                                              "__",
                                              overlapArea$index.1)

              ## calculate the overlap between each parent poly and each child
              overlapArea$overlappingArea <- sf::st_area(overlapArea$geometry)

              overlapArea <- sf::st_set_geometry(overlapArea, NULL)

              ## aggregate the overlaps by rows (by 'parents')
              overlaps <- stats::aggregate(
                overlapArea$overlappingArea, by = list(
                  overlapArea$trackID, overlapArea$genetID.1), FUN = sum)


              names(overlaps) <- c("parentTrackID", "childGenetID",
                                   "overlappingArea")

              overlaps$childGenetID <- paste0(
                "genet__", overlaps$childGenetID, "__",seq(
                  from = 500, by = 1, length.out = nrow(overlaps)))

              ## transform the overlaps dataframe into a 'wide' dataframe, so
              # that every row is a unique parent trackID, and the columns are
              # children (not completely aggregated yet). The value is the
              # overlap between each parent/child pair
              overlapsTemp <- stats::reshape(overlaps,
                                             v.names = "overlappingArea",
                                             idvar = "parentTrackID",
                                             timevar = "childGenetID",
                                             direction = "wide")

              ## correct the column names for this data.frame
              ## remove old "parentTrackID" column
              temp <- as.data.frame(overlapsTemp[,2:ncol(overlapsTemp)])
              ## add parentTrackID as rownames
              rownames(temp) <- as.character(overlapsTemp$parentTrackID)
              ## add childGenetID as the column names
              names(temp)  <- sapply(strsplit(names(overlapsTemp)[2:ncol(
                overlapsTemp)], "overlappingArea."), unlist)[2,]

              overlaps <- temp

              ## transpose the 'overlaps' matrix so that we can aggregate by
              # genet (genetID are rows, trackID are columns)
              overlaps <- t(overlaps)
              ## simplify the row names so they only have the genetID number
              overlaps <- as.data.frame(overlaps)
              overlaps$genetID <- as.numeric(sapply(strsplit(rownames(overlaps),
                                                             "__"),
                                                    unlist)[2,])
              overlaps[is.na(overlaps)==TRUE] <- 0

              ## aggregate the overlap data by genetID
              overlapsTemp <- stats::aggregate(
                overlaps[,1:(ncol(overlaps)-1)], by = list(overlaps$genetID),
                FUN = sum)

              temp <- as.data.frame(overlapsTemp[,2:ncol(overlapsTemp)])
              ## fix rownames (childGenetID)
              rownames(temp) <- paste0("genet_",overlapsTemp$Group.1)
              ## fix colnames (parentTrackID)
              names(temp) <- names(overlaps)[names(overlaps)!="genetID"]

              overlaps <- temp

              ## transpose the matrix back so that each row is a parent and each
              # column is a child
              overlaps <- t(overlaps)

              ## each parent can only have one child, and each child can only
              # have one parent (since we've already clustered the polygons by
              # genet). Each parent-child pair will be determined by the
              # greatest amount of overlap. This is done by going through the
              # 'overlaps' matrix and finding the largest overlap value. Then,
              # the trackID from that parent will be assigned to the child.

              whileOverlaps <- overlaps
              done <- FALSE
              counter <- 0
              while (!done) {
                ## get the row and column indices of the maximum value
                maxInds <- which(whileOverlaps == max(whileOverlaps,
                                                      na.rm = TRUE),
                                 arr.ind = TRUE)
                ## get the trackID of the parent
                maxParent <- rownames(whileOverlaps)[maxInds[1,1]]
                ## get the numeric genetID of the child
                maxChild <- strsplit(colnames(
                  whileOverlaps)[maxInds[1,2]],"_")[[1]][2]

                ## put the trackID from the parent into the tempNextYear d.f
                # rows that correspond to the genetID of the child
                tempNextYear[tempNextYear$genetID==maxChild,
                             "trackID"] <- maxParent

                ## overwrite the 'max' value with an NA, so we can find the next
                # largest value
                whileOverlaps[maxInds[1,1],maxInds[1,2]] <- NA
                ## overwrite all of the other values in the parent row and the
                # child column with NAs also (since each parent can only have
                # one child, and each child can only have one parent)
                whileOverlaps[maxInds[1,1],] <- NA
                whileOverlaps[,maxInds[1,2]] <- NA


                counter <- counter + 1 ## update the 'counter'
                if (counter > 500) {
                  stop("tracking 'while' loop is running out of control!")
                } ## end of 'if' thats checking that the counter isn't too large

                if (sum(whileOverlaps, na.rm = TRUE)==0) {
                  ## if the sum of the  matrix is empty (is all NAs), then stop
                  # the while loop
                  done <- TRUE
                } ## end of 'if' that's redefining 'done' if matrix is all NAs
              } ## end of 'while' that's finding the max values in the matrix
              ## end of 'if' that contains the work if there ARE overlaps
            }

            ## ORPHANS: deal with 'child' polygons that don't have parents
            ## give them new trackIDs
            orphans <- tempNextYear[is.na(tempNextYear$trackID)==TRUE,]
            ## make sure that there are orphans
            if (nrow(orphans)>0) {
              ## make a unique trackID for each genetID in the 'orphan' dataset
              orphanIDs <- data.frame(
                "genetID" = unique(tempNextYear[is.na(
                  tempNextYear$trackID)==TRUE, "genetID"]$genetID), ## get the
                # unique genetIDs of the 'orphans'
                "trackID" = paste0(unique(
                  tempNextYear$sp_code_6), ## get the unique 6-letter species
                  # code
                  "_",unique(tempNextYear$Year), ## get the unique year
                  "_",c(1:length(unique(tempNextYear[is.na(
                    tempNextYear$trackID)==TRUE,"genetID"]$genetID)))))

              ## add the orphan trackIDs to the 'orphan's data.frame
              orphans$trackID <- orphanIDs[match(orphans$genetID,
                                                 orphanIDs$genetID),"trackID"]

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
            parents <- tempCurrentYear[tempCurrentYear$trackID %in%
                                         tempNextYear$trackID,]
            ## define the children data.frame
            children <- tempNextYear[tempNextYear$trackID %in%
                                       tempCurrentYear$trackID,]

            ## ASSIGN DEMOGRAPHIC DATA TO CHILDREN
            if (nrow(children)>0) {
              ## give children a 0 in the recruit column, since they have a
              # parent (as long as the parent doesn't have an NA--was recruited
              # in a year when we couldn't know when it was recruited)
              if (inv[i] - inv[i-1] <= 1) {
                children$recruit <- 0
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
              ## assign size in the next year (From 'children' df) to the
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
              # prior to year i, then they don't get saved to the next year. If
              # they were observed in a year that is >= to dorm+1 years prior to
              # year i, then they get added to the 'tempNextYear' data.frame,
              # which both go into the 'tempCurrentYear' data.frame before the
              # next iteration of the loop
              ## get the ghost individuals
              ghostsTemp <- tempCurrentYear[!(tempCurrentYear$trackID %in%
                                                tempNextYear$trackID),]

              ## is 'i' the last year of sampling? If not, then do the
              # following:
              if (inv[i] != max(inv)) {
                ## check that these individuals can be 'ghosts' in the next year
                # (if the gap between the year of their observation and year i+1
                # is greater than the dormancy argument (+1), then) they are not
                # ghosts, and get a 0 for survival
                ghosts <- ghostsTemp[((inv[i+1] - ghostsTemp$Year) <=
                                        (dorm + 1)),]
                ## put the 'ghosts' that exceed the dormancy argument into their
                # own data.frame
                deadGhosts <- ghostsTemp[((inv[i+1] - ghostsTemp$Year) >
                                            (dorm + 1)),]
                ## if the 'i' year is the last year of sampling:
              } else { #if (inv[i] == max(inv)) {
                ## get the 'ghost' data (have to use a different 'i+1' value if
                # this is the next year)
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
              if (nrow(deadGhosts)>0){
                deadGhosts$survives_tplus1 <- 0
                deadGhosts$ghost <- 0
              }
              ## get the data.frames in a consistent format
              ghosts <- ghosts[,names(children)]
              deadGhosts <- deadGhosts[,names(children)]
            } else { ## if (dorm==0); any trackID that doesn't have a child in
              # year 'i' is 'dead'
              deadGhosts <- tempCurrentYear[!(tempCurrentYear$trackID %in%
                                                tempNextYear$trackID),]
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
            # become the data for the current year in the next i of the loop
            tempNextYear <- rbind(children, orphans, ghosts)

            if (exists("assignOut") == TRUE) {
              ## if this is not the first year, then add demographic data to the
              # output d.f
              assignOut <- rbind(assignOut, parents, deadGhosts)
            } else { #if (exists("assignOut") == FALSE) {
              ## if the assignOut df is empty
              assignOut <- rbind(parents, deadGhosts)
            }

            ## if this is the last year of sampling, put the 'tempNextYear' data
            # into the 'assignOut' d.f, but only after making sure that they
            # have an 'NA' in the 'survives_tplus1' and 'size_tplus1' columns
            if (inv[i] == max(inv)) {
              assignOut <- rbind(assignOut, tempNextYear)
            }

            ## assign the data from the current i (tempNextYear) to be the past
            # year data (tempCurrentYear) for the next i
            tempCurrentYear <- tempNextYear
          } ## end of 'if' statement that determines if the tempNextYear data
          # exists
        } ## end of 'if' that determines if the tempCurrentYear data exists
      } ## end of 'if' statement that determines if gap between inv[i-1] and
      # inv[i] is less than or equal to  dorm+1
    } ## end of loop i
  } else {
    ## if there are only observations in the last year of 'inv', then there
    # cannot be any demographic data assigned.
    ## make sure there are 'NA's' in the appropriate columns
    tempCurrentYear[, c('size_tplus1', 'survives_tplus1')] <- NA
    assignOut <- tempCurrentYear
  }
  ## populate the 'nearEdge' column
  ## make an empty column
  assignOut$nearEdge <- FALSE
  ## make a boundary box that is within the 'buff' argument of the actual quad
  buffEdgeOutside <- st_as_sfc(st_bbox(assignOut))
  buffEdgeInside <- st_as_sfc(st_bbox(assignOut) + c(buff, buff, -buff,
                                                     -buff))
  buffEdge <-  st_difference(buffEdgeOutside, buffEdgeInside)
  ## find out which polys intersect with the buffered quad
  assignOut[st_intersects(assignOut, buffEdge, sparse = FALSE),
            "nearEdge"] <- TRUE

  ## clean up output assignOuta.frame (remove NAs and unneeded columns)
  assignOut <- assignOut[is.na(assignOut$Species)==FALSE,
                         !(names(assignOut) %in% c("ghost","genetID", "index",
                                                   "sp_code_6"))]
  ## output ---------------------------------------------------------------
return(assignOut)
  }

# testing -----------------------------------------------------------------
 # example input data ------------------------------------------------------
#
# prepares the dataset to feed into the 'assign' function (the 'Assign'
# function will do this ahead of time when the user calls it)
# sampleDat <- grasslandData[grasslandData$Site == "KS"
#                            & grasslandData$Quad == "q33"
#                            & grasslandData$Species == "Ambrosia psilostachya",]
# # this should be a data.frame
# dat <- sampleDat
# #
# # # get the appropriate grasslandInventory data for the "unun_11" quadrat,
# # # to tell the 'assign' function when the quadrat was sampled
# sampleInv<- grasslandInventory[["q33"]]
# # this should be an integer vector
# inv <- sampleInv
#
# testOutput <- assign(dat = dat,
#                      inv = inv,
#                      dorm = 1,
#                      buff = .05,
#                      buffGenet = .001,
#                      clonal =  1)


# # ggplot(testOutput) +
# #   geom_sf(aes(fill = as.factor(trackID)), alpha = 0.5) +
# #   scale_fill_discrete(guide = FALSE) +
# #   theme_classic()
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
