#' tracks individual plants through time
#'
#' @return
#' @export
#'
#' @examples

# required packages -------------------------------------------------------
library(sf)
library(mapview) #don't actually need for function, just for checking
library(dplyr)

# example input data ------------------------------------------------------
# grasslandData (or exact same format), subset to a unique site, quad,
# species
load("./data/grasslandData.rda")
load("./data/grasslandInventory.rda")
# prepares the dataset to feed into the 'assign' function (the 'Assign'
# function will do this ahead of time when the user calls it)
sampleDat <- grasslandData[grasslandData$Site == "CO"
                           & grasslandData$Quad == "unun_11"
                           & grasslandData$Species == "Bouteloua gracilis",]
# this should be a data.frame

# get the appropriate grasslandInventory data for the "unun_11" quadrat,
# to tell the 'assign' function when the quadrat was sampled
sampleInv <- grasslandInventory[["unun_11"]]
# this should be an integer vector

# 'assign' function ---------------------------------------------------------

assign <- function(sampleDat, inv, dorm, buff, buffGenet, overlap, clonal,...){
  # arguments ---------------------------------------------------------------
  ## double check the format of the inputs, and add additional columns required
  dat <- st_sf(sampleDat) # data in 'grasslandData' format, must be an sf object
    ## add columns for trackID, age, size_t+1, and recruit
    dat$trackID <- NA
    dat$age <- NA
    dat$size_tplus1 <- NA
    dat$recruit <- NA
    dat$survives_tplus1 <- NA
    dat$index <- c(1:nrow(dat))## assign an arbitrary, unique index number to
    # each row in the dataset
    dat$ghost = NA
    dat$trackID = NA
  inv <- sort(sampleInv) ## integer vector of quadrat sampling years in
  # sequential order

  ## user-defined arguments
  dorm <- 1 ## dormancy allowed by the function
  buff <- .05 ## buffer of movement allowed from year to year, in meters
  buffGenet <- 0.001 ## buffer between polygons (/2) that is the maximum allowed for
  # them to be considered genets
  overlap <- .50 ## the percentage of overlap (in decimal form) between focal
  # indiv. and next year indiv. that will be required to consider them both
  # the same individual (not 100% sure on this one yet...)
  clonal <- 1 ## binary option that indicates whether this species is allowed to
  # break into multiple polygons for one individual

  #conalBuff <- .01 ## the buffer of overlap that indicates polygons are genets
  ## of the same individual (if genets are allowed), in meters ##probably don't
  ##need...

  ## work -------------------------------------------------------------------

  datFirst <- dat[dat$Year==inv[1],] ## get the data just for the first year of
  # sampling
  if (clonal==1) { ## if this species is 'clonal' (clonal argument == 1), use
    # Dave's groupByGenet function to give a singleTrackID to a group of
    # polygons if they are close enough to one another (user-defined)
    datFirst$genetID <- groupByGenet(datFirst, buffGenet) # put these temporary genet
    # ID's in a temporary 'genetID' column, which will tell us to give these
    # genets the same track ID

  }
  if (clonal == 0) { ## if genets are not allowed (clonal == 0), then each
    # individual polygon gets a unique value in the 'genetID' column
    datFirst$genetID <- 1:nrow(datFirst)
  }

  ## assign a unique trackID to every unique genetID in this first-year dataset
  IDs <- data.frame( "genetID" = sort(unique(datFirst$genetID)), ## get a vector
                     # of all of the unique genetIDs
            "trackID" = paste0(unique(dat$sp_code_6), ## get the unique 6-letter
                               # species code
                               "_",unique(datFirst$Year), ## get the unique year
                  "_",c(1:length(unique(datFirst$genetID))))) ## get a vector of
  # unique numbers that is the same length as the genetIDs in this quad/year

  datFirst<- merge(datFirst[,names(datFirst) != "trackID"], IDs, by = "genetID")

  ## assign the first-year data to the 'tempCurrentYear' data.frame
  tempCurrentYear <- datFirst ## this data.frame will get redefined for
  # each iteration of the for-loop below

  ## give all individuals in year #1 a '0' in the ghost column
  tempCurrentYear$ghost <- 0

  ##  i = year in inventory
  for (i in 2:length(inv)) {

    ## CHECK IF YEARS ARE CONTINUOUS -- check to see if the sampling years of
    # 'tempCurrentYear' and 'tempNextYear' are not far enough apart to exceed
    # the 'dormancy' argument. If dormancy is not exceeded, then go ahead with
    # this loop. If it is, then freshly redefine 'tempCurrentYear' and proceed
    # to the next 'i'
    if (inv[i] - inv[i-1] <= (dorm+1)) {
      ## 'tempCurrentYear' is the sf data.frame of the 'current' year
      tempCurrentBuff <- st_buffer(tempCurrentYear, buff) ## need to add a buffer
      # to this data.frame

      ## need to get the sf data.frame of the 'next' year (year 'i')
      tempNextYear <- dat[dat$Year==inv[i],]

      ## AGGREGATE BY GENET for year i (if clonal = 1)
      if(clonal==1) {
        tempNextYear$genetID <- groupByGenet(tempNextYear, buffGenet)
      }
      ## assign unique genetIDs for every polygon (if clonal = 0)
      if(clonal==0) {
        tempNextYear$genetID <- 1:nrow(tempNextYear)
      }

      ## FIND OVERLAPPING POLYGONS
      ## trying to get the amount of overlap between each polygon
      overlapArea <- st_intersection(tempCurrentBuff, tempNextYear)

      ## get the trackID names merged with the index value (same values used
      # for names of rows in 'overlaps' matrix)
      overlapArea$parentName <- paste0(overlapArea$trackID, "__",
                                       overlapArea$index)
      ## get the genetID names merged with the index value (same values used for
      # names of columns in 'overlaps' matrix)
      overlapArea$childName <- paste0("genet__", overlapArea$genetID.1, "__",
                                      overlapArea$index.1)

      ## calculate the overlap between each parent poly and each child poly
      overlapArea$overlappingArea <- st_area(overlapArea$geometry)

      overlapArea <- st_set_geometry(overlapArea, NULL)

      ## aggregate the overlaps by rows (by 'parents')
      overlaps <- overlapArea %>%
        group_by(trackID, genetID.1) %>%
        summarize(overlappingArea = sum(overlappingArea))

      names(overlaps) <- c("parentTrackID", "childGenetID",
                              "overlappingArea")

      overlaps$childGenetID <- paste0("genet__", overlaps$childGenetID,
                "__",seq(from = 500, by = 1, length.out = nrow(overlapTemp)))

      ## transform the overlaps dataframe into a 'wide' dataframe, so that
      # every row is a unique parent trackID, and the columns are children
      # (not completely aggregated yet)
      overlaps <- as.data.frame(pivot_wider(overlaps,
                                      id_cols = c(parentTrackID, childGenetID),
                                            names_from = childGenetID,
                                            values_from = overlappingArea))
      rownames(overlaps) <- overlaps$parentTrackID
      overlaps <- overlaps[,2:ncol(overlaps)]

      ## transpose the 'overlaps' matrix so that we can aggregate by genet
      # (genetID are rows, trackID are columns)
      overlaps <- t(overlaps)
      ## simplify the row names so they only have the genetID number
      overlaps <- as.data.frame(overlaps)
      overlaps$genetID <- as.numeric(sapply(strsplit(rownames(overlaps),"__"),
                                            unlist)[2,])
      overlaps[is.na(overlaps)==TRUE] <- 0

      overlaps <- overlaps %>%
        group_by(genetID) %>%
        summarize(across(names(overlaps)[1:(ncol(overlaps)-1)], sum))

      overlaps<- as.data.frame(overlaps)
      rownames(overlaps) <- paste0("genet_",overlaps$genetID)
      overlaps <- overlaps[,2:ncol(overlaps)]

      ## transpose the matrix back so that each row is a parent and each
      # column is a child
      overlaps <- t(overlaps)

      ## each parent can only have one child, and each child can only have one
      # parent (since we've already clustered the polygons by genet). Each
      # parent-child pair will be determined by the greatest amount of overlap.
      # This is done by going through the 'overlaps' matrix and finding the
      # largest overlap value. Then, the trackID from that parent will be
      # assigned to the child.

      whileOverlaps <- overlaps
      done <- FALSE
      counter <- 0
      while (!done) {
        ## get the row and column indices of the maximum value
        maxInds <- which(whileOverlaps == max(whileOverlaps, na.rm = TRUE),
                        arr.ind = TRUE)
        ## get the trackID of the parent
        maxParent <- rownames(whileOverlaps)[maxInds[1,1]]
        ## get the numeric genetID of the child
        maxChild <- strsplit(colnames(whileOverlaps)[maxInds[1,2]],"_")[[1]][2]

        ## put the trackID from the parent into the tempNextYear dataframe rows
        # that correspond to the genetID of the child
        tempNextYear[tempNextYear$genetID==maxChild,"trackID"] <- maxParent

        ## overwrite the 'max' value with an NA, so we can find the next
        # largest value
        whileOverlaps[maxInds[1,1],maxInds[1,2]] <- NA
        ## overwrite all of the other values in the parent row and the child
        # column with NAs also (since each parent can only have one child, and
        # each child can only have one parent)
        whileOverlaps[maxInds[1,1],] <- NA
        whileOverlaps[,maxInds[1,2]] <- NA


        counter <- counter + 1 ## update the 'counter'
        if (counter > 500) {
          stop("tracking while loop is running out of control!")
        } ## end of 'if' thats checking that the counter isn't too large

        if (sum(whileOverlaps, na.rm = TRUE)==0) {
          ## if the sum of the 'whileOverlaps' matrix is all 0 (is all NAs),
          # then stop the while loop
        done <- TRUE
          } ## end of 'if' thats redefining 'done' if matrix is all NAs
      } ## end of 'while' thats finding the max values in the matrix

      ## ORPHANS: deal with 'child' polygons that don't have parents
      ## give them new trackIDs
      orphans <- tempNextYear[is.na(tempNextYear$trackID)==TRUE,] %>%
        select(!trackID)

      ## make a unique trackID for each unique genetID in the 'orphan' dataset
      orphanIDs <- data.frame(
        "genetID" = unique(tempNextYear[is.na(tempNextYear$trackID)==TRUE,
            "genetID"]$genetID), ## get the unique genetIDs of the 'orphans'
        "trackID" = paste0(unique(
          tempNextYear$sp_code_6), ## get the unique 6-letter species code
               "_",unique(tempNextYear$Year), ## get the unique year
               "_",c(1:length(unique(tempNextYear[is.na(
                 tempNextYear$trackID)==TRUE,"genetID"]$genetID)))))

      ## add demographic data for the orphans! (but first check that this isn't
      # the first year of sampling after a gap) ###AES###

      ## add the orphan trackIDs to the 'orphan's data.frame
      orphans <- left_join(orphans, orphanIDs, by = "genetID") %>%
        select(names(tempNextYear))
      ## overwrite the orphan rows in the 'tempNextYear' data.frame w/ the data
      # from the 'orphans' data.frame that contains trackIDs
      tempNextYear <- rbind(tempNextYear[!(tempNextYear$genetID %in%
                                             orphanIDs$genetID),],orphans)

      ## GHOSTS: parent polygons that don't have 'children' in year i. If
      # they were observed in a year that is more than dorm+1 years prior to
      # year i, then they don't get saved to the next year. If they were
      # observed in a year that is >= to dorm+1 years prior to year i, then
      # they get added to the 'tempNextYear' data.frame, which both go into the
      # 'tempCurrentYear' data.frame before the next iteration of the loop

      ## get the ghost individuals
      ghosts <- tempCurrentYear[!(tempCurrentYear$trackID %in%
                                    tempNextYear$trackID),]
      ## check that these indivduals can be 'ghosts' in the next year (if the
      # gap between the year of their observation and year i+1 is greater than
      # the dormancy argument (+1), then) they are not ghosts, and get a 0 for
      # survival
      ghosts <- ghosts[((inv[i+1] - ghosts$Year) <= (dorm + 1)),]
      ##give these ghosts a '1' in the 'ghost' column
      ghosts[,"ghost"] <- 1

      ## redefine the tempCurrentYear data.frame, but w/out the 'ghosts'
      tempCurrentYear <- tempCurrentYear[!(tempCurrentYear$index %in%
                                           ghosts$index),]

      ggplot() +
        #geom_sf(data = st_buffer(tempCurrentYear,.01), fill = "blue", alpha = 0.5) +
        geom_sf(data = tempCurrentYear, fill = "blue", alpha = 0.8) +
        geom_sf(data = ghosts, fill = "red", alpha = 0.8) +
        theme_classic()

      ###AES### need to figure out what to do next

      ## ASSIGN DEMOGRAPHIC DATA to the individuals that survived from
      # year i-1 to year i
      tempCurrentYear

    } ## end of 'if' statement that determines if gap between inv[i-1] and
    # inv[i] is less than or equal to  dorm+1

    ##PREPARE FOR NEXT i
    ## get data for year i and put it in 'tempCurrentYear' for the next 'i'
    tempCurrentYear <- tempNextYear
    ## add a 'ghost' column and put in an 'NA' that indicates these polygons
    # were actually observed in the 'current' year
    tempCurrentYear$ghost <- 0
    ## add the 'ghost' data from year i-1 (parents w/ no children)
    tempCurrentYear <- rbind(tempCurrentYear, tempGhosts)
  } #end of loop i

# output ---------------------------------------------------------------
}
