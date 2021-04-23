#' tracks individual plants through time
#'
#' @return
#' @export
#'
#' @examples

# required packages -------------------------------------------------------
library(sf)
library(mapview) #don't actually need for function, just for checking

# example input data ------------------------------------------------------
# grasslandData (or exact same format), subset to a unique site, quad,
# species
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

assign <- function(sampleDat, inv, dorm, buff, overlap, clonal,...){
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
  inv <- sort(sampleInv) ## integer vector of quadrat sampling years in
  # sequential order

  ## user-defined arguments
  dorm <- 1 ## dormancy allowed by the function
  buff <- .05 ## buffer of movement allowed from year to year, in meters
  buffGenet <- 0 ## buffer between polygons (/2) that is the maximum allowed for
  # them to be considered genets
  overlap <- .50 ## the percentage of overlap (in decimal form) between focal
  # indiv. and next year indiv. that will be required to consider them both
  # the same individual (not 100% sure on this one yet...)
  clonal <- 1 ## binary option that indicates whether this species is allowed to
  # break into multiple polygons for one individual
  conalBuff <- .01 ## the buffer of overlap that indicates polygons are genets
  # of the same individual (if genets are allowed), in meters

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
  IDs <- data.frame( "genetID" = sort(unique(datFirst$genetID)), ## get a vector of unique numbers
                       "trackID" = paste0(unique(dat$sp_code_6),"_",unique(datFirst$Year),
                                          "_",c(1:length(unique(datFirst$genetID)))))
  datFirst<- merge(datFirst[,names(datFirst) != "trackID"], IDs, by = "genetID")

  ## put this trackID information in the master data.frame for this species
  dat[dat$index %in% datFirst$index, "trackID"] <- datFirst$trackID
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

      ## need to get the sf data.frame of the 'next' year
      tempNextYear <- dat[dat$Year==inv[i],]

      ## AGGREGATE BY GENET for year i (if clonal = 1)
      if(clonal==1) {
        tempNextYear$genetID <- groupByGenet(tempNextYear, buffGenet)
      }

      ## FIND OVERLAPPING POLYGONS
      ## see if there is any overlap between the tempNextYear data and the
      # tempCurrentYear data (buffered)
      overlaps <- st_intersects(tempNextYear, tempCurrentBuff)
      # returns a list, where each object in the list is the row index of the
      # polygon in tempNextYear, and the contents of the list are row indices of
      # the overlapping polygons from tempCurrentBuff. Ultimately, we must chose

      ###AES### need to deal with the fact that now we can have multiple parent
      #polygons!! (still can only have one parent trackID-wise)

      ## UNAMBIGUOUS PARENT: assign trackIDs to polygons from year t+1 that #
      # overlap only one polygon from year t
      ## get the row index, index , and trackID numbers of the tempNextYear
      # (child) polygons that only have one parent, and teh tempCurrentYear
      # (parent) polygons that only have one child
      oneParentRowNums <- which(sapply(overlaps, length)==1)
      unambigParent <- data.frame("childRowIndex" = oneParentRowNums,
                                  "childIndex" = tempNextYear[oneParentRowNums,"index"]$index,
                                  "parentTrackID" = sapply(overlaps[oneParentRowNums], function(x)
                                    tempCurrentYear[unlist(x),"trackID"]$trackID
                                  ),
                                  "parentRowIndex" = sapply(overlaps[oneParentRowNums], unlist
                                  ),
                                  "parentIndex" = sapply(overlaps[oneParentRowNums], function(x)
                                    tempCurrentYear[unlist(x),"index"]$index))

      mapview(tempCurrentBuff, col.regions = "pink") + mapview(tempNextYear, col.regions = "yellow") + mapview(tempNextYear[tempNextYear$index %in% unambigParent$childIndex,], col.regions = "orange")

      ## put the trackIDs from the parents into the appropriate spot in the
      # big data.frame for each child
      dat[dat$index %in% unambigParent$childIndex, "trackID"] <-
        unambigParent$parentTrackID

      ## AMBIGUOUS PARENT / TIEBREAKERS: are there any polygons from tempNextYear
      # that overlap with
      # more than one polygon from tempCurrentYear?
      ## get the row index numbers of the tempCurrentYear (parent) polygons that
      # overlap with the same t+1 polygons (child) -- the parents that share a
      # child polygon, which isn't allowed (two different parents can't have
      # the same child)
      currentOverlapRowNums <-
        sort(unique(unlist(
          overlaps[which(
            sapply(overlaps, length)  ## add up how many year t polys overlap each
            # year t+1 poly
            >1)## identify the year t+1 polys that are overlapped by >1 year t poly
            ] ## get the elements of the list that correspond to each of the year
          # t+1 polys w/ >1 year t poly overlap
        ) ## get a numeric vector of the row index for each of those multiple
        # overlap year t polys
        ))

      currentYearOverlaps <- tempCurrentBuff[currentOverlapRowNums,]

      ## get the row index numbers of the tempNextYear (child) polygons that are
      # overlapped by  more than one polygon from year t (parent) -- children
      # than have >1 parent (not allowed)
      nextOverlapRowNums <- which(sapply(overlaps, length) >1)
      nextYearOverlaps <- tempNextYear[nextOverlapRowNums,]

      ## compare the overlap between current year and next year for each combo of
      # polygons
      for(j in nextOverlapRowNums) {
        ## get the data for the first polygon in year t+1 that overlaps with >1
        # polygon from year t
        overlapNext <- tempNextYear[j,]
        ## get the data for the polygons from year t that overlap this one
        overlapCurrents <- tempCurrentBuff[overlaps[[j]],]

        ## compare the extent of overlap between year t+1 poly and each of the
        # year t polys
        tempOverlap <- st_intersection(overlapNext, overlapCurrents)
        ## determine which of the year t polys overlaps the year t poly the most
        overlapWinner <- tempOverlap[tempOverlap$Area.1 ==
                                       max(tempOverlap$Area.1),]
        ## put the trackID of the 'overlapWinner' poly in the complete dataset
        #row for the year t+1 polygon. Match by using the 'index' (unique number
        # assigned to each row in the raw dataset)
        dat[dat$index==overlapNext$index,"trackID"] <- overlapWinner$trackID.1
      } #end of loop j

      ## NO PARENT (overlaps list index will not have any information in it)
      ## get the row index of child polygons that are 'new' (i.e. no parent)
      noParentRowNums <- which(sapply(overlaps, length)==0)
      ## asign these child polygons a new trackID in the master data.frame
      dat[dat$index %in% tempNextYear[noParentRowNums,"index"]$index, #get the
          # unique index for each parent-less child polygon
          "trackID"] <-  paste0(unique(dat$sp_code_6), ## get the 6 letter code
                                # for this species
                                "_", inv[i], ## get the year in which child poly recruited
                                "_", c(1:length(noParentRowNums))) # get a vector
      # of numbers as long as the number of parentless child polys)

      # #test
      # testDat <- dat[is.na(dat$trackID)==FALSE &
      #                  dat$trackID %in% unique(dat$trackID)[1:20] &
      #                  dat$trackID %in% dat$trackID[duplicated(dat$trackID)]
      #                ,]
      #
      # ggplot(dat = testDat) +
      #   geom_sf(aes(lty = as.factor(testDat$Year), fill = trackID, alpha = Year)) +
      #   #scale_fill_discrete(guide = FALSE) +
      #   scale_alpha(range = c(.2,.8), guide = FALSE) +
      #   xlim(c(0,1)) +
      #   ylim(c(0,1))


      ## NO CHILD--becomes a ghost
      ## get the data for the current year (parent) polygons that don't have
      # any 'children' (get the row index numbers of every parent poly that is
      # included in the 'overlaps' dataset, then find the row index numbers from
      # tempCurrentYear that are not included in that list)
      tempGhosts <- tempCurrentYear[!1:nrow(tempCurrentYear) %in%
                                      unique(unlist(overlaps)),]
      ## add a '1' in a 'ghost' column, indicating that these polys are 'ghosts
      tempGhosts$ghost <- 1
      ## CHECK DORMANCY -- need to check that the gap in years between 'parent'
      # and 'child' does not exceed the dormancy argument. If the difference
      # between the year of the 'ghost' and the year to which we will be comparing
      # in the next iteration of the loop ('child' polygons in year i+1), then we
      # do not include it as a 'ghost'
      ## remove those 'ghosts' for which there is a gap between 'parent' year and 'child' year in next loop that exceeds the dormancy argument
      tempGhosts <- tempGhosts[!(inv[i+1] - tempGhosts$Year > (dorm + 1)),]

      ###AES### include an option for tie-breaking (if statements that are initiated by a user-defined option) for area vs. distance

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
