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

  ## put this trackID information in the master data.frame for this species
  # using the unique index information
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
      ## see if there is any overlap between the tempNextYear data and the
      # tempCurrentYear data (buffered)
      overlaps <- st_intersects(tempCurrentBuff, tempNextYear, sparse = FALSE)
      # returns a matrix, where the rows are the parent row index numbers, and
      # the columns are the row index numbers of the child polygons. A 'TRUE'
      # value indicates that those polygons overlap

      ## get the trackID of each 'parent' polygon, and put them as 'names' onto
      # the 'overlaps' list, but first combine it with the 'unique index', so we
      # don't have multiple list elements with the same name
      rownames(overlaps) <- paste0(tempCurrentYear$trackID, "__",
                                   tempCurrentYear$index)
      ## get the genetIDs + unique index of each 'child' polygon
      colnames(overlaps) <-paste0("genet",tempNextYear$genetID,"__",
                                  tempNextYear$index)


      ## trying to get the amount of overlap between each polygon
      overlapArea <- st_intersection(tempCurrentBuff, tempNextYear)
      ## get the trackID names merged with the index value (same values used
      # for names of rows in 'overlaps' matrix)
      overlapArea$parentName <- paste0(overlapArea$trackID, "__",
                                       overlapArea$index)
      ## get the genetID names merged with the index value (same values used for
      # names of columns in 'overlaps' matrix)
      overlapArea$childName <- paste0("genet", overlapArea$genetID.1, "__",
                                      overlapArea$index.1)

      ## calculate the overlap between each parent poly and each child poly
      overlapArea$overlappingArea <- st_area(overlapArea$geometry)

      ## aggregate the matrix by unique trackID (over rows), and then by unique
      # genetID (over columns)

      for(l in tempCurrentYear$trackID) {
        duplicateTrackIDrows <- overlaps[which(sapply(strsplit(rownames(overlaps),
                 "__"), unlist)[1,]==l),] ## get those rows that have the same
        # trackID (are genets)
        overlappingCols <- ifelse(is.vector(duplicateTrackIDrows)==TRUE, ## test
                               names(which(duplicateTrackIDrows==TRUE)),   ## yes

          names(which(apply(duplicateTrackIDrows,1, ## margin of the apply (rows)
                                             function(x)
            sum(x) >= 1)
            ==1)) ## no
          ) ## if there is at least 1 overlap, return a 'TRUE'

        overlaps
      }

      unlist(strsplit(rownames(overlaps), "__"))



      ggplot()+
        geom_sf(data = tempCurrentBuff[c(1,3,4,5,6),], aes(fill = trackID), lty = 1, alpha = .5) +
        scale_fill_discrete(guide = FALSE) +
        theme_classic() +
        geom_sf(data = tempNextYear[c(1,2,5,14,16,18,4,7,8),], aes(fill = as.factor(genetID)), lty = 2, alpha = .8) +
        xlim(c(0,1)) +
        ylim(c(0,1)) +
        geom_sf(data = tempNextYear, aes(), fill = "orange", alpha = .3)




      ## UNAMBIGUOUS PARENT: 1 parent polygon only has 1 child polygon
      unambigParentRowNums <- which(sapply(overlaps, length)==1) ## get the row
      # index numbers of 'parent' polygons that have only one 'child' polygon

      ## get the row index, index , and trackID numbers of the tempNextYear
      # (child) polygons that only have one parent, and the tempCurrentYear
      # (parent) polygons that only have one child
      unambig <- data.frame("parentRowIndex" = oneParentRowNums,
                                  "parentIndex" = tempCurrentYear[unambigParentRowNums,"index"]$index,
                                  "parentTrackID" = tempCurrentYear[unambigParentRowNums, "trackID"]$trackID,
                            "parentGenetID" = tempCurrentYear[unambigParentRowNums, "genetID"]$genetID,
                                  "childRowIndex" = sapply(overlaps[unambigParentRowNums], unlist
                                  ),
                                  "childIndex" = sapply(overlaps[unambigParentRowNums], function(x)
                                    tempNextYear[unlist(x),"index"]$index),
                            "childGenetID" = sapply(overlaps[unambigParentRowNums], function(x)
                              tempNextYear[unlist(x),"genetID"]$genetID))

      unambigParentDat <- tempCurrentYear[oneParentRowNums,]  ## get the data for these unambiguous parents
      unambigChildDat <- tempNextYear[unambigParent$childRowIndex,]
      ## get the data for the unambiguous children

      mapview(tempCurrentYear) + mapview(tempNextYear, col.regions = "red") +
      mapview(unambigParentDat, col.regions = "blue") + mapview(unambigChildDat, col.regions = "orange")






      # polygon assign trackIDs to polygons from year t+1 that #
      # overlap only one polygon from year t

      oneParentRowNums <- which(sapply(overlaps, length)==1)

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
