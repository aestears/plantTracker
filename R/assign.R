#' tracks individual plants through time
#' @param dat
#' @param inv
#' @param dorm
#' @param buff
#' @param buffGenet
#' @param clonal
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @import tidyr sf dplyr


# example input data ------------------------------------------------------
# grasslandData (or exact same format), subset to a unique site, quad,
# species
load("./data/grasslandData.rda")
load("./data/grasslandInventory.rda")
# prepares the dataset to feed into the 'assign' function (the 'Assign'
# function will do this ahead of time when the user calls it)
sampleDat <- grasslandData[grasslandData$Site == "KS"
                           & grasslandData$Quad == "q33"
                           & grasslandData$Species == "Solidago mollis",]
# this should be a data.frame
dat <- sampleDat

# get the appropriate grasslandInventory data for the "unun_11" quadrat,
# to tell the 'assign' function when the quadrat was sampled
sampleInv <- grasslandInventory[["q33"]]
# this should be an integer vector
inv <- sampleInv

assign <- function(dat, inv, dorm, buff, buffGenet, clonal,...){
  # arguments ---------------------------------------------------------------
  ## double check the format of the inputs, and add additional columns required
  dat <- sf::st_sf(dat) # data in 'grasslandData' format, must be an sf object
    ## add columns for trackID, age, size_t+1, and recruit
    dat$trackID <- NA
    dat$age <- NA
    dat$size_tplus1 <- NA
    dat$recruit <- NA
    dat$survives_tplus1 <- NA
    dat$index <- c(1:nrow(dat))## assign an arbitrary, unique index number to
    # each row in the dataset
    dat$ghost = NA
    dat$rametArea = NA
  ## inv is a vector of the years in which the target quadrat is sampled
  inv <- sort(inv) ## integer vector of quadrat sampling years in
  # sequential order
  ## make sure that the 'assignOut' data.frame that contains the output is empty
  if(exists("assignOut")){
    rm(assignOut)
  }
  # ## test user-defined args
  # ## user-defined arguments
  # dorm <- 0 ## dormancy allowed by the function
  # buff <- .05 ## buffer of movement allowed from year to year, in meters
  # buffGenet <- 0.001 ## buffer between polygons (/2) that is the maximum allowed for
  # # them to be considered genets
  # clonal <- 1 ## binary option that indicates whether this species is allowed to
  # # break into multiple polygons for one individual

  ## work -------------------------------------------------------------------
  ## find the first year in the dataset that was actually sampled
  firstDatYear <- min(dat$Year)
  ## find the index of the first year in teh quadrat inventory
  firstYearIndex <- which(inv==firstDatYear)
  ## get the dataset for the first year of actually sampling
  tempCurrentYear <- dat[dat$Year==firstDatYear,]
  ## assign genetIDs and trackIDs to the first year of sampling dataset
  if (clonal==1) { ## if this species is 'clonal' (clonal argument == 1), use
    # Dave's groupByGenet function to give a singleTrackID to a group of
    # polygons if they are close enough to one another (user-defined)
    tempCurrentYear$genetID <- groupByGenet(tempCurrentYear, buffGenet) # put
    # these temporary genet ID's in a temporary 'genetID' column, which will
    # tell us to give these genets the same track ID
    ## aggregate size by genetID (total size for each ramet)
    tempCurrentYear <- sf::st_join(tempCurrentYear, (tempCurrentYear %>%
                                                  dplyr::group_by(genetID) %>%
                                               dplyr::summarize(rametAreaTemp = sum(Area))), by = "genetID") %>%
      dplyr::select(-c(genetID.x, rametArea)) %>%
      dplyr::rename("genetID" = "genetID.y", "rametArea" = "rametAreaTemp")
  }
  if (clonal == 0) { ## if genets are not allowed (clonal == 0), then each
    # individual polygon gets a unique value in the 'genetID' column
    tempCurrentYear$genetID <- 1:nrow(tempCurrentYear)
    tempCurrentYear$rametArea <- tempCurrentYear$Area
  }

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

  tempCurrentYear<- merge(tempCurrentYear[,names(tempCurrentYear) != "trackID"], IDs, by = "genetID")

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
  for (i in (firstYearIndex+1):28) {
       #length(inv)) {
    ## CHECK IF YEARS ARE CONTINUOUS -- check to see if the sampling years of
    # 'tempCurrentYear' and 'tempNextYear' are not far enough apart to exceed
    # the 'dormancy' argument. If dormancy is not exceeded, then go ahead with
    # this loop. If it is, then freshly redefine 'tempCurrentYear' and proceed
    # to the next 'i'
    if (inv[i] - inv[i-1] <= (dorm+1)) {
      ## 'tempCurrentYear' is the sf data.frame of the 'current' year
      ## need to get the sf data.frame of the 'next' year (year 'i')
      tempNextYear <- dat[dat$Year==inv[i],]

      ## MAKE SURE THERE IS DATA IN YEAR i-1 (tempCurrentYear isn't empty) If
      # not, then go to the next i (only after replacing 'tempCurrentYear' with
      # the tempNextYear data.frame)
      if (nrow(tempCurrentYear)>0) {
        ## add a buffer to the current year data
        tempCurrentBuff <- sf::st_buffer(tempCurrentYear, buff)

        ## MAKE SURE THERE IS DATA IN YEAR i (tempNextYear isn't empty) If not,
        # then go to next i
        if (nrow(tempNextYear)>0) {
          ## AGGREGATE BY GENET for year i (if clonal = 1)
          if(clonal==1) {
            tempNextYear$genetID <- groupByGenet(tempNextYear, buffGenet)
            ## aggregate size by genetID (total size for each ramet)
            tempNextYear <- sf::st_join(tempNextYear, (tempNextYear %>%
                                  dplyr::group_by(genetID) %>%
                                  dplyr::summarize(rametAreaTemp = sum(Area))),
                                    by = "genetID") %>%
              dplyr::select(-c(genetID.x, rametArea)) %>%
              dplyr::rename("genetID" = "genetID.y",
                            "rametArea" = "rametAreaTemp")
          }
          ## assign unique genetIDs for every polygon (if clonal = 0)
          if(clonal==0) {
            tempNextYear$genetID <- 1:nrow(tempNextYear)
            tempNextYear$rametArea <- tempNextYear$Area
          }

          ## FIND OVERLAPPING POLYGONS
          ## trying to get the amount of overlap between each polygon
          overlapArea <- sf::st_intersection(tempCurrentBuff, tempNextYear)

          ## get the trackID names merged with the index value (same values used
          # for names of rows in 'overlaps' matrix)
          overlapArea$parentName <- paste0(overlapArea$trackID, "__",
                                           overlapArea$index)
          ## get the genetID names merged with the index value (same values used
          # for names of columns in 'overlaps' matrix)
          overlapArea$childName <- paste0("genet__",
                                          overlapArea$genetID.1,
                                          "__",
                                          overlapArea$index.1)

          ## calculate the overlap between each parent poly and each child poly
          overlapArea$overlappingArea <- sf::st_area(overlapArea$geometry)

          overlapArea <- sf::st_set_geometry(overlapArea, NULL)

          ## aggregate the overlaps by rows (by 'parents')
          overlaps <- overlapArea %>%
            dplyr::group_by(trackID, genetID.1) %>%
            dplyr::summarize(overlappingArea = sum(overlappingArea))

          names(overlaps) <- c("parentTrackID", "childGenetID",
                               "overlappingArea")

          overlaps$childGenetID <- paste0("genet__", overlaps$childGenetID,
                                          "__",seq(from = 500, by = 1,
                                                   length.out = nrow(overlaps)))

          ## transform the overlaps dataframe into a 'wide' dataframe, so that
          # every row is a unique parent trackID, and the columns are children
          # (not completely aggregated yet). The value is the overlap between
          # each parent/child pair
          overlaps <- as.data.frame(tidyr::pivot_wider(overlaps,
                                      id_cols = c(parentTrackID, childGenetID),
                                                names_from = childGenetID,
                                                values_from = overlappingArea))
          temp <- as.data.frame(overlaps[,2:ncol(overlaps)])
          rownames(temp) <- as.character(overlaps$parentTrackID)
          names(temp) <- names(overlaps)[2:ncol(overlaps)]
          overlaps <- temp

          ## transpose the 'overlaps' matrix so that we can aggregate by genet
          # (genetID are rows, trackID are columns)
          overlaps <- t(overlaps)
          ## simplify the row names so they only have the genetID number
          overlaps <- as.data.frame(overlaps)
          overlaps$genetID <- as.numeric(sapply(strsplit(rownames(overlaps),
                                                         "__"),
                                                unlist)[2,])
          overlaps[is.na(overlaps)==TRUE] <- 0

          overlaps <- as.data.frame(overlaps %>%
            dplyr::group_by(genetID) %>%
            dplyr::summarize(across(names(overlaps)[1:(ncol(overlaps)-1)],
                                    sum)))

          temp <- as.data.frame(overlaps[,2:ncol(overlaps)])
          rownames(temp) <- paste0("genet_",overlaps$genetID)
          names(temp) <- names(overlaps)[2:ncol(overlaps)]
          overlaps <- temp

          ## transpose the matrix back so that each row is a parent and each
          # column is a child
          overlaps <- t(overlaps)

          ## each parent can only have one child, and each child can only have
          # one parent (since we've already clustered the polygons by genet).
          # Each parent-child pair will be determined by the greatest amount of
          # overlap. This is done by going through the 'overlaps' matrix and
          # finding the largest overlap value. Then, the trackID from that
          # parent will be assigned to the child.

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
            maxChild <- strsplit(colnames(
              whileOverlaps)[maxInds[1,2]],"_")[[1]][2]

            ## put the trackID from the parent into the tempNextYear dataframe
            # rows that correspond to the genetID of the child
            tempNextYear[tempNextYear$genetID==maxChild,"trackID"] <- maxParent

            ## overwrite the 'max' value with an NA, so we can find the next
            # largest value
            whileOverlaps[maxInds[1,1],maxInds[1,2]] <- NA
            ## overwrite all of the other values in the parent row and the child
            # column with NAs also (since each parent can only have one child,
            # and each child can only have one parent)
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
          orphans <- tempNextYear[is.na(tempNextYear$trackID)==TRUE,]
          ## make sure that there are orphans
          if (nrow(orphans)>0) {
            ## make a unique trackID for each genetID in the 'orphan' dataset
            orphanIDs <- data.frame(
              "genetID" = unique(tempNextYear[is.na(tempNextYear$trackID)==TRUE,
                                              "genetID"]$genetID), ## get the
              # unique genetIDs of the 'orphans'
              "trackID" = paste0(unique(
                tempNextYear$sp_code_6), ## get the unique 6-letter species code
                "_",unique(tempNextYear$Year), ## get the unique year
                "_",c(1:length(unique(tempNextYear[is.na(
                  tempNextYear$trackID)==TRUE,"genetID"]$genetID)))))

            ## add the orphan trackIDs to the 'orphan's data.frame
            orphans <- dplyr::left_join(orphans, orphanIDs, by = "genetID") %>%
              dplyr::select(-trackID.x) %>%
              dplyr::rename("trackID" = "trackID.y")
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
          ## (first have to assign the parents first)
          ## define the PARENTS data.frame
          parents <- tempCurrentYear[tempCurrentYear$trackID %in%
                                       tempNextYear$trackID,]
          ## define the children data.frame
          children <- tempNextYear[tempNextYear$trackID %in%
                                     tempCurrentYear$trackID,]

          ## ASSIGN DEMOGRAPHIC DATA TO CHILDREN
          if (nrow(children)>0) {
            ## give children a 0 in the recruit column, since they have a parent
            # (as long as the parent doesn't have an NA--was recruited in a year
            # when we couldn't know when it was recruited)
            if(inv[i] - inv[i-1] <= 1) {
              children$recruit <- 0
            } else {
              children$recruit <- NA
            }
            ## give the children the appropriate age column (only if the parents
            # don't have an 'NA' for age) (parent's age + 1)
            tempParents <- sf::st_set_geometry(parents, NULL) %>%
              dplyr::select(trackID, age)  %>%
              dplyr::mutate(age = age +(inv[i] - inv[i-1])) %>%
              dplyr::rename("trackIDtemp" = "trackID")

            children <- children %>%
              dplyr::select(-c(age)) %>% ## remove the 'age' column
              dplyr::left_join(tempParents, by = c("trackID"="trackIDtemp"))
            ## add the'age'column from the appropriate parents (+1), joined by
            # trackID
          }

          ## PARENTS
          ## ASSIGN DEMOGRAPHIC DATA FOR PARENTS
          ## make sure that there is actually a parents data.frame
          if (nrow(parents)>0) {
            ## PARENTS ('parents' data.frame + 'deadGhosts' data.frame)
            ## give the 'parents' a '1' for survival
            parents[,"survives_tplus1"] <- 1
            ## give the 'parents' a 0 in the 'ghosts' column
            parents[,"ghost"] <- 0
            ## assign size in the next year for parents data.frame
            childSizeTemp <- sf::st_set_geometry(
              children[,c("trackID", "rametArea")],NULL) %>%
              dplyr::rename("size_tplus1" = "rametArea",
                            "trackIDtemp" = "trackID") %>%
              dplyr::group_by(trackIDtemp) %>% ## aggregate size by trackID
              dplyr::summarize("size_tplus1" = mean(size_tplus1))

            parents <- parents %>%
              dplyr::select(-size_tplus1) %>%
              dplyr::left_join(childSizeTemp, by = c("trackID" = "trackIDtemp"))
          }

          ## ONLY ALLOW GHOSTS IF DORMANCY > 1
          if (dorm>0) {
            ## GHOSTS: parent polygons that don't have 'children' in year i. If
            # they were observed in a year that is more than dorm+1 years prior
            # to year i, then they don't get saved to the next year. If they were
            # observed in a year that is >= to dorm+1 years prior to year i, then
            # they get added to the 'tempNextYear' data.frame, which both go into the
            # 'tempCurrentYear' data.frame before the next iteration of the loop
            ## get the ghost individuals
            ghostsTemp <- tempCurrentYear[!(tempCurrentYear$trackID %in%
                                              tempNextYear$trackID),]
            ## check that these indivduals can be 'ghosts' in the next year (if the
            # gap between the year of their observation and year i+1 is greater than
            # the dormancy argument (+1), then) they are not ghosts, and get a 0 for
            # survival
            ghosts <- ghostsTemp[((inv[i+1] - ghostsTemp$Year) <= (dorm + 1)),]
            ## put the 'ghosts' that exceed the dormancy argument into their
            # own data.frame
            deadGhosts <- ghostsTemp[((inv[i+1] - ghostsTemp$Year) > (dorm + 1)),]
            ## ASSIGN DEMOGRAPHIC DATA TO GHOSTS
            ##give these survived ghosts a '1' in the 'ghost' column
            if (nrow(ghosts)>0) {
              ghosts$ghost <- 1
            }
            ## give the deadGhosts a 0 for survival (if there are any deadGhosts)
            if(nrow(deadGhosts)>0){
              deadGhosts[,"survives_tplus1"] <- 0
            }
            ## get the data.frames in a consistent format
            ghosts <- ghosts[,names(children)]
            deadGhosts <- deadGhosts[,names(children)]
          } else { ## any trackID that doesn't have a child in year 'i' is 'dead'
            deadGhosts <- tempCurrentYear[!(tempCurrentYear$trackID %in% tempNextYear$trackID),]
            deadGhosts$survives_tplus1 <- 0
            deadGhosts <- deadGhosts[,names(children)]
            ghosts <- NULL
          }
          ## PREPARE FOR NEXT i
          ## arrange columns of children, orphans, and ghosts into the same order
          orphans <- orphans[, names(children)]
          parents <- parents[,names(children)]
          ## bind children, orphans, and ghosts into one data.frame, that will become
          # the data for the current year in the next iteration of the loop
          tempNextYear <- rbind(children, orphans, ghosts)

          ## STORE PARENTS DEMOGRAPHIC DATA
          if (exists("assignOut")==FALSE) { ## if the assignOut data.frame is empty
            assignOut <- rbind(parents, deadGhosts)
          }
          if (exists("assignOut")==TRUE) { ## if this is not the first year, then add demographic data
            assignOut <- rbind(assignOut, parents, deadGhosts)
          }
          tempCurrentYear <- tempNextYear
        } ## end of 'if' statement that determines if the tempNextYear data
        # exists
        else { ## if the tempNextYear data does not exist, then keep the
          # tempCurrentYear data.frame.
          ## Check if the dormancy is exceeded
          ifelse (test = inv[i+1] - inv[i] <= (dorm+1), ## check if dorm. is ok
            yes = (tempCurrentYear$ghost <- 1), ## if dormancy is okay, put a 1
            # in the ghost column of tempCurrentYear
            no = (tempCurrentYear <- NA) ## if dormancy is exceeded, then make
            # the 'tempCurrentYear' data.frame empty
          )
        } ## end of 'else' that has steps if tempNextYear is empty
      } ## end of 'if' that determines if the tempCurrentYear data exists
      else { ## what to do if 'tempCurrentYear' doesn't exist
        tempCurrentYear <- tempNextYear
        }
      } ## end of 'if' statement that determines if gap between inv[i-1] and
    # inv[i] is less than or equal to  dorm+1
    else { ## if the gap between years exceeds the dormancy argument
      tempCurrentYear <- dat[dat$Year==inv[i],]
      ## group by genet
      if (clonal==1) { ## if this species is 'clonal' (clonal argument == 1), use
        # Dave's groupByGenet function to give a singleTrackID to a group of
        # polygons if they are close enough to one another (user-defined)
        tempCurrentYear$genetID <- groupByGenet(tempCurrentYear, buffGenet) # put these temporary
        # genet ID's in a temporary 'genetID' column, which will tell us to give
        # these genets the same track ID
        ## aggregate size by genetID (total size for each ramet)
        tempCurrentYear <- sf::st_join(tempCurrentYear, (tempCurrentYear %>%
                                         dplyr::group_by(genetID) %>%
                                         dplyr::summarize(rametAreaTemp = sum(Area))), by = "genetID") %>%
          dplyr::select(-c(genetID.x, rametArea)) %>%
          dplyr::rename("genetID" = "genetID.y", "rametArea" = "rametAreaTemp")
      }
      if (clonal == 0) { ## if genets are not allowed (clonal == 0), then each
        # individual polygon gets a unique value in the 'genetID' column
        tempCurrentYear$genetID <- 1:nrow(tempCurrentYear)
        tempCurrentYear$rametArea <- tempCurrentYear$Area
      }
      ## assign a unique trackID to every unique genetID in this first-year dataset
      IDs <- data.frame( "genetID" = sort(unique(tempCurrentYear$genetID)), ## get a vector
                         # of all of the unique genetIDs
                         "trackID" = paste0(unique(dat$sp_code_6), ## get the unique 6-letter
                                            # species code
                                            "_",unique(tempCurrentYear$Year), ## get the unique year
                                            "_",c(1:length(unique(tempCurrentYear$genetID))))) ## get a vector of
      # unique numbers that is the same length as the genetIDs in this quad/year

      tempCurrentYear<- merge(tempCurrentYear[,names(tempCurrentYear) != "trackID"], IDs, by = "genetID")
      } ## end of exceeded dormancy 'else'
    } ## end of loop i
  ## clean up output data.frame
  assignOut <- assignOut[is.na(assignOut$Species)==FALSE,] %>%
    dplyr::select(-ghost, -genetID, -index )
# output ---------------------------------------------------------------
return(assignOut)
}




# testing -----------------------------------------------------------------
# ggplot() +
#   #geom_sf(data = deadGhosts, fill = "red") +
#   geom_sf(data = parents, aes(fill = as.factor(trackID)), alpha = 0.8) +
#   geom_sf(data = children, aes(fill = as.factor(trackID)), alpha = 0.6) +
#   geom_sf(data = orphans, fill = "grey") +
#   scale_fill_discrete(guide = FALSE) +
#   theme_classic()
#
# multiples <- assignOut[assignOut$trackID %in% names(which(table(assignOut$trackID)>1)),]
#
# ggplot(st_buffer(multiples,0.02)) +
#   geom_sf(aes(fill = trackID), alpha = 0.6) +
#   scale_fill_discrete(guide = FALSE) +
#   theme_classic()
#
# ggplot() +
#   geom_sf(data = st_buffer(assignOut[assignOut$trackID==unique(assignOut$trackID)[2],],.02), aes(fill = as.factor(Year)), alpha = 0.6) +
#   geom_sf(data = st_buffer(dat[dat$Year==2007,], .02)) +
#   lims(x = c(0,1), y = c(0,1)) +
#   theme_classic()

testOutput <- assign(sampleDat, sampleInv, 1, .05, .001, 0)

ggplot(st_buffer(testOutput,.02)) +
  geom_sf(aes(fill = as.factor(Year)), alpha = 0.5) +
  #scale_fill_discrete(guide = FALSE) +
  theme_classic()

