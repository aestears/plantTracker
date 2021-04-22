#' tracks individual plants through time
#'
#' @return
#' @export
#'
#' @examples

# required packages -------------------------------------------------------
library(sf)

# example input data ------------------------------------------------------
# grasslandData (or exact same format), subset to a unique site, quad,
# species
load("~/PlantTracker/data/grasslandData.rda")
load("~/PlantTracker/data/grasslandInventory.rda")
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
  inv <- sort(sampleInv) ## integer vector of quadrat sampling years in
  # sequential order

  ## user-defined arguments
  dorm <- 1 ## dormancy allowed by the function
  buff <- .05 ## buffer of movement allowed from year to year, in meters
  overlap <- .50 ## the percentage of overlap (in decimal form) between focal
  # indiv. and next year indiv. that will be required to consider them both
  # the same individual (not 100% sure on this one yet...)
  clonal <- 1 ## binary option that indicates whether this species is allowed to
  # break into multiple polygons for one individual

  ## work -------------------------------------------------------------------

  datFirst <- dat[dat$Year==inv[1],] ## get the data just for the first year of
  # sampling

  ## assign a unique trackID to every polygon in this first-year dataset
  numbs <- c(1:nrow(datFirst)) ## get a vector of unique numbers
  datFirst$trackID <-  paste0(unique(dat$sp_code_6),"_",unique(datFirst$Year),
                     "_",numbs) ## add data for species and year of recruitment

  ## assign the first-year data to the 'tempCurrentYear' data.frame
  tempCurrentYear <- datFirst ## this data.frame will get redefined for
  # each iteration of the for-loop below

  ##  i = year in inventory
  for (i in 2:length(inv)) {

    ## 'tempCurrentYear' is the sf data.frame of the 'current' year
    tempCurrentBuff <- st_buffer(tempCurrentYear, buff) ## need to add a buffer
    # to this data.frame

    ## need to get the sf data.frame of the 'next' year
    tempNextYear <- dat[dat$Year==inv[i],]

    ## see if there is any overlap between the tempNextYear data and the
    # tempCurrentYear data (buffered)
    overlaps <- st_intersects(tempNextYear, tempCurrentBuff)
    # returns a list, where each object in the list is the row index of the
    # polygon in tempNextYear, and the contents of the list are row indices of
    # the overlapping polygons from tempCurrentBuff. Ultimately, we must chose

    ## are there any polygons from tempNextYear that overlap with more than one
    # polygon from tempCurrentYear?
    ## get the 'index' numbers of the tempCurrentYear polygons that overlap with
    # the same t+1 polygons
    currentOverlapRowNums <- as.numeric(names(which(table(unlist(
      overlaps[which(sapply(overlaps, length)>1)]))>1)))
    currentYearOverlaps <- tempCurrentBuff[currentOverlapRowNums,]

    ## get the 'index' numbers of the tempNextYear polygons that overlap with
    # more than one polygon from year t
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
      #row for the year t+1 polygon

    } #end of loop j
    ## then deal with the polygons from year t+1 that don't overlap at all with
    # the buffered polygons from year t
    ###AES### what the heck now?
    ## j = individual unique trackID in the dataset in that year (only or
    # individuals that already have a trackID assigned)
    for (j in dat[dat$Year==inv[i],'trackID']) { ## loop through each of the
      # unique trackIDs. If there are more than 1 polygon w/ the same trackID,
      # then it should include both of them. Ideally this should ignore
      # those individuals w/out a trackID? (###AES###)
      plant <- dat[dat$Year==inv[i] & dat$trackID %in% j,] ## getting data for
      # individual(s) j, which we will then send through the tracking function
      plantBuff <- dataBuff[dataBuff$Area %in% plant$Area,] ## select
      # the buffered polygon(s) for the appropriate individual
      for (k in 1:(dorm+1)) { # this loop compares the focal plant to
        # the quadrat data in the the next year. This loop will repeat 1 +
        # as many times as the 'dorm' argument
        # (compare to year t+1, ..., t+dorm).

        ## data for next sampled year: dat[dat$Year==inv[i + k],]

        ## is the quadrat sampled in year inv[i] + k ???
        if ((inv[i] + k) %in% inv) { ## quadrat IS measured in year inv[i] + k
          ## get overlapping area between unbuffered plant and next year
          # (excluding next year polygons that already have a trackID)
          overlapsNoBuff <- st_intersection(plant, dat[dat$Year==inv[i + k]
                                              & is.na(dat$trackID)==TRUE,])
          ## get the complete raw data for the overlapping polygons
          overlapsNoBuffRaw <- dat[dat$Year==inv[i + k] &
                                     is.na(dat$trackID)==TRUE ,][(st_intersects(
                                       plant, dat[dat$Year==inv[i + k] &
                                                    is.na(dat$trackID)==TRUE,],
                                       sparse = FALSE)[1,]),]

          ## is there a match w/ polygons in the next year w/out focal plant
          # buffer? That is above the threshold of the 'overlaps' argument? That
          # also doesn't already have a trackID?
          if (sum(overlapsNoBuff$Area.1)/plant$Area > overlap) {
            ## if YES,  give to focal individual: surv = 1, size_tplus1 = #, trackID for
            # next year polygons, age for next year polygons, recruit for next
            # year polygons
            dat[dat$Year==inv[i] & dat$trackID %in% j, 'size_tplus1'] <-
              sum(overlapsNoBuffRaw$Area)
            dat[dat$Year==inv[i] & dat$trackID %in% j, 'surv'] <- 1

            dat[dat$Year==inv[i + k] & dat$Area %in% overlapsNoBuffRaw$Area,
                'trackID'] <- dat[dat$Year==inv[i] & dat$trackID %in% j,
                    'trackID'] ## give the trackID from current individual to
                      # the overlapping polys from the next year
            dat[dat$Year==inv[i + k] & dat$Area %in% overlapsNoBuffRaw$Area,
                'age'] <- dat[dat$Year==inv[i] &dat$trackID %in% j, 'age'] + k
                      ## give the age of the focal individual + k to the
                      # overlapping polys from the next year
            dat[dat$Year==inv[i + k] & dat$Area %in% overlapsNoBuffRaw$Area,
                'recruit'] <- 0 ## give a 0 in the recruit column for the
                      # overlapping polys, since they aren't recruits
          } else {
            ## if there is NOT an overlap between unbuffered individual and
            # polys in the next year...
            ## get overlapping area between buffered plant and next year
            # (excluding next year polygons that already have a trackID)
            overlapsBuff <- st_intersection(plantBuff, dat[dat$Year==inv[i + k]
                                               & is.na(dat$trackID)==TRUE,])
            ## get the complete raw data for the overlapping polygons
            overlapsBuffRaw <- dat[dat$Year==inv[i + k] &
                                       is.na(dat$trackID)==TRUE ,][(
                          st_intersects(plantBuff, dat[dat$Year==inv[i + k] &
                          is.na(dat$trackID)==TRUE,], sparse = FALSE)[1,]),]

            ## is there a match w/ polygons in the next year w/ focal
            # plant buffer? That is above the threshold of the 'overlaps'
            # argument? That also doesn't already have a trackID?
            if (sum(overlapsBuff$Area.1)/plantBuff$Area > overlap) {
              ## if YES, give to focal individual: surv = 1, size_tplus1 = #,
              # trackID for next year polygons, age for next year polygons,
              # recruit for next year polygons
              dat[dat$Year==inv[i] & dat$trackID %in% j, 'size_tplus1'] <-
                sum(overlapsBuffRaw$Area)
              dat[dat$Year==inv[i] & dat$trackID %in% j, 'surv'] <- 1

              dat[dat$Year==inv[i + k] & dat$Area %in% overlapsBuffRaw$Area,
                  'trackID'] <- dat[dat$Year==inv[i] & dat$trackID %in% j,
                      'trackID'] ## give the trackID from current individual to
              # the overlapping polys from the next year
              dat[dat$Year==inv[i + k] & dat$Area %in% overlapsBuffRaw$Area,
                  'age'] <- dat[dat$Year==inv[i] &dat$trackID %in% j, 'age'] + k
              ## give the age of the focal individual + k to the
              # overlapping polys from the next year
              dat[dat$Year==inv[i + k] & dat$Area %in% overlapsBuffRaw$Area,
                  'recruit'] <- 0 ## give a 0 in the recruit column for the
              # overlapping polys, since they aren't recruits
            } else { ## if there is NOT an overlap between the buffered
              # invididual and polys in the next year
              if (k == dorm +1) { ## if YES (we've reached the maximum allowable
                # k for the dorm arg.), then plant is dead in next year.
                ## Give to focal individual: surv = 0, size_tplus1 = NA
                dat[dat$Year==inv[i] & dat$trackID %in% j, 'size_tplus1'] <- NA
                dat[dat$Year==inv[i] & dat$trackID %in% j, 'surv'] <- 0
              }
            }
          }
        } else { ## quadrat is NOT measured in year inv[i] + k
          ## does k = dorm + 1 (have we reached the maximum value of k according
          # to the the dormancy that's allowed by the user-defined arg?)
          if (k == dorm+1) { ## if TRUE, then this is the last year inv[i] + k
            # that is allowed by the dorm argument, but there is not data, so
            # the focal individual gets an 'NA' in both the surv and
            # size_tplus1 columns
            dat[dat$Year==inv[i] & dat$trackID %in% j, 'size_tplus1'] <- NA
            dat[dat$Year==inv[i] & dat$trackID %in% j, 'surv'] <- NA
          }
        }
      }
    } ###AES###
    # Now do the same tracking loop, but for individuals that don't yet
    # have a trackID !
    if (nrow(dataCurrent[is.na(dataCurrent$trackID)==TRUE,])>0) { #make sure
      # that there actually are individuals w/out trackIDs before moving
      # on to the next loop
      for (l in 1:nrow(dataCurrent[is.na(dataCurrent$trackID)==TRUE,])) {
        ## loop through each of the unique individuals that don't have trackIDs.
        plant <- dataCurrent[l,] # getting data for
        # one individual j, which we will then send through
        # the tracking function
        plantBuff <- dataCurrentBuff[l,] # select
        # the buffered polygon for the appropriate individual

        ## If this is the first year of sampling after a break, then give this
        # plant[l] an NA is the 'age' and 'recruit' columns
        if (inv[i] - inv[i-1] > 1) { # if the difference between year[i] and
          # year[i-1] is greater than 1, then there was a measurement gap in
          # the previous year
          ## give this individual an NA for age and recruit (likely already
          # is, but just in case...)
          plant$age <- NA
          plantBuff$age <- NA

          plant$recruit <- NA
          plantBuff$recruit <- NA
        } else { ## If this is NOT the first year af sampling after a break,
          # then give this plant[l] a '0' for 'age', and a '1' or 'recruit'
          plant$age <- 0
          plantBuff$age <- 0

          plant$recruit <- 1
          plantBuff$recruit <- 1
        }
        ## Assign this new obs. a unique trackID
        plant$trackID <-  IDs[1] # populate trackID
        plantBuff$trackID <- plant$trackID

        # redefine the master trackID vector so that the trackIDs don't repeat
        IDs <- IDs[-1]

        ## will insert the same tracking loops as above! (maybe an entire
        # different function? ) ###AES###
      }
    }
  } # end loop 'i'
}
# output ---------------------------------------------------------------


