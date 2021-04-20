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

assign <- function(){
  # arguments ---------------------------------------------------------------
  dat <- st_sf(sampleDat) # data in 'grasslandData' format, must be an sf object
    ## add columns for trackID, age, size_t+1, and recruit
    dat$trackID <- NA
    dat$age <- NA
    dat$size_tplus1 <- NA
    dat$recruit <- NA
    dat$survives_tplus1 <- NA
  inv <- sort(sampleInv) # integer vector of quadrat sampling years in
  # sequential order
  dorm <- 1 # dormancy allowed by the function
  buff <- .05 # buffer of movement allowed from year to year, in meters
  overlap <- .50 # the percentage of overlap (in decimal form) between focal
  # indiv. and next year indiv. that will be required to consider them both
  # the same individual (not 100? sure on this one yet...)
  # work -------------------------------------------------------------------
  # generate a vector of integer numbers to use as track IDs (this
  # variable will be continually re-defined? through the for-loop
  # so that trackIDs won't be repeated?)
  IDs <- c(1:1000000)


  # buffer the data for all of the polygons (NOTE: this doesn't have trackID
  # info! Or any data about recruit or age!)
  dataBuff <- st_buffer(dat, buff)

  # probably use a for loop ...
  # i = year in inventory
  for (i in seq_along(inv)) {
    # check if this the first year of sampling, or if this is the year
    # immediately following a gap in sampling.
    if (i == 1) { #is this the first year of sampling?
      #if this is the first year of sampling, then get the raw data for
      # year i, from which we will select one individual to focus on
      dat[dat$Year==inv[i],] #data for current [i] year
      #If this is the first year (i = 1), then all individuals get an NA
      # in the 'age' and 'recruit' columns.
      dat[dat$Year==inv[i],'age'] <- NA
      dat[dat$Year==inv[i],'recruit'] <- NA
      dat[dat$Year==inv[i],'trackID'] <- IDs[1:nrow(dataCurrent)] ## populate
      # trackIDs here redefine the master trackID vector so that the trackIDs
      # don't repeat
      IDs <- IDs[-(1:nrow(dataCurrent))]
    }

    ## j = individual unique trackID in the dataset in that year (only or
    # individuals that already have a trackID assigned)
    for (j in dat[dat$Year==inv[i],'trackID']) { ## loop through each of the
      # uniquetrackIDs. If there are more than 1 polygon w/ the same trackID,
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
  }
}
# output ---------------------------------------------------------------


