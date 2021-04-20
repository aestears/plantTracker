#' tracks indivudal plants through time
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
      dat[dat$Year==inv[i],'age'] <- -9999
      dat[dat$Year==inv[i],'recruit'] <- -9999
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

        #getting the data for year
        # along the sequence of sampled years i (year in which focal plant
        # is in) + k (the number along a sequence of values up to (1+ the
        # value of dorm))]. This selects only values for the sampled years, not
        # necessarily in a numerical order of year + 1!! We will then compare
        # the buffered focal plant to this dataset to see if there
        # is any overlap (i.e. there should be no empty dataNext data.frames)

        ## is there adequate overlap between the unbuffered polygon and
        # the next year?
        if (sum(st_intersection(plant, dat[dat$Year==inv[i + k],])$Area.1)/plant$Area > overlap) {
          # is the total area of the overlapping next year polygons
          # >= to the 'overlap' argument?
        } else {
          if (sum(st_intersection(plantBuff, dataNext)$Area.1)/
              plantBuff$Area > overlap) { ## if not, is there adequate overlap
            # between the buffered polygon and the next year?
          } else { ## if there isn't overlap between either the buffered
            # or the un-buffered polygon and any polygons in the next year,
            # can we go to the next year[k]?
            ## Has the quadrat been sampled in year k+1 (along a sequential
            # vector of years, not necessarily along the inv. vector of sampled years)?

            if((inv[i] + (k+1)) %in% inv) { ## is the next year in numerical
              # sequence: plant[j]$yr + (k+1) present in the 'inventory' vector?

              ## If yes, this quad was sampled. Iss the gap between year of the
              # focal plant (inv[i]) and year k+1 inside the 'dormancy' arg?
              if () { #if this is outside the dormacy argument, then give the
                ## NO: then give the focal plant an NA in the 'surv' and
                # 'size_t+1' columns
              } else {
                ## YES:  then go to the next k
              }
            } else {
              ## if not, then this quad was NOT sampled in this year
              ## if k = dorm+1 (if the )
            }


            ## If it was NOT sampled,
###AES###
          }


        }

        overlapNoBuff <-  # see if there
          overlapNoBuffArea <- st_intersection(plant, dataNext)
        # is overlap of the un-buffered polygon with anything
        ## st_intersection gives the area of the focal individual that
        # is overlapped... doesn't indicate which ones in year k+1 are
        # doing the overlapping ... problem

        #use st_overlaps? or st_intersects?
        #st_intersects gives you the index of the overlapping polygons!
        mapview(plant) + mapview(dataNext[c(30,36,37),])
        # what is the percentage of overlap?
        overlapNoBuff
        if(nrow(overlapNoBuff)>0) { #if there IS overlap between unbuffered polygon and something in year k, then assign those w/ the overlap the same track ID

          ###AES###
        } else {

        }
        ## if not, see if there is overlap of the buffered polygon with anything
        mapview(plant, col.regions = "green") + mapview(plantBuff, col.regions = "purple") + mapview(dataNext)
        ## if there is a solution in k=1, then STOP; but if there is not a solution in k =1, then continue the loop. How do we do this??
        ## is there a way to stop a for-loop if a certain condition is met??
      }

    }
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


