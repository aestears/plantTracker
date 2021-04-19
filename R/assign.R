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

function(){
  # arguments ---------------------------------------------------------------
  dat <- st_sf(sampleDat) # data in 'grasslandData' format, must be an sf object
  inv <- sort(sampleInv) # integer vector of quadrat sampling years in
  # sequential order
  dorm <- 1 # dormancy allowed by the function
  buff <- .05 # buffer of movement allowed from year to year, in meters
  overap <- .50 # the percentage of overlap (in decimal form) between focal indiv. and next year indiv. that will be required to consider them both the same individual (not 100? sure on this one yet...)
   # work -------------------------------------------------------------------
   # adding a buffer of the distance specified above in the argument 'buff',
   # which will allow limited movement of this individual between years
   dataBuff <- st_buffer(dat, dist = buff)
     # get the data for the first year of sampling
     # probably use a for loop ...
     # i = year in inventory
     for (i in seq_along(inv)) {
       dataCurrent <- dat[dat$Year==inv[i],] # getting the data for year i, from
       # which we will select one individual to focus on

       # check if this the first year of sampling, or if this is the year
       # immediately following a gap in sampling.
       if (i == 1) { #is this the first year of sampling?
         #If this is the first year (i = 1), then all individuals get an NA
         # in the 'age' and 'recruit' columns.
         dataCurrent$age <- -9999
         dataCurrent$recruit <- -9999
         }


       ###AES### deal with track IDs.
       # get the buffered data for this year
       dataCurrentBuff <- dataBuff[dataBuff$Year==inv[i],]

       # j = individual in the dataset in that year
       for (j in seq_along(dataCurrent$Species)) {
         plant <- dataCurrent[j,] #getting data for one individual j, which we
         # will then compare to individuals in the next
         # year to see if it survives

         plantBuff <- dataCurrentBuff[j,] # select the buffered polygon for
         # the appropriate individual

          ## need to account for the fact that the buffer might be >1 ??? How ??

           for (k in 1:(dorm+1)) { # this loop compares the focal plant to
             # the quadrat data in the the next year. This loop will repeat 1 +
             # as many times as the 'dorm' argument
             # (compare to year t+1, ..., t+dorm).

             dataNext <- dat[dat$Year==inv[i + k],] #getting the data for year
             # i (year in which focal plant is in) + k (the number along a
             # sequence of values up to (1+ the value of dorm)). We will then
             # compare the buffered focal plant to this dataset to see if
             # there is any overlap

             overlapNoBuff<-  # see if there
             overlapNoBuffArea <- st_intersection(plant, dataNext)
             # is overlap of the un-buffered polygon with anything
             ## st_intersection gives the area of the focal individual that is overlapped... doesn't indicate which ones in year k+1 are doing the overlapping ... problem
             #use st_overlaps? or st_intersects?
             #st_intersects gives you the index of the overlapping polygons!
             mapview(plant) + mapview(dataNext[c(30,36,37),])
             # what is the percentage of overlap?
             overlapNoBuff
              if(nrow(overlapNoBuff)>0) { #if there IS overlap between unbuffered polygon and something in year k, then assign those w/ the overlap the same track ID

##########ALICE#########
              } else {

              }

                         ## if not, see if there is overlap of the buffered polygon with anything
             mapview(plant, col.regions = "green") + mapview(plantBuff, col.regions = "purple") + mapview(dataNext)
             ## if there is a solution in k=1, then STOP; but if there is not a solution in k =1, then continue the loop. How do we do this??
            ## is there a way to stop a for-loop if a certain condition is met??
           }

       }

       }
    }
   # output ---------------------------------------------------------------

 }
