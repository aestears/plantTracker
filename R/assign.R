#' tracks indivudal plants through time
#'
#' @return
#' @export
#'
#' @examples

# # required packages -------------------------------------------------------
# library(sf)
#
# # example input data ------------------------------------------------------
# # grasslandData (or exact same format), subset to a unique site, quad,
# # species
# # prepares the dataset to feed into the 'assign' function (the 'Assign'
# # function will do this ahead of time when the user calls it)
# sampleDat <- grasslandData[grasslandData$Site == "CO"
#                            & grasslandData$Quad == "unun_11"
#                            & grasslandData$Species == "Bouteloua gracilis",]
# # this should be a data.frame
#
# # get the appropriate grasslandInventory data for the "unun_11" quadrat,
# # to tell the 'assign' function when the quadrat was sampled
# sampleInv <- grasslandInventory[["unun_11"]]
# # this should be an integer vector
#
# # 'assign' function ---------------------------------------------------------
#
# function(){
#   # arguments
#   dat <- st_sf(sampleDat) # data in 'grasslandData' format, must be an sf object
#   inv <- sort(sampleInv) # integer vector of quadrat sampling years in
#   # sequential order
#   dorm <- 1 # dormancy allowed by the function
#   buff <- .05 # buffer of movement allowed from year to year, in meters
#   # work
#     # first, get the data for the first year of sampling
#     # probably use a for loop ...
#     # i = year in inventory
#     i = 1
#     temp1 <- dat[dat$Year==inv[i],]
#
#   # output
#
# }
