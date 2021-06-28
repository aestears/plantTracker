#' getNeighbors
#' @description
#'
#'
#' @param dat
#' @param radius
#' @param method
#' @param compType
#' @param focal
#' @param trackID
#' @param species
#' @param quad
#' @param year
#' @param site
#' @param geometry
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @details
#'
#' @return
#' @export
#'
#' @examples
getNeighbors <- function (dat, radius, method,
                          compType = 'allSpp',
                          focal = 'genet',
                          trackID = 'trackID',
                          species = "Species",
                          quad = "Quad",
                          year = "Year",
                          site = "Site",
                          geometry = "geometry",
                          ...
                          ) {
  #dat ## a df that has been assigned demographic data by the 'trackSpp'
  # function. If compType == 'allSpp', 'dat' must include ALL species present in
  # the quadrat
  #radius # a numeric value that has the desired radius around each individual
  # for which neighborhood area/count will be calculated.
  #method # must equal either 'count' or 'area'. If 'count', then the number of
  # genets within the specified radius will be counted. If 'area', then the area
  # within the specified radius that is occupied will be calculated.
  ###AES question: for count method--should I use ramets or genets??
  #focal ## an argument that indicates whether the focal indidvidual should be
  # each ramet ('ramet') or each genet ('genet'). Defualt is 'genet'
  #compType ## an argument that indicates whether you want to calculate only
  # conspecific local neighborhood (only look at individuals of the same
  # species -- 'oneSpp'), or heterospecific local neighborhood (look at
  # individuals of all species -- 'allSpp'). Default is 'allSpp'.
  #trackID ## the name of the column that contains the 'trackID'
  # (unique genet identifier) data

# argument checks ---------------------------------------------------------
## use the checkDat function to check the 'dat' function
checkedDat <- checkDat(dat = dat, species = species, site = site, quad = quad,
                       year = year, geometry = geometry, reformatDat = TRUE)

dat <- checkedDat$dat

## check trackID arg.
trackIDuserName <- paste0(trackID,"_USER")
trackIDdefaultName <- "trackID"
## change the name of the 'trackID' column to have the default column name
names(dat)[names(dat) == trackIDuserName] <- trackIDdefaultName

# work --------------------------------------------------------------------
## assign a unique index to each row to simplify the looping process
dat$index <- 1:nrow(dat)


## subset the 'dat' d.f so that it only has the required columns
# (and dstore the rest)
## data to 'store'
datStore <- st_drop_geometry(dat[,!names(dat) %in% c("Species", "Site", "Quad", "Year",
                                    "trackID", "geometry")])
## trimmed 'dat' to use in the function
dat <- dat[,names(dat) %in% c("Species", "Site", "Quad", "Year",
                              "trackID", "geometry", "index")]


## make an empty column in 'dat' to contain the output neighborhood data
dat$neighbors <- NA

for (i in unique(dat$Site)) { ## loop through each site
  for (j in unique(dat[dat$Site== i ,"Quad"]$Quad)) { ## loop through each
    # quadrat
    for (k in unique(dat[dat$Site == i & dat$Quad == j, "Year"]$Year)) { ## loop
      # through each year
      if (compType == 'oneSpp') { ## calculating conspecific neighborhood (need to
        # subset by species)
        for (l in unique(dat[dat$Site == i & dat$Quad == j & dat$Year == k,
                             "Species"]$Species)) {
          ## get the data for this site/quad/year/species combo
          datOneSp <- dat[dat$Site == i & dat$Quad == j &
                         dat$Year == k & dat$Species == l,]
            if (focal == 'ramet') { ## calculate the neighborhood for each
              # individual ramet
              for (n in 1:nrow(datOneSp)) { ## loop through each unique row
                ## make a buffer around the focal individual of the specified
                # 'radius'
                rametBuff <- suppressWarnings(st_buffer(x = datOneSp[n,],
                                                        dist = radius))
                ## subtract the focal ramet from the radius
                rametNeigh <- suppressWarnings(st_difference(x = rametBuff,
                                                             y = datOneSp[n,]))
                ## compare 'rametNeigh' (without the focal ramet) to the rest of
                # the ramets in this quad (of this spp.) (make sure to remove
                # the overlap with the focal individual)
                if ( method == 'count') {
                  ## get the number of genets that overlap w/ the focal radius
                  overlaps <- datOneSp[st_intersects(rametNeigh, datOneSp,
                                                     sparse = FALSE),]
                  ## remove focal individuals
                  overlaps <- overlaps[!overlaps$trackID %in%
                                         datOneSp[n,]$trackID,]

                  ## get the number of genets that are w/in the radius
                  count <- length(unique(overlaps$trackID))
                  ## put this value into 'dat' in the 'neighbors' column
                  dat[dat$index %in%
                        datOneSp[n,]$index,
                      ]$neighbors <- count
                } else if (method == 'area') {
                  ## get the overlapping polygons
                  overlapArea <- suppressWarnings(
                    st_intersection(rametNeigh, datOneSp)[st_is(st_intersection(
                      rametNeigh, datOneSp), type = c("MULTIPOLYGON", "POLYGON"))
                      ,])
                  ## get the area of the ramets that are w/in the radius
                  area <- sum(st_area(overlaps))
                  ## put this value into 'dat' in the 'neighbors' column for
                  # each ramet that belongs to the focal genet
                  dat[dat$index %in% datOneSp[n,]$index,
                      ]$neighbors <- area
                }
              }
            } else if (focal == 'genet') { ## calculate the neighborhood for
              # each genet
              ## loop through each unique trackID
              for (m in unique(datOneSp$trackID)) {
                ## make a buffer around the focal individual of the specified
                # 'radius'
                rametBuff <- suppressWarnings(st_buffer(
                  x = datOneSp[datOneSp$trackID == m,],
                  dist = radius))
                ## make into one large polygon
                rametBuff <- st_union(rametBuff)
                ## subtract the focal ramet from the radius
                rametNeigh <- suppressWarnings(
                  st_difference(x = rametBuff,
                          y = st_union(datOneSp[datOneSp$trackID == m,])))
                ## compare 'rametNeigh' (without the focal ramet) to the rest of
                # the ramets in this quad (of this spp.) (make sure to remove
                # the overlap with the focal individual)
                if ( method == 'count') {
                  ## get the number of genets that overlap w/ the focal radius
                  overlaps <- datOneSp[st_intersects(rametNeigh, datOneSp,
                                                   sparse = FALSE),]
                  ## remove focal individuals
                  overlaps <- overlaps[!overlaps$trackID %in%
                                         datOneSp[datOneSp$trackID == m,]$trackID,]

                  ## get the number of genets that are w/in the radius
                  count <- length(unique(overlaps$trackID))
                  ## put this value into 'dat' in the 'neighbors' column
                  dat[dat$index %in%
                        datOneSp[datOneSp$trackID == m,]$index,
                      ]$neighbors <- count
                } else if (method == 'area') {
                  ## get the overlapping polygons
                  overlapArea <- suppressWarnings(
                    st_intersection(rametNeigh, datOneSp)[st_is(st_intersection(
                      rametNeigh, datOneSp), type = c("MULTIPOLYGON", "POLYGON"))
                      ,])
                  ## get the area of the ramets that are w/in the radius
                  area <- sum(st_area(overlaps))
                  ## put this value into 'dat' in the 'neighbors' column for
                  # each ramet that belongs to the focal genet
                  dat[dat$index %in% datOneSp[datOneSp$trackID == m,]$index,
                      ]$neighbors <- area
                }
              }
            }
        } ## end of loop to get data by species
      } else if (compType == 'allSpp') { ## calculating heterospecific
        # neighborhood (don't need to subset by species)
        ## get the data for this site/quad/year combo
        datSpp<- dat[dat$Site == i & dat$Quad == j &
                          dat$Year == k ,]
        if (focal == 'ramet') { ## calculate the neighborhood for each
          # individual ramet
          for (n in 1:nrow(datSpp)) { ## loop through each unique row
            ## make a buffer around the focal individual of the specified
            # 'radius'
            rametBuff <- suppressWarnings(st_buffer(x = datSpp[n,],
                                                    dist = radius))
            ## subtract the focal ramet from the radius
            rametNeigh <- suppressWarnings(st_difference(x = rametBuff,
                                                         y = datSpp[n,]))
            ## compare 'rametNeigh' (without the focal ramet) to the rest of
            # the ramets in this quad (of this spp.) (make sure to remove
            # the overlap with the focal individual)
            if ( method == 'count') {
              ## get the number of genets that overlap w/ the focal radius
              overlaps <- datSpp[st_intersects(rametNeigh, datSpp,
                                               sparse = FALSE),]
              ## remove focal individuals
              overlaps <- overlaps[!overlaps$trackID %in%
                                     datSpp[n,]$trackID,]

              ## get the number of genets that are w/in the radius
              count <- length(unique(overlaps$trackID))
              ## put this value into 'dat' in the 'neighbors' column
              dat[dat$index %in%
                    datSpp[n,]$index,
                  ]$neighbors <- count
            } else if (method == 'area') {
              ## get the overlapping polygons
              overlapArea <- suppressWarnings(
                st_intersection(rametNeigh, datSpp)[st_is(st_intersection(
                  rametNeigh, datSpp), type = c("MULTIPOLYGON", "POLYGON"))
                  ,])
              ## get the area of the ramets that are w/in the radius
              area <- sum(st_area(overlapArea))
              ## put this value into 'dat' in the 'neighbors' column for
              # each ramet that belongs to the focal genet
              dat[dat$index %in% datSpp[n,]$index,
                  ]$neighbors <- area
            }
          }
        } else if (focal == 'genet') { ## calculate the neighborhood for
          # each genet
          ## loop through each unique trackID
          for (m in unique(datSpp$trackID)) {
            ## make a buffer around the focal individual of the specified
            # 'radius'
            rametBuff <- suppressWarnings(st_buffer(
              x = datSpp[datSpp$trackID == m,],
              dist = radius))
            ## make into one large polygon
            rametBuff <- st_union(rametBuff)
            ## subtract the focal ramet from the radius
            rametNeigh <- suppressWarnings(
              st_difference(x = rametBuff,
                            y = st_union(datSpp[datSpp$trackID == m,])))
            ## compare 'rametNeigh' (without the focal ramet) to the rest of
            # the ramets in this quad (of this spp.) (make sure to remove
            # the overlap with the focal individual)

            if ( method == 'count') {
              ## get the number of genets that overlap w/ the focal radius
              overlaps <- datSpp[st_intersects(rametNeigh, datSpp,
                                               sparse = FALSE),]
              ## remove focal individuals
              overlaps <- overlaps[!overlaps$trackID %in%
                                     datSpp[datSpp$trackID == m,]$trackID,]

              ## get the number of genets that are w/in the radius
              count <- length(unique(overlaps$trackID))
              ## put this value into 'dat' in the 'neighbors' column
              dat[dat$index %in%
                    datSpp[datSpp$trackID == m,]$index,
                  ]$neighbors <- count
            } else if (method == 'area') {
              ## get the overlapping polygons
              overlapArea <- suppressWarnings(
                st_intersection(rametNeigh, datSpp)[st_is(st_intersection(
                  rametNeigh, datSpp), type = c("MULTIPOLYGON", "POLYGON"))
                  ,])
              ## get the area of the ramets that are w/in the radius
              area <- sum(st_area(overlapArea))
              ## put this value into 'dat' in the 'neighbors' column for
              # each ramet that belongs to the focal genet
              dat[dat$index %in% datSpp[datSpp$trackID == m,]$index,
                  ]$neighbors <- area
            }
          }
        }
      } ## end of loop to get heterospecific dataset (all species present)
    } ## end of loop to get data by year
    } ## end of loop to get data by quadrat
  } ## end of loop to get data by site

# output ------------------------------------------------------------------
## prepare the output
## change names of 'dat' back to the user-input values
## get the user-provided column names
userColNames <- c(checkedDat$userColNames[c(1:4)],
                  trackIDuserName,
                  checkedDat$userColNames[c(5)])

## rejoin the trackSppOut d.f with the 'extra' data stored in 'datStore'
dat <- merge(dat, datStore, by = "index")
## remove the 'indexStore' value
dat <- dat[,names(dat) != "index"]

## re-name the appropriate columns in 'dat' data.frame with the
# user-provided names of 'dat'
## from above, user-provided names are stored in 'userColNames'
## make a vector of default column names
defaultNames <- c("Species", "Site", "Quad", "Year",  "geometry", "trackID")

## reset the names for the columns that we changed to 'default' values
names(dat)[which(names(dat) %in% defaultNames)] <- userColNames

## remove the '_USER' from the 'extra' column names
names(dat) <- gsub(names(dat),
                           pattern = "_USER", replacement = "")

## return
return(dat)
} ## end of 'getNeighbors()'


# testing -----------------------------------------------------------------

dat <- grasslandData[grasslandData$Site == "CO" &
                    grasslandData$Quad == "ungz_5a",]
datIDs<- trackSpp(dat = dat, inv = grasslandInventory, dorm = 1, buff = 0.05,
                  buffGenet = 0.005,
                  clonal = data.frame("Species" = c("Bouteloua gracilis",
                                                    "Agropyron smithii",
                                                    "Sphaeralcea coccinea"),
                  "clonal" = c(1,1,0)))
names(datIDs)[c(1,6)] <- c("speciesName", "uniqueID")

getNeighbors(dat = datIDs, radius = .15, method = "area",
          compType = 'allSpp',
          focal = 'genet',
          trackID = 'uniqueID',
          species = "speciesName",
          quad = "Quad",
          year = "Year",
          site = "Site",
          geometry = "geometry")

#### load packages ####
#require(tidyverse) #v1.3.0
#require(sf) #v0.9-7
#require(mapview) #v2.9.0
#require(lwgeom) #v0.2-5
#
# #make a vector for years in the dataset
# year <- sort(unique(poly$year))
# #make a vector for quadrats in the dataset
# quad <- unique(poly$quad)
# #make a vector for species in the dataset
# species <- unique(poly$species)
#
# #### get the data for SP_ID from tracking files and merge this with the poly dataset (for individuals with only one observation--haven't been clumped at genet scale) ####
#
# #load all tracking data files with the following loop
#
# #set the working directory to the location of the folder "PolygonTrackingResults"
# trackingDatWD <- "/Users/Alice/Dropbox/Grad School/Research/Trait Project/CO_sgs Analysis/trackingData/InputData/PolygonTrackingResults" #change to appropriate working directory
# setwd(trackingDatWD)
#
# for(k in 1:length(species)) { #loop through all species
#   trackTemp <- read.csv(paste("./",
#                               str_to_upper(paste(str_sub(species[k],1,3),
#                                                  str_sub(species[k],str_locate(species[k]," ")[1]+1,str_locate(species[k]," ")[1]+3),sep = ""))
#                               ,"_buf5_dorm1.csv", sep = ""))
#   if (k == 1) {
#     trackSP <- trackTemp
#   } else {
#     trackSP <- rbind(trackSP, trackTemp)
#   }
# }
# #trackSP data file has all raw tracking results
#
# ## aggregate this data.frame by quad, year, species, and trackID to get the trackIDs of individuals that are singletons (i.e. not mapped as multiple polygons)
# #will use this data to identify the SP_ID's that are 'good' in the trackSP data.set
# trackGOOD <- aggregate(trackSP, by = list(trackSP$quad, trackSP$year, trackSP$Species, trackSP$trackID), FUN = length)
# # drop all IDs that have more than one individual (would have >1 SP_ID)
# trackGOOD <- filter(trackGOOD,quad==1)
# #add a column so we know the individual record is 'good'
# trackGOOD$need <- "need"
# #fix the names of the data.frame
# trackGOOD <- trackGOOD[,c("Group.1", "Group.2", "Group.3", "Group.4", "need")]
# names(trackGOOD) <- c("quad","year","Species","trackID","need")
# #merge with trackSP dataset to ID 'good' records (denoted by "need" in the need column)
# trackSP <- left_join(trackSP, trackGOOD, by = c("quad", "year", "Species", "trackID"))
# #filter the trackSP data for those that we "need"
# trackSP <- filter(trackSP, need=="need")
# #merge with poly dataset to get SP_IDs
# poly<- left_join(poly, trackSP[,c("quad","year","SP_ID","Species","trackID")], by = c("quad", "year_t"="year", "species"="Species", "trackID"))
#
# #### calculate nearest neighbor for polygons with only one SP_ID ####
# shpWD <- "/Users/Alice/Dropbox/Grad School/Research/Trait Project/Data/Adler Dowloaded Datasets/Adler_CO_Downloaded Data/CO_shapefiles" #change to your file that contains the CO shapefiles
# setwd(shpWD)
#
# #make a bounding box that is the shape and size of the quadrat
# baseShape <- st_read(dsn = "gzgz_5a", layer = "poly_gzgz_5a_1997")
# box <- st_as_sfc(st_bbox(baseShape)) #make a box that is the size of a quadrat (1mx1m)
# boxBuffer <- st_buffer(box, dist = -0.05)#make a 5cm buffer inside the box
# boxBuffer <- st_difference(box, boxBuffer)
#
# ## Calculate nearest neighbor for 5cm radius
# for(j in 1:length(quad)) { #loop throught the quadrats
#   temp1 <- poly[poly$quad==quad[j],]
#   for (i in 1:length(year)) { #loop through the years
#     temp2 <- temp1[temp1$year==year[i],]
#     for (k in 1:length(species)) { #loop through the species
#       temp3 <- temp2[temp2$species == species[k],]
#       if (nrow(temp3)>0) {
#         #now find the appropriate shapefile
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("poly_",quad[j],"_",year[i], sep = ""), stringsAsFactors = FALSE)
#         shape$SP_ID <- as.numeric(shape$SP_ID)
#         shapeFocal <- drop_na(shape[shape$Species == as.character(species[k]),],"Species") #all of the polygons for the target species
#         shapeFocal <- st_make_valid(shapeFocal)
#         if(nrow(shapeFocal)>0){
#           for (m in 1:nrow(shapeFocal)) {
#             buffer <- st_buffer(shapeFocal[m,],.05)
#             area <- st_intersection(buffer, shapeFocal)
#             area2 <- st_difference(area, shapeFocal[m,])
#             area3 <- sum(st_area(area2))
#             poly[poly$quad==quad[j] &
#                    poly$year==year[i] &
#                    poly$SP_ID == shapeFocal[m,]$SP_ID &
#                    is.na(poly$SP_ID)==FALSE ,
#                  "neighbor_area_5"] <- area3
#           }
#         }
#       }
#     }
#   }
# }
#
# #for 10cm nearest neighbor radius
# for(j in 1:length(quad)) {
#   temp1 <- poly[poly$quad==quad[j],]
#   for (i in 1:length(year)) {
#     temp2 <- temp1[temp1$year==year[i],]
#     for (k in 1:length(species)) {
#       temp3 <- temp2[temp2$species == species[k],]
#       if (nrow(temp3)>0) {
#         #now find the appropriate shapefile
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("poly_",quad[j],"_",year[i], sep = ""), stringsAsFactors = FALSE)
#         shape$SP_ID <- as.numeric(shape$SP_ID)
#         shapeFocal <- drop_na(shape[shape$Species == as.character(species[k]),],"Species") #all of the polygons
#         shapeFocal <- st_make_valid(shapeFocal)
#         if(nrow(shapeFocal)>0){
#           for (m in 1:nrow(shapeFocal)) {
#             buffer <- st_buffer(shapeFocal[m,],.1)
#             area <- st_intersection(buffer, shapeFocal)
#             area2 <- st_difference(area, shapeFocal[m,])
#             area3 <- sum(st_area(area2))
#             poly[poly$quad==quad[j] &
#                    poly$year==year[i] &
#                    poly$SP_ID == shapeFocal[m,]$SP_ID &
#                    is.na(poly$SP_ID)==FALSE ,
#                  "neighbor_area_10"] <- area3
#           }
#         }
#       }
#     }
#   }
# }
#
# #for 15cm nearest neighbor radius
# for(j in 1:length(quad)) {
#   temp1 <- poly[poly$quad==quad[j],]
#   for (i in 1:length(year)) {
#     temp2 <- temp1[temp1$year==year[i],]
#     for (k in 1:length(species)) {
#       temp3 <- temp2[temp2$species == species[k],]
#       if (nrow(temp3)>0) {
#         #now find the appropriate shapefile
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("poly_",quad[j],"_",year[i], sep = ""), stringsAsFactors = FALSE)
#         shape$SP_ID <- as.numeric(shape$SP_ID)
#         shapeFocal <- drop_na(shape[shape$Species == as.character(species[k]),],"Species") #all of the polygons
#         shapeFocal <- st_make_valid(shapeFocal)
#         if(nrow(shapeFocal)>0){
#           for (m in 1:nrow(shapeFocal)) {
#             buffer <- st_buffer(shapeFocal[m,],.15)
#             area <- st_intersection(buffer, shapeFocal)
#             area2 <- st_difference(area, shapeFocal[m,])
#             area3 <- sum(st_area(area2))
#             poly[poly$quad==quad[j] &
#                    poly$year==year[i] &
#                    poly$SP_ID == shapeFocal[m,]$SP_ID &
#                    is.na(poly$SP_ID)==FALSE ,
#                  "neighbor_area_15"] <- area3
#           }
#         }
#       }
#     }
#   }
# }
#
# #for 20cm nearest neighbor radius
# for(j in 1:length(quad)) {
#   temp1 <- poly[poly$quad==quad[j],]
#   for (i in 1:length(year)) {
#     temp2 <- temp1[temp1$year==year[i],]
#     for (k in 1:length(species)) {
#       temp3 <- temp2[temp2$species == species[k],]
#       if (nrow(temp3)>0) {
#         #now find the appropriate shapefile
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("poly_",quad[j],"_",year[i], sep = ""), stringsAsFactors = FALSE)
#         shape$SP_ID <- as.numeric(shape$SP_ID)
#         shapeFocal <- drop_na(shape[shape$Species == as.character(species[k]),],"Species") #all of the polygons
#         shapeFocal <- st_make_valid(shapeFocal)
#         if(nrow(shapeFocal)>0){
#           for (m in 1:nrow(shapeFocal)) {
#             buffer <- st_buffer(shapeFocal[m,],.2)
#             area <- st_intersection(buffer, shapeFocal)
#             area2 <- st_difference(area, shapeFocal[m,])
#             area3 <- sum(st_area(area2))
#             poly[poly$quad==quad[j] &
#                    poly$year==year[i] &
#                    poly$SP_ID == shapeFocal[m,]$SP_ID &
#                    is.na(poly$SP_ID)==FALSE ,
#                  "neighbor_area_20"] <- area3
#           }
#         }
#       }
#     }
#   }
# }
#
# #### calculate nearest neighbor area for individuals with multiple SP ID's (mapped as multiple polygons) ####
# #set up a buffer box
# baseShape <- st_read(dsn = "gzgz_5a", layer = "poly_gzgz_5a_1997")
# box <- st_as_sfc(st_bbox(baseShape)) #make a box that is the size of a quadrat (1mx1m)
# boxBuffer <- st_buffer(box, dist = -0.05)#make a 5cm buffer inside the box
# boxBuffer <- st_difference(box, boxBuffer)
#
# #reminders:
# # trackingDatWD is an object that stores the path of the "PolygonTrackingResults" directory
# # shpWD is an object that stores the path of the "CO_shapefiles" directory
#
# #for 5cm radius
# for (k in 1:length(species)){ ##loop through species
#   setwd(trackingDatWD)
#   tempSP <- read.csv(paste(str_to_upper(paste(str_sub(species[k],1,3),
#                                               str_sub(species[k],str_locate(species[k]," ")[1]+1,str_locate(species[k]," ")[1]+3),sep = ""))
#                            ,"_buf5_dorm1.csv", sep = ""))
#   tempSP$area <- round(tempSP$area, 3)
#   #load datafile that has nearest neighbor data (data after genets have been grouped)
#   #filter for species
#   neSP <- poly[poly$species==species[k],]
#   neSP$area <- round(neSP$area_t, 3)
#   neSP <- neSP[is.na(neSP$neighbor_area_5)==TRUE,]
#   for(j in 1:length(quad)) {
#     ##loop through quadrats
#     quadNeSP <- neSP[neSP$quad==paste(quad[j]),]
#     quadtempSP <- tempSP[tempSP$quad==paste(quad[j]),]
#     for(i in 1:length(year)) {
#       ##loop through years
#       yearNeSP <- quadNeSP[quadNeSP$year==year[i],]
#       yeartempSP <- quadtempSP[quadtempSP$year==year[i],]
#       #get the SP_IDs of the individuals that do not have nearest neighbor info
#       # (the SP_ID is in the 'tracking' data file, but the species ID is in the 'ne' data file)
#       yearSP <- yeartempSP[yeartempSP$trackID %in% unique(yearNeSP$trackID),]
#       if(nrow(yearSP)>length(unique(yearSP$trackID))) {
#         #load appropriate shape file
#         setwd(shpWD)
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("poly_",quad[j],"_",year[i], sep=""))
#         #get focal species from shape file
#         shapeFocal <- drop_na(shape[shape$Species == paste(species[k]),],"Species")
#         shapeFocal <- st_make_valid(shapeFocal)
#         if(nrow(shape)>0) {
#           for(m in 1:length(unique(yearSP$trackID))) {
#             ##loop through the groups of individuals that have the same track ID
#             #identify all of the indviduals that have the same track ID
#             oneSP <- yearSP[yearSP$trackID==unique(yearSP$trackID)[m],]
#             #identify the individuals that have the same track ID (by their SP_ID in the shape file)
#             if(nrow(oneSP)>1) {
#               oneShape <- shapeFocal[shapeFocal$SP_ID %in% (oneSP$SP_ID) ,]
#               # draw buffers around all of the shapes
#               buffer <- st_buffer(oneShape,.05)
#               # remove the overlap and combine into one shapefile
#               buffer2 <- st_union(buffer)
#               finalBuffer<- st_difference(buffer2, st_union(oneShape))
#               # calculate overlap with other 'individuals'
#               area <- st_intersection(st_make_valid(finalBuffer), st_make_valid(shapeFocal))
#               area2 <- sum(st_area(area))
#               #insert this value into the appropriate place in the dataframe for nearest neighbor values
#               poly[poly$species==species[k] &
#                      poly$quad==quad[j] &
#                      poly$year==year[i] &
#                      poly$trackID==unique(yearSP$trackID)[m],
#                    "neighbor_area_5"]<- area2
#             }
#           }
#         }
#       }
#     }
#   }
# }
#
# #for 10cm radius
# for (k in 1:length(species)){
#   setwd(trackingDatWD)
#   tempSP <- read.csv(paste(
#     str_to_upper(paste(str_sub(species[k],1,3),
#                        str_sub(species[k],str_locate(species[k]," ")[1]+1,str_locate(species[k]," ")[1]+3),sep = ""))
#     ,"_buf5_dorm1.csv", sep = ""))
#   tempSP$area <- round(tempSP$area, 3)
#   #load datafile that has nearest neighbor data (data after genets have been grouped)
#   #filter for species
#   neSP <- poly[poly$species==species[k],]
#   neSP$area <- round(neSP$area_t, 3)
#   neSP <- neSP[is.na(neSP$neighbor_area_10)==TRUE,]
#   for(j in 1:length(quad)) {
#     ##loop through quadrats
#     quadNeSP <- neSP[neSP$quad==paste(quad[j]),]
#     quadtempSP <- tempSP[tempSP$quad==paste(quad[j]),]
#     for(i in 1:length(year)) {
#       ##loop through years
#       yearNeSP <- quadNeSP[quadNeSP$year==year[i],]
#       yeartempSP <- quadtempSP[quadtempSP$year==year[i],]
#       #get the SP_IDs of the individuals that do not have nearest neighbor info
#       # (the SP_ID is in the 'tracking' data file, but the species ID is in the 'ne' data file)
#       yearSP <- yeartempSP[yeartempSP$trackID %in% unique(yearNeSP$trackID),]
#       if(nrow(yearSP)>length(unique(yearSP$trackID))) {
#         #load appropriate shape file
#         setwd(shpWD)
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("poly_",quad[j],"_",year[i], sep=""))
#         #get focal species from shape file
#         shapeFocal <- drop_na(shape[shape$Species == paste(species[k]),],"Species")
#         shapeFocal <- st_make_valid(shapeFocal)
#         if(nrow(shape)>0) {
#           for(m in 1:length(unique(yearSP$trackID))) {
#             ##loop through the groups of individuals that have the same track ID
#             #identify all of the indviduals that have the same track ID
#             oneSP <- yearSP[yearSP$trackID==unique(yearSP$trackID)[m],]
#             #identify the individuals that have the same track ID (by their SP_ID in the shape file)
#             if(nrow(oneSP)>1) {
#               oneShape <- shapeFocal[shapeFocal$SP_ID %in% (oneSP$SP_ID) ,]
#               # draw buffers around all of the shapes
#               buffer <- st_buffer(oneShape,.1)
#               # remove the overlap and combine into one shapefile
#               buffer2 <- st_union(buffer)
#               finalBuffer<- st_difference(buffer2, st_union(oneShape))
#               # calculate overlap with other 'individuals'
#               area <- st_intersection(st_make_valid(finalBuffer), st_make_valid(shapeFocal))
#               area2 <- sum(st_area(area))
#               #insert this value into the appropriate place in the dataframe for nearest neighbor values
#               poly[poly$species==species[k] &
#                      poly$quad==quad[j] &
#                      poly$year==year[i] &
#                      poly$trackID==unique(yearSP$trackID)[m],
#                    "neighbor_area_10"]<- area2
#             }
#           }
#         }
#       }
#     }
#   }
# }
#
# #for 15cm radius
# for (k in 1:length(species)){
#   setwd(trackingDatWD)
#   tempSP <- read.csv(paste(str_to_upper(paste(str_sub(species[k],1,3),
#                                               str_sub(species[k],str_locate(species[k]," ")[1]+1,str_locate(species[k]," ")[1]+3),sep = ""))
#                            ,"_buf5_dorm1.csv", sep = ""))
#   tempSP$area <- round(tempSP$area, 3)
#   #load datafile that has nearest neighbor data (data after genets have been grouped)
#   #filter for species
#   neSP <- poly[poly$species==species[k],]
#   neSP$area <- round(neSP$area_t, 3)
#   neSP <- neSP[is.na(neSP$neighbor_area_15)==TRUE,]
#   for(j in 1:length(quad)) {
#     ##loop through quadrats
#     quadNeSP <- neSP[neSP$quad==paste(quad[j]),]
#     quadtempSP <- tempSP[tempSP$quad==paste(quad[j]),]
#     for(i in 1:length(year)) {
#       ##loop through years
#       yearNeSP <- quadNeSP[quadNeSP$year==year[i],]
#       yeartempSP <- quadtempSP[quadtempSP$year==year[i],]
#       #get the SP_IDs of the individuals that do not have nearest neighbor info
#       # (the SP_ID is in the 'tracking' data file, but the species ID is in the 'ne' data file)
#       yearSP <- yeartempSP[yeartempSP$trackID %in% unique(yearNeSP$trackID),]
#       if(nrow(yearSP)>length(unique(yearSP$trackID))) {
#         #load appropriate shape file
#         setwd(shpWD)
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("poly_",quad[j],"_",year[i], sep=""))
#         #get focal species from shape file
#         shapeFocal <- drop_na(shape[shape$Species == paste(species[k]),],"Species")
#         shapeFocal <- st_make_valid(shapeFocal)
#         if(nrow(shape)>0) {
#           for(m in 1:length(unique(yearSP$trackID))) {
#             ##loop through the groups of individuals that have the same track ID
#             #identify all of the indviduals that have the same track ID
#             oneSP <- yearSP[yearSP$trackID==unique(yearSP$trackID)[m],]
#             #identify the individuals that have the same track ID (by their SP_ID in the shape file)
#             if(nrow(oneSP)>1) {
#               oneShape <- shapeFocal[shapeFocal$SP_ID %in% (oneSP$SP_ID) ,]
#               # draw buffers around all of the shapes
#               buffer <- st_buffer(oneShape,.15)
#               # remove the overlap and combine into one shapefile
#               buffer2 <- st_union(buffer)
#               finalBuffer<- st_difference(buffer2, st_union(oneShape))
#               # calculate overlap with other 'individuals'
#               area <- st_intersection(st_make_valid(finalBuffer), st_make_valid(shapeFocal))
#               area2 <- sum(st_area(area))
#               #insert this value into the appropriate place in the dataframe for nearest neighbor values
#               poly[poly$species==species[k] &
#                      poly$quad==quad[j] &
#                      poly$year==year[i] &
#                      poly$trackID==unique(yearSP$trackID)[m],
#                    "neighbor_area_15"]<- area2
#             }
#           }
#         }
#       }
#     }
#   }
# }
#
# #for 20cm radius
# for (k in 1:length(species)){
#   setwd(trackingDatWD)
#   tempSP <- read.csv(paste(str_to_upper(paste(str_sub(species[k],1,3),
#                                               str_sub(species[k],str_locate(species[k]," ")[1]+1,str_locate(species[k]," ")[1]+3),sep = ""))
#                            ,"_buf5_dorm1.csv", sep = ""))
#   tempSP$area <- round(tempSP$area, 3)
#   #load datafile that has nearest neighbor data (data after genets have been grouped)
#   #filter for species
#   neSP <- poly[poly$species==species[k],]
#   neSP$area <- round(neSP$area_t, 3)
#   neSP <- neSP[is.na(neSP$neighbor_area_20)==TRUE,]
#   for(j in 1:length(quad)) {
#     ##loop through quadrats
#     quadNeSP <- neSP[neSP$quad==paste(quad[j]),]
#     quadtempSP <- tempSP[tempSP$quad==paste(quad[j]),]
#     for(i in 1:length(year)) {
#       ##loop through years
#       yearNeSP <- quadNeSP[quadNeSP$year==year[i],]
#       yeartempSP <- quadtempSP[quadtempSP$year==year[i],]
#       #get the SP_IDs of the individuals that do not have nearest neighbor info
#       # (the SP_ID is in the 'tracking' data file, but the species ID is in the 'ne' data file)
#       yearSP <- yeartempSP[yeartempSP$trackID %in% unique(yearNeSP$trackID),]
#       if(nrow(yearSP)>length(unique(yearSP$trackID))) {
#         #load appropriate shape file
#         setwd(shpWD)
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("poly_",quad[j],"_",year[i], sep=""))
#         #get focal species from shape file
#         shapeFocal <- drop_na(shape[shape$Species == paste(species[k]),],"Species")
#         shapeFocal <- st_make_valid(shapeFocal)
#         if(nrow(shape)>0) {
#           for(m in 1:length(unique(yearSP$trackID))) {
#             ##loop through the groups of individuals that have the same track ID
#             #identify all of the indviduals that have the same track ID
#             oneSP <- yearSP[yearSP$trackID==unique(yearSP$trackID)[m],]
#             #identify the individuals that have the same track ID (by their SP_ID in the shape file)
#             if(nrow(oneSP)>1) {
#               oneShape <- shapeFocal[shapeFocal$SP_ID %in% (oneSP$SP_ID) ,]
#               oneShape <- st_union(oneShape, by_feature = FALSE)
#               # draw buffers around all of the shapes
#               buffer <- st_buffer(oneShape,.20)
#               # remove the overlap and combine into one shapefile
#               buffer2 <- st_union(buffer)
#               finalBuffer<- st_difference(buffer2, oneShape)
#               # calculate overlap with other 'individuals'
#               area <- st_intersection(st_make_valid(finalBuffer), st_make_valid(shapeFocal))
#               area2 <- sum(st_area(area))
#               #insert this value into the appropriate place in the dataframe for nearest neighbor values
#               poly[poly$species==species[k] &
#                      poly$quad==quad[j] &
#                      poly$year==year[i] &
#                      poly$trackID==unique(yearSP$trackID)[m],
#                    "neighbor_area_20"]<- area2
#             }
#           }
#         }
#       }
#     }
#   }
# }
#
#
# # #write  to csv file
# # write.csv(poly,
# #           "/Users/Alice/Dropbox/Grad School/Research/Trait Project/Data/CO Analysis Data Files/Intermediate Analysis Files/polygon_demo_2_12_18.csv",
# #           row.names = FALSE)
#
#
# #### Calculate Nearest Neighbor for Points Dataset ####
# pointsWD <- "/Users/Alice/Dropbox/Grad School/Research/Trait Project/CO_sgs Analysis/trackingData/SurvivalData" #change the location of your 'point_species_survD.csv" file
# setwd(pointsWD)
#
# #read in point survival data
# points <- read.csv("point_species_survD.csv", stringsAsFactors = FALSE)
# #add a column for 'site'
# points$Site <- "CO"
#
# #create an empty column for nearest neighbor density
# points$neighbors_5 <- NA
# points$neighbors_10 <- NA
# points$neighbors_15 <- NA
# points$neighbors_20 <- NA
#
# #make a vector for years in the dataset
# year <- sort(unique(points$year))
# quad <- unique(points$quad) #make a vector of quads in the dataset
# species <- unique(points$species) #make a vector of species in the dataset
#
# setwd(shpWD) #change WD to the path for the CO_shapefiles folder (which is the same as "shpWD", defined above)
#
# #calculate nearest neighbor density for 5cm radius
# for(j in 1:length(quad)) {
#   temp1 <- points[points$quad==quad[j],]
#   for (i in 1:length(year)) {
#     temp2 <- temp1[temp1$year==year[i],]
#     for (k in 1:length(species)) {
#       temp3 <- temp2[temp2$species == species[k],]
#       if (nrow(temp3)>0) {
#         shape <- st_read(dsn = paste(quad[j]), layer = paste("pnt_",quad[j],"_",year[i], sep = ""))
#         shape <- shape[shape$Species == as.character(species[k]),]
#         for (m in 1:nrow(shape)) {
#           buffer <- st_buffer(shape[m,],.05)
#           count <- ifelse(length(unlist(st_intersects(buffer, shape)))==0,
#                           yes = (0), no = (length(unlist(st_intersects(buffer, shape)))-1))
#           points[points$x < (as.data.frame(shape[m,"coords_x1"])[1,1]+.0001) &
#                    points$x > (as.data.frame(shape[m,"coords_x1"])[1,1]-.0001) &
#                    points$quad == quad[j] &
#                    points$year == year[i] &
#                    points$species == species[k], "neighbors_5"] <- count
#         }
#       }
#     }
#   }
# }
#
# # remove individuals from the 'points' dataset that are <5cm from the edges of a plot (won't include in the analysis, since we cannot accurately estimate nearest neighbor density)
# points$edgeAS <- NA
# points[points$x<0.05 | points$x>0.95 | points$y<0.05 | points$y>0.95,"edgeAS"] <- TRUE
# points[points$x>=0.05 & points$x<=0.95 & points$y>=0.05 & points$y<=0.95,"edgeAS"] <- FALSE
#
# #### remove files that aren't needed for further analysis ####
# rm(list = ls()[!(ls() %in% c('points','poly'))])
#
# #### for next script, need 'points' and 'poly' data.frames ####
# #save as an .RData file
# path = "/Users/Alice/Dropbox/Grad School/Research/Trait Project/CO_sgs Analysis/CO-Sgs-paper/scripts" #location where you'll put the environment data file
# setwd(path)
# save.image('script0_output.RData')
