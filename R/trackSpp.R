#' Tracks genets through time for multiple species and sites
#'
#' @description This function tracks individual plants in a mapped quadrat
#' through time to generate a demographic data dataset that includes survival
#' and growth.
#'
#' @details
#' This function is a wrapper function that applies the \code{\link{assign}}
#' function accross multiple species, quadrats, or sites. For each species and
#' quadrat, [trackSpp()] loads a spatially referenced data.frame ('dat'), and
#' then uses the \code{\link{groupByGenet}} function to assign genetIDs to
#' polygons (if 'clonal' = 1) such that polygons that form the same genet have
#' the same genetID. A buffer of a distance defined by 'buff' is applied around
#' each genet polygon. Then, the spatial data for each genet from the current
#' year (year *t*) is compared to individuals in the next year (year *t+1*).
#' Then [trackSpp()] calculates the amount of overlapping area between polygons
#' of each year *t* genet and polygons of each year *t+1* genet (using
#' \code{\link[sf]{st_intersection}}). If there is unambiguous overlap between a
#' 'parent' genet from year *t* and a 'child' genet from year *t+1*, then that
#' 'child' gets the same identifying trackID as the parent. If there is a 'tie,'
#' where more than one parent overlaps the same child or more than one child
#' overlaps the same parent, the parent-child pair with the greatest amount of
#' overlap receives the same trackID. Polygons in year *t+1* that do not have a
#' parent are given new trackIDs and are identified as new recruits. If dormancy
#' is not allowed, then polygons in year *t* that do not have child polygons get
#' a '0' in the 'survival' column. If dormancy is allowed, parent polygons
#' without child polygons are stored as 'ghosts' and are then compared to data
#' from year *t+1+i* to find potential child polygons, where *i*='dorm'
#' argument.
#'
#' @param dat An sf data.frame of the same format as
#' \code{\link{grasslandData}}. It must have columns that contains a unique
#' identification for each research site (default name is "Site"), species name
#' (default name is "Species"), quadrat identifier (default name is "Quad"),
#' year of data collection (default name is "Year") and an s.f 'geometry' column
#' that contains a polygon or multipolygon data type for each
#' individual observation.
#' This function will add columns called "basalArea", "trackID",
#' "age", "size_tplus1", "recruit" and "survives_tplus1", so 'dat' should not
#' contain columns with these names.
#' @param inv A named list. The name of each element of the list is a quadrat
#' name in 'dat', and the contents of that list element is a numeric vector of
#' all of the years in which that quadrat (or other unique spatial area) was
#' sampled.
#' @param dorm A numeric vector of length 1, indicating the number of years this
#' species is allowed to go dormant, i.e. be absent from the map but be
#' considered the same individual when it reappears. This must be an integer
#' greater than or equal to 0. OR dorm can be a data.frame with
#' the columns "Species" and "dorm". This data.frame must have a row for each
#' unique species present in 'dat', with species name as a character string in
#' the "Species" column, and a numeric value greater than or equal to 0 in the
#' 'dorm' column that indicates the number of years this species is allowed to
#' go dormant.
#' @param buff A numeric vector of length 1 that is greater than or equal to
#' zero, indicating how far (in the same units as spatial values in 'dat') a
#' polygon can move from year \code{i} to year \code{i}+1 and still be
#' considered the same individual. OR buff can be a data.frame with
#' the columns "Species" and "buff". This data.frame must have a row for each
#' unique species present in 'dat', with species name as a character string in
#' the "Species" column, and a numeric value in the 'buff' column that indicates
#' how far a polygon can move from year to year and still be considered the same
#' genet.
#' @param buffGenet A numeric vector of length 1 that is greater than or equal
#' to zero, indicating how close (in the same units as spatial values in 'dat')
#' polygons must be to one another in the same year to be grouped as a genet
#' (if 'clonal' argument = 1). OR buffGenet can be a data.frame with
#' the columns "Species" and "buffGenet". This data.frame must have a row for
#' each unique species present in 'dat', with species name as a character string
#' in the "Species" column, and a numeric value greater than or equal to 0 in
#' the 'buffGenet' column that indicates how close polygons of that species must
#' be to one another to be considered the same genet. This argument is passed to
#' the \code{\link{groupByGenet}} function, which is used inside
#' \code{\link{assign}} function.
#' @param clonal A numeric Boolean vector of length 1, indicating whether a
#' species is allowed to be clonal or not (i.e. if multiple polygons (ramets)
#' can be grouped as one individual (genet)). OR clonal can be a data.frame with
#' the columns "Species" and "clonal". This data.frame must have a row for each
#' unique species present in 'dat', with species name as a character string in
#' the "Species" column, and a Boolean value in the 'clonal' column indicating
#' whether that species is allowed to be clonal (1) or not (0).
#' @param species An optional single character string argument that indicates
#' the name of the column in 'dat' that contains species name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Species".
#' @param site An optional single character string argument that indicates
#' the name of the column in 'dat' that contains site name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Site".
#' @param quad An optional single character string argument that indicates
#' the name of the column in 'dat' that contains quadrat name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Quad".
#' @param year An optional single character string argument that indicates
#' the name of the column in 'dat' that contains data for year of sampling. It
#' is unnecessary to include a value for this argument if the column name is
#' "Year".
#' @param geometry An optional single character string argument that indicates
#' the name of the column in 'dat' that contains sf geometry data. It is
#' unnecessary to include a value for this argument if the column name is
#' "geometry".
#'
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return An sf data.frame with the same columns as 'dat,' but with the
#' following additional columns:
#'
#' \item{trackID}{A unique value for each individual genet, consisting of the
#' 6-letter species code, the year in which this individual was recruited, and a
#' unique index number, all separated by a "_".}
#' \item{age}{An integer indicating the age of this individual in year *t*.
#' Values of NA indicate that an accurate age cannot be calculated because this
#' individual was observed either in the first year of sampling or in a year
#' following a gap in sampling, so the exact year of recruitment is not known.}
#' \item{size_tplus1}{The  size of this genet in year *t+1*, in the same units
#' as the 'area' column in 'dat'.}
#' \item{recruit}{A Boolean integer indicating whether this individual is a new
#' recruit in year *t* (1), or existed in a previous year (0). Values of NA
#' indicate that this individual was observed either in the first year of
#' sampling or in a year following a gap in sampling, so it is not possible to
#' accurately determine whether or not it is a new recruit in year *t*.}
#' \item{survives_plus1}{A Boolean integer indicating whether this individual
#' survived (1), or died (0) in year *t+1*.}
#' \item{genetArea}{The size of this entire genet in year *t*, in the same units
#' as the 'area' column in 'dat.' If the 'clonal' argument =0, then this number
#' will be identical to the 'area' column in 'dat'. }
#'
#' @seealso [assign()], which is used inside the [trackSpp()] function.
#' [trackSpp()] applies the [assign()] function over multiple species and
#' quadrats. The [assign()] function uses the [groupByGenet()] function to group
#' ramets into genets (if 'clonal' argument = 1).
#'
#' @examples
#' dat <- grasslandData[grasslandData$Site %in% c("CO", "AZ"),]
#' names(dat)[1] <- "speciesName"
#' inv <- grasslandInventory[unique(dat$Quad)]
#' out_dat <- trackSpp(dat = dat,
#'  inv = inv,
#'  dorm = 1,
#'  buff = .05,
#'  buffGenet = 0.005,
#'  clonal = data.frame("Species" = unique(dat$speciesName),
#'  "clonal" = c(1,1,0,0,1,1,0,0)),
#'  species = "speciesName"
#'  )
#'
#' @export
#' @import sf

trackSpp <- function(dat, inv, dorm , buff , buffGenet , clonal,
                     species = "Species",
                     site = "Site",
                     quad = "Quad",
                     year = "Year",
                     geometry = "geometry",
                     ...) {
  # argument checks ---------------------------------------------------------
  ## arguments

  #dat ## an sf d.f that contains all of the raw digitized map data
  # (in grasslandData format) for as many sites and quads as you'd like. Subset
  # by spp. and quad. before being passed to assign()

  #inv ## a list of the sampling years for each quadrat included in dat (in the
  # same format as grasslandInventory). Subset by quad before being passed to
  # assign()

  ## check the 'dat' and 'inv' arguments using the 'checkDat' function
  checkData <- checkDat(dat = dat, inv = inv,
                        species = species,
                        site = site,
                        quad = quad,
                        year = year,
                        geometry = geometry,
                        reformatDat = TRUE)

  dat <- checkData$dat
  inv <- checkData$inv
  usrNames <- checkData$userColNames

  # dorm ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates the
  # dormancy (in years) that is allowed. If multiple values, is subset by spp.
  # before being passed to assign()
  ## check dorm argument
  if(is.null(dorm)==TRUE) {
    stop("The 'dorm' argument must have a value.")
  } else {
    if (is.numeric(dorm)) { ## is the value of dorm a single numeric integer?
      if (dorm < 0 | ## dorm must be greater than or equal to 0
          round(dorm) != dorm | ## dorm must be a whole number
          length(dorm)!=1) { ## dorm must be a vector of length 1
stop("If 'dorm' is not specified for every species, it must be a single numeric
value that is a whole number greater than or equal to 0")
      }
    } else if (is.data.frame(dorm)) {
      if (sum(!names(dorm) %in% c("Species", "dorm")) == 0) {
        if(sum(!unique(dat$Species) %in% dorm$Species) > 0 | ## dorm must have
           # data for all species
           sum(is.na(dat$dorm)) > 0 | ## can't have NA values in dorm
           !is.numeric(dorm$dorm) | ## can't have non-numeric values for
           # dorm$dorm
           sum(dorm$dorm < 0) > 0 | ## can't be less than 0
           round(dorm$dorm) != dorm$dorm ## must be whole numbers
        ) {
stop("If the 'dorm' argument is specified by species, it must be a data.frame
that includes a 'Species' column with a row for every species in 'dat', and a
'dorm' column that contains positive, whole number values for each species with
no NAs.")
        }
      } else {
stop("If the 'dorm' argument is specifed by species, the column names must be
'Species' and 'dorm'")
      }
    } else {
stop("The 'dorm' argument must be either a single numeric value that is a whole
number greater than or equal to 0, OR a data.frame that has a 'Species' column
with values for each species in 'dat', and a 'dorm' column with numeric,
positive whole number values for each species.")
    }
  }

  #buff ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates the
  # buffer distance-- i.e. the distance a genet can move from year to year (in
  # the same units as distances in dat). If multiple values, is subset by spp.
  # before being passed to assign()
  ## check buff argument
  if(is.null(buff)==TRUE) {
    stop("The 'buff' argument must have a value.")
  } else {
    if (is.numeric(buff)) { ## is the value of buff a single numeric?
      if (buff < 0 | ## buff must be greater than or equal to 0
          buff > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buff must
          # not be larger than the dimensions of the quadrat
          length(buff)!=1) { ## buff must be a vector of length 1
stop("If 'buff' is not specified for every species, it must be a single numeric
value that is greater than or equal to 0")
      }
    } else if (is.data.frame(buff)) {
      if (sum(!names(buff) %in% c("Species", "buff")) == 0) {
        if(sum(!unique(dat$Species) %in% buff$Species) > 0 | ## buff must have
         # data for all species
         sum(is.na(dat$buff)) > 0 | ## can't have NA values in buff
         !is.numeric(buff$buff) | ## can't have non-numeric values for buff$buff
         sum(buff$buff < 0) > 0 | ## can't be less than 0
         round(buff$buff) != buff$buff ## must be whole numbers
      ) {
stop("If the 'buff' argument is specified by species, it must be a data.frame
that includes a 'Species' column with a row for every species in 'dat', and a
'buff' column that contains positive, numeric values for each species with no
NAs.")
        }
      } else {
stop("If the 'buff' argument is specifed by species, the column names must be
'Species' and 'buff'")
      }
    } else {
stop("The 'buff' argument must be either a single numeric value that is greater
than or equal to 0, OR a data.frame that has a 'Species' column with values for
each species in 'dat', and a 'buff' column with numeric values for each
species.")
    }
  }

  #buffGenet ## either a single value (applied to all spp.)or a data.frame with
  # the same number of rows as the number of species in dat that indicates how
  # close together ramets must be to be considered the same genet (in the same
  # units as distances in dat). If multiple values, is subset by spp. before
  # being passed to assign() (and then passed to groupByGenet())
  ## check buffGenet argument
  if(is.null(buffGenet)==TRUE) {
    stop("The 'buffGenet' argument must have a value.")
  } else {
    if (is.numeric(buffGenet)) { ## is the value of buffGenet a single numeric?
      if (buffGenet < 0 | ## buffGenet must be greater than or equal to 0
          buffGenet > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buffGenet
          # must not be larger than the dimensions of the quadrat
          length(buffGenet)!=1) { ## buffGenet must be a vector of length 1
stop("If 'buffGenet' is not specified for every species, it must be a single
numeric value that is greater than or equal to 0")
      }
    } else if (is.data.frame(buffGenet)) {
      if (sum(!names(buffGenet) %in% c("Species", "buffGenet")) == 0) {
      if(sum(!unique(dat$Species) %in% buffGenet$Species) > 0 | ## buffGenet
         # must have data for all species
         sum(is.na(dat$buffGenet)) > 0 | ## can't have NA values in buffGenet
         !is.numeric(buffGenet$buffGenet) | ## can't have non-numeric values for
         # buffGenet$buffGenet
         sum(buffGenet$buffGenet < 0) > 0 | ## can't be less than 0
         round(buffGenet$buffGenet) != buffGenet$buffGenet ## must be whole
         # numbers
      ) {
stop("If the 'buffGenet' argument is specified by species, it must be a
data.frame that includes a 'Species' column with a row for every species in
'dat', and a 'buffGenet' column that contains positive, numeric values for each
species with no NAs.")
      }
      } else {
stop("If the 'buffGenet' argument is specifed by species, the column names must
be 'Species' and 'buffGenet'")
      }
    } else {
stop("The 'buffGenet' argument must be either a single numeric value that is
greater than or equal to 0, OR a data.frame that has a 'Species' column with
values for each species in 'dat', and a 'buffGenet' column with numeric values
for each species.")
    }
  }

  #clonal ## either a single value (applied to all spp.) or a data.frame with
  # the same number of rows as the number of species in dat that indicates
  # whether or not a species is allowed to be clonal. One column contains the
  # species names, he second column contains clonal args. If multiple values, is
  # subset by spp. before being passed to assign()
  ## check clonal argument
  if(is.null(clonal)==TRUE) {
    stop("The 'clonal' argument must have a value.")
  } else {
    if (is.numeric(clonal)) { ## is the value of clonal a single numeric?
      if (clonal != 1 & clonal != 0 | ## clonal must be either 0 or 1
        !is.numeric(clonal) | ## clonal must be numeric
        length(clonal)!=1){ ## clonal must be a vector of length = 1
stop("If 'clonal' is not specified for every species, it must be a single
numeric value that is either 0 or 1.")
      }
    } else if (is.data.frame(clonal)) {
      if (sum(!names(clonal) %in% c("Species", "clonal")) == 0) {
      if(sum(!unique(dat$Species) %in% clonal$Species) > 0 | ## clonal
         # must have data for all species
         sum(is.na(clonal$clonal)) > 0 | ## can't have NA values in clonal
         !is.numeric(clonal$clonal) | ## can't have non-numeric values for
         # clonal$clonal
         sum((clonal$clonal != 1 & clonal$clonal != 0)) > 0 ## clonal values
         # must be either 0 or 1
      ) {
stop("If the 'clonal' argument is specified by species, it must be a data.frame
that includes a 'Species' column with a row for every species in 'dat', and a
'clonal' column that contains numeric values of either 0 or 1 for each species
with no NAs.")
      }
      } else {
stop("If the 'clonal' argument is specifed by species, the column names must be
'Species' and 'clonal'")
      }
    } else {
stop("The 'clonal' argument must be either a single numeric value that is
greater than or equal to 0, OR a data.frame that has a 'Species' column with
values for each species in 'dat', and a 'clonal' column that contains numeric
values of either 0 or 1 for each species with no NAs.")
    }
  }

  # work --------------------------------------------------------------------
  ## first, break the 'dat' d.f into two pieces, one with columns we need, and
  # another with columns that are 'extra'--will rejoin at the end of the
  # function
  ## assign an arbitary index number so we can re-join the data at the end
  dat$indexStore <- c(1:nrow(dat))
  ## put the 'extra' data into a separate d.f to 'store'
  datStore <- dat[, !names(dat) %in% c("Site", "Quad", "Year", "Species",
                                       "geometry")]
  datStore <- st_drop_geometry(datStore)
  ## put the data we actually need into the 'dat' object
  dat <- dat[,names(dat) %in% c("Site", "Quad", "Year", "Species",
                                "geometry", "indexStore")]

  ## get the 6-letter species code for each observation
  ## make a column in the d.f with the 6-letter species code for each row
  dat$sp_code_6  <- sapply(strsplit(dat$Species, " "), function(x)
    paste0(substr(toupper(x[1]), 1, 3), ## species name
           substr(toupper(x[2]), 1, 3)) ## genus name
  )

  ## get the basal area for each observation
  dat$basalArea <- st_area(dat)

  ## get the site(s)
  for(i in unique(dat$Site)) { ## i = the name of the site
    print(paste0("Site: ",i))
    ## get the quadrats w/in that site
    for (j in unique(dat[dat$Site==i,]$Quad)) { ## j = the name of the quad
      print(paste0("-- Quadrat: ",j))
      ## get the quadratInventory data for this quad
      if (is.list(inv)==TRUE) { ## if there is inv data for >1 quadrat
        invQuad <- inv[[j]]
      } else if (is.vector(inv)) { ## if there is inv data for only 1 quadrat
        invQuad <- inv
      }
      ## get the species w/in that quad
      for (k in unique(dat[dat$Site==i & dat$Quad==j,]$Species)) {
        ## k = the name of the species
        ## get the data for this site/quad/species
        datSp <- dat[dat$Site == i &
                       dat$Quad == j &
                       dat$Species == k,]

        ## get dorm value
        if(is.numeric(dorm)){
          dormK <- dorm
        } else if (is.data.frame(dorm)) {
          dormK <- dorm[dorm$Species==k,"dorm"]
        }

        ## get clonal value
        if(is.numeric(clonal)){
          clonalK <- clonal
        } else if (is.data.frame(clonal)) {
          clonalK <- clonal[clonal$Species==k,"clonal"]
        }

        ## get buff value
        if(is.numeric(buff)){
          buffK <- buff
        } else if (is.data.frame(buff)) {
          buffK <- buff[buff$Species==k,"buff"]
        }

        ## get buffGenet value
        if(is.numeric(buffGenet)){
          buffGenetK <- buffGenet
        } else if (is.data.frame(buffGenet)) {
          buffGenetK <- buffGenet[buffGenet$Species==k,"buffGenet"]
        }

        ## put this dataset into the 'assign' function
        datOut <- assign(dat = datSp,
                         inv = invQuad,
                         dorm = dormK,
                         buff = buffK,
                         buffGenet = buffGenetK,
                         clonal = clonalK
        )
        ## see if the output d.f exists yet (trackSpOut)
        ## if it does exist, then add datOut for the current spp. to the output
        if (exists("trackSppOut")==TRUE) {
          trackSppOut <- rbind(trackSppOut, datOut)
        }
        ## if not, then put datOut into trackSppOut b/c it is the first spp.
        if (exists("trackSppOut")==FALSE) {
          trackSppOut <- datOut
        }
        ## print the name of the species that was just finished
        if (k == unique(dat[dat$Site==i & dat$Quad==j,]$Species)[1]) {
          cat(paste0("---- Species: ",k))
        } else {
          cat(paste0("; ",k))
        }
      }
      ## notify user of last year of sampling (or last year of sampling before a
      # gap)
cat(paste0("Note: Individuals in year ", max(invQuad)," have a value of 'NA' in
the 'survives_tplus1' and 'size_tplus1' columns because ",max(invQuad), " is the
last year of sampling in this quadrat."))
      ## find years that exceed the 'dorm' gap
      invComp <- data.frame(inv = c(NA, invQuad), invNext = c(invQuad, NA))
      invComp$diff <- invComp$invNext - invComp$inv
      gapYears <- invComp[invComp$diff>dorm &
                            is.na(invComp$diff) == FALSE,"inv"]
      if (length(gapYears) > 0) {
cat(paste0("Note: Individuals in year(s) ", gapYears," have a value of 'NA' in
the 'survives_tplus1' and 'size_tplus1' columns because ", gapYears," is the
last year of sampling in this quadrat before a gap that exceeds the 'dorm'
argument."))
      }
    }
  }

  ## rejoin the trackSppOut d.f with the 'extra' data stored in 'datStore'
  trackSppOut <- merge(trackSppOut, datStore, by = "indexStore")
  ## remove the 'indexStore' value
  trackSppOut <- trackSppOut[,names(trackSppOut) != "indexStore"]

  ## re-name the appropriate columns in 'trackSppOut' data.frame with the
  # user-provided names of 'dat'
  ## from above, user-provided names are stored in 'usrNames'
  ## make a vector of default column names
  defaultNames <- c("Species", "Site", "Quad", "Year",  "geometry")

  ## reset the names for the columns that we changed to 'default' values
  names(trackSppOut)[which(names(trackSppOut) %in% defaultNames)] <-
  usrNames

  ## remove the '_USER' from the 'extra' column names
  names(trackSppOut) <- gsub(names(trackSppOut),
                             pattern = "_USER", replacement = "")

# output ------------------------------------------------------------------
return(trackSppOut)
}

# Testing -----------------------------------------------------------------
# dat <- grasslandData#[grasslandData$Site == "CO"
#                      #& grasslandData$Quad %in% c("unun_11","ungz_5a")
#                      #& grasslandData$Species == "Bouteloua gracilis",]
# names(dat)[1]<- "Species_Name"
# names(dat)[8] <- "location"
# inv <- grasslandInventory
# dorm <- 1
# buff <- 0.05
# buffGenet <- 0.005
# clonal <- data.frame(Species = unique(dat$Species),
#                      clonal = c(1,1,0,0,0,0,1,1,1,0,1,1,0,0))
#
# testOut <- trackSpp(dat, inv, dorm, buff, buffGenet, clonal,
#                     species = "Species_Name",
#                     quad = "location")

#
# testDat <- st_drop_geometry(dat)
# testDat$test <- "old"
# testOutputTest <- st_drop_geometry(testOut)
# testOutputTest$test <- "new"
#
# testTest <- full_join(testDat,testOutputTest, by = c("Species", "Clone",
#                             "Seedling", "Stems", "Basal", "Type", "Site",
#                             "Quad", "Year", "sp_code_4", "sp_code_6", "Area"))
# testBad <- testTest[is.na(testTest$test.y),]
#
# testBadSmall <- testTest[testTest$Site=="CO" & testTest$Quad == "unun_11" &
#                            testTest$Species == "Bouteloua gracilis",]

### AES make an example in the documentation that specifies all args as numeric,
# and another example where they specify all four arguments as data.frames

### be careful about warnings: because people freak out about them

###AES make fake errors in data to make sure that the argument checks are
# working correctly  (i.e. a random NA in a column, character for dorm, etc.)

###AES maybe also try to practice on larger subset of data.frame
