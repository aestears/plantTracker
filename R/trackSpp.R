#' Tracks genets through time for multiple species and sites
#'
#' @return
#' @export
#'
#' @examples

trackSpp <- function(dat, inv, dorm , buff , buffGenet , clonal,
                     species = "Species",
                     site = "Site",
                     quad = "Quad",
                     year = "Year",
                     geometry = "geometry",
                     ...) {

  ###AES start working on documentation for this function
  # -- can send assign and trackSpp both to the same help page! make help more
  # generalized, then add details about the differences to the 'details' section

  # argument checks ---------------------------------------------------------
  ## arguments

  #dat ## an sf d.f that contains all of the raw digitized map data
  # (in grasslandData format) for as many sites and quads as you'd like. Subset
  # by spp. and quad. before being passed to assign()

  ## check the 'dat' and 'inv' arguments using the 'checkDat' function
  dat <- checkDat(dat = dat, inv = inv, datNames = datNames,
                  trackerFormat = TRUE, inheritFromTrackSpp = FALSE,
                  printGoAhead = FALSE)$dat

  #inv ## a list of the sampling years for each quadrat included in dat (in the
  # same format as grasslandInventory). Subset by quad before being passed to
  # assign()

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
        stop("If 'dorm' is not specified for every species, it must be a single
             numeric value that is a whole number greater than or equal to 0")
      }
    } else if (is.data.frame(dorm)) {
      if (sum(!names(dorm) %in% c("Species", "dorm")) == 0) {
        if(sum(!unique(dat$Species) %in% dorm$Species) > 0 | ## dorm must have
           # data for all species
           sum(is.na(dat$dorm)) > 0 | ## can't have NA values in dorm
           !is.numeric(dorm$dorm) | ## can't have non-numeric values for dorm$dorm
           sum(dorm$dorm < 0) > 0 | ## can't be less than 0
           round(dorm$dorm) != dorm$dorm ## must be whole numbers
        ) {
          stop("If the 'dorm' argument is specified by species, it must be a
             data.frame that includes a 'Species' column with a row for every
             species in 'dat', and a 'dorm' column that contains positive, whole
             number values for each species with no NAs.")
        }
      } else {
        stop("If the 'dorm' argument is specifed by species, the column names
             must be 'Species' and 'dorm'")
      }
    } else {
      stop("The 'dorm' argument must be either a single numeric value that is a
           whole number greater than or equal to 0, OR a data.frame that has a
           'Species' column with values for each species in 'dat', and a 'dorm'
           column with numeric, positive whole number values for each species.")
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
        stop("If 'buff' is not specified for every species, it must be a single
             numeric value that is greater than or equal to 0")
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
        stop("If the 'buff' argument is specified by species, it must be a
             data.frame that includes a 'Species' column with a row for every
             species in 'dat', and a 'buff' column that contains positive,
             numeric values for each species with no NAs.")
        }
      } else {
        stop("If the 'buff' argument is specifed by species, the column names
             must be 'Species' and 'buff'")
      }
    } else {
      stop("The 'buff' argument must be either a single numeric value that is
      greater than or equal to 0, OR a data.frame that has a
           'Species' column with values for each species in 'dat', and a 'buff'
           column with numeric values for each species.")
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
        stop("If 'buffGenet' is not specified for every species, it must be a
        single numeric value that is greater than or equal to 0")
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
             data.frame that includes a 'Species' column with a row for every
             species in 'dat', and a 'buffGenet' column that contains positive,
             numeric values for each species with no NAs.")
      }
      } else {
        stop("If the 'buffGenet' argument is specifed by species, the column
        names must be 'Species' and 'buffGenet'")
      }
    } else {
      stop("The 'buffGenet' argument must be either a single numeric value that
      is greater than or equal to 0, OR a data.frame that has a 'Species' column
      with values for each species in 'dat', and a 'buffGenet' column with
           numeric values for each species.")
    }
  }

  #clonal ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates whether
  # or not a species is allowed to be clonal. One column contains the species
  # names, he second column contains clonal args. If multiple values, is subset
  # by spp. before being passed to assign()
  ## check clonal argument
  if(is.null(clonal)==TRUE) {
    stop("The 'clonal' argument must have a value.")
  } else {
    if (is.numeric(clonal)) { ## is the value of clonal a single numeric?
      if (clonal != 1 & clonal != 0 | ## clonal must be either 0 or 1
        !is.numeric(clonal) | ## clonal must be numeric
        length(clonal)!=1){ ## clonal must be a vector of length = 1
        stop("If 'clonal' is not specified for every species, it must be a
        single numeric value that is either 0 or 1.")
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
        stop("If the 'clonal' argument is specified by species, it must be a
             data.frame that includes a 'Species' column with a row for every
             species in 'dat', and a 'clonal' column that contains numeric
             values of either 0 or 1 for each species with no NAs.")
      }
      } else {
        stop("If the 'clonal' argument is specifed by species, the column
        names must be 'Species' and 'clonal'")
      }
    } else {
      stop("The 'clonal' argument must be either a single numeric value that
      is greater than or equal to 0, OR a data.frame that has a 'Species' column
      with values for each species in 'dat', and a 'clonal' column that contains numeric
             values of either 0 or 1 for each species with no NAs.")
    }
  }

  #Species
  #Site
  #Quad
  #Year
  #sp_code_6
  #geometry

  # work --------------------------------------------------------------------
  ## get the site(s)
  for(i in unique(dat$Site)) { ## i = the name of the site
    print(paste0("Site: ",i))
    ## get the quadrats w/in that site
    for (j in unique(dat[dat$Site==i,]$Quad)) { ## j = the name of the quad
      print(paste0("--- Quadrat: ",j))
      ## get the quadratInventory data for this quad
      if (is.list(inv)==TRUE) { ## if there is inv data for >1 quadrat
        invQuad <- inv[[j]]
      } else if (is.vector(inv)) { ## if there is inv data for only 1 quadrat
        invQuad <- inv
      }
      ## get the species w/in that quad
      for (k in unique(dat[dat$Site==i & dat$Quad==j,]$Species)) {
        ## k = the name of the species
        print(paste0("------- Species: ",k))

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
      }
      ## notify user of last year of sampling (or last year of sampling before a
      # gap)
      print(paste0("Note: Individuals in year ", max(invQuad)," have a value of 'NA' in the 'survives_tplus1' and 'size_tplus1' columns because ", max(invQuad),
                   " is the last year of sampling in this quadrat."))
      ## find years that exceed the 'dorm' gap
      invComp <- data.frame(inv = c(NA, invQuad), invNext = c(invQuad, NA))
      invComp$diff <- invComp$invNext - invComp$inv
      gapYears <- invComp[invComp$diff>dorm & is.na(invComp$diff) == FALSE,"inv"]
      if (length(gapYears) > 0) {
        print(paste0("Note: Individuals in year(s) ", gapYears," have a value of 'NA' in the 'survives_tplus1' and 'size_tplus1' columns because ", gapYears,
                     " is the last year of sampling in this quadrat before a gap that exceeds the 'dorm' argument."))
      }
    }
  }

  ## re-name the appropriate columns in 'trackSppOut' data.frame with the
  # user-provided names of 'dat'
userDatNames <- checkDat(dat, inv, trackerFormat = FALSE,
                         inheritFromTrackSpp = TRUE,
           printGoAhead = FALSE)$userDatNames

defaultDatNames <- c("Species", "Site", "Quad", "Year", "sp_code_6", "geometry")

names(trackSppOut)[which(names(trackSppOut) %in% defaultDatNames)] <-
  userDatNames


# output ------------------------------------------------------------------
return(trackSppOut)
}

# Testing -----------------------------------------------------------------
# dat <- grasslandData#[grasslandData$Site == "CO"
#                      #& grasslandData$Quad %in% c("unun_11","ungz_5a")
#                      #& grasslandData$Species == "Bouteloua gracilis",]
# inv <- grasslandInventory
# dorm <- 1
# buff <- 0.05
# buffGenet <- 0.005
# clonal <- data.frame(Species = unique(dat$Species),
#                      clonal = c(1,1,0,0,0,0,1,1,1,0,1,1,0,0))
#
# testOut <- trackSpp(dat, inv, dorm, buff, buffGenet, clonal)
#
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
