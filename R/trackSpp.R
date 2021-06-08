#' Tracks genets through time for multiple species and sites
#'
#' @return
#' @export
#'
#' @examples

trackSpp <- function(dat, inv, dorm = NULL, buff = NULL, buffGenet = NULL,
                     clonal = NULL, sppArgs = NULL, ...) {

  # argument checks ---------------------------------------------------------
  ## arguments

  #dat ## an sf d.f that contains all of the raw digitized map data
  # (in grasslandData format) for as many sites and quads as you'd like. Subset
  # by spp. and quad. before being passed to assign()

  #inv ## a list of the sampling years for each quadrat included in dat (in the
  # same format as grasslandInventory). SUbset by quad before being passed to
  # assign()

  #sppArgs ## a data.frame that contains species-specific arguments for dorm,
  # clonal, buff, and buffGenet. Does not need to contain values for every
  # single species. If a species is not present in this d.f, then default values
  # of args. are used. This argument is optional--if not present, then default
  # arg. values are used. If there is a row for a species, than all columns must
  # have a value for that species.
  ## checking the sppArgs argument
  ## if the the sppArgs has all value for all species present in dat, then
  # defaults for dorm, clonal, buff, and buffGenet are not required!
  ## is there a sppArgs d.f?
  if (is.null(sppArgs)==FALSE) {
    ## check to make sure that the spp. listed in sppArgs are actually present
    # in dat
    if (is.null(sppArgs$Species)==TRUE) {
      warning("the sppArgs data.frame must contain a 'Species' column with at
              least one species name. If it s not populated, then default
              argument values will be used.")
    }
    if (sum(!(sppArgs$Species %in% unique(dat$Species))) > 0) {
      stop("One or more species names in the 'sppArgs' data.frame is not present
      in the 'dat' data.frame. Check for spelling errors!")
    }
    ## see if sppArgs has rows for all species in dat, and values for all
    # arguments!
    if (sum(!(unique(dat$Species) %in% sppArgs$Species)) == 0 & ## are values
        # for all species in dat present in sppArgs?
        sum(!(c("buff", "buffGenet", "clonal", "dorm") %in% names(sppArgs)))
        == 0 ## are
        # there columns with values for all possible args?
        ) {
      ## if there ARE sppArgs for each species in dat and values for each
      # argument, then don't need other default args (and don't need to check
      # them)--actually, set them to NA!
      dorm <- NA
      buff <- NA
      buffGenet <- NA
      clonal <- NA

      ## but DO need to check the values of each of the arg. values in sppArgs
      ## check clonal args
      if(sum(sppArgs$clonal != 1 & sppArgs$clonal != 0 | ## clonal must be 0 or 1
         !is.numeric(sppArgs$clonal)) > 0) ## clonal must be numeric
        {
        stop("'clonal' argument must be a numeric column in the
             sppArgs data.frame containing boolean values")
      }
      ## check dorm args
      if(sum(sppArgs$dorm < 0 | ## dorm must be greater than or equal to 0
         !is.numeric(sppArgs$dorm) | ## dorm must be numeric
         round(sppArgs$dorm) != dorm) > 0) ## dorm must be a whole number)
         {
        stop("'dorm' argument must be a a numeric column in the sppArgs
             data.frame, containing only positive, whole numbers")
      }
      ## check buff args
      if(sum(!is.numeric(sppArgs$buff) | ## buff must be numeric
         sppArgs$buff > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buff must not
         # be larger than the dimensions of the quadrat
         sppArgs$buff < 0 ## buff must be greater than or equal to zero
      ) > 0) {
        stop("'buff' argument must be anumeric column in the sppArgs data.frame,
        and cannot exceed the maximum dimensions of the quadrat")
      }
      ## check buffGenet args
      if(sum(!is.numeric(sppArgs$buffGenet) | ## buffGenet must be numeric
         sppArgs$buffGenet > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buffGenet must not
         # be larger than the dimensions of the quadrat
         sppArgs$buffGenet < 0 ## buffGenet must be greater than or equal to zero
      ) > 0) {
        stop("'buffGenet' argument must be a numeric column in the sppArgs
             data.frame, and cannot exceed the maximum dimensions of the
             quadrat")
      }

    } else { ## if there are NOT sppArgs for each species in dat, then you DO
      # need other default args, and need to check them! (do that later)
      ## still need to check the values for sppArgs that do exist!
      ## check clonal args
      if(is.null(sppArgs$clonal) == FALSE) { ## does the 'clonal' column exist?
        # if not, don't need to check it!
        if(sum(sppArgs$clonal != 1 & sppArgs$clonal != 0 | ## clonal  = 0 or 1
               !is.numeric(sppArgs$clonal)) > 0) ## clonal must be numeric
        {
          stop("'clonal' argument must be a numeric column in the
             sppArgs data.frame containing boolean values")
        }
      }

      ## check dorm args
      if(is.null(sppArgs$dorm) == FALSE) { ## does the 'dorm' column exist?--
        # if not, don't need to check it!
        if(sum(sppArgs$dorm < 0 | ## dorm must be greater than or equal to 0
               !is.numeric(sppArgs$dorm) | ## dorm must be numeric
               round(sppArgs$dorm) != dorm) > 0) ## dorm must be a whole number)
        {
          stop("'dorm' argument must be a a numeric column in the sppArgs
             data.frame, containing only positive, whole numbers")
        }
      }
      ## check buff args
      if(is.null(sppArgs$buff) == FALSE){ ## does the 'buff' column exist?--
        # if not, don't need to check it!
        if(sum(!is.numeric(sppArgs$buff) | ## buff must be numeric
               sppArgs$buff > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buff
               # must not be larger than the dimensions of the quadrat
               sppArgs$buff < 0 ## buff must be greater than or equal to zero
        ) > 0) {
          stop("'buff' argument must be anumeric column in the sppArgs
          data.frame, and cannot exceed the maximum dimensions of the quadrat")
        }
      }
      ## check buffGenet args
      if(is.null(sppArgs$buffGenet) == FALSE) { ## does the 'buffGenet' column
        # exist?-- if not, don't need to check it!
        if(sum(!is.numeric(sppArgs$buffGenet) | ## buffGenet must be numeric
               sppArgs$buffGenet > max(st_bbox(dat)[c("xmax", "ymax")]) |
               ## buffGenet must notbe larger than the dimensions of the quadrat
               sppArgs$buffGenet < 0 ## buffGenet must be >= to zero
        ) > 0) {
          stop("'buffGenet' argument must be a numeric column in the sppArgs
             data.frame, and cannot exceed the maximum dimensions of the
             quadrat")
        }
      }
    }
  }  ## if there is NOT a sppArgs d.f, then defaults for other args are
    # required (and need to be checked) (do that later, outside this 'else')

  # dorm ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates the
  # dormancy (in years) that is allowed. If multiple values, is subset by spp.
  # before being passed to assign()
  ## check dorm argument
  if(is.null(dorm)==TRUE) {
    stop("The 'dorm' argument must have a value, since dorm values for all
         species are not specified in the 'sppArg' data.frame.")
  } else {
    if(is.na(dorm)==FALSE){
      if(dorm < 0 | ## dorm must be greater than or equal to 0
         !is.numeric(dorm) | ## dorm must be numeric
         round(dorm) != dorm | ## dorm must be a whole number
         length(dorm)!=1){ ## dorm must be a vector of length = 1
        stop("'dorm' argument must be a a numeric vector of length = 1, containing
      a positive, whole number")
      }
    }
  }

  #buff ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates the
  # buffer distance-- i.e. the distance a genet can move from year to year (in
  # the same units as distances in dat). If multiple values, is subset by spp.
  # before being passed to assign()
  ## check buff argument
  if(is.null(buff)==TRUE) {
    stop("The 'buff' argument must have a value, since buff values for all
         species are not specified in the 'sppArg' data.frame.")
  } else {
    if (is.na(buff)==FALSE) {
      if(!is.numeric(buff) | ## buff must be numeric
         buff > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buff must not be larger
         # than the dimensions of the quadrat
         buff < 0 ## buff must be greater than or equal to zero
      ) {
        stop("'buff' argument must be numeric and cannot exceed the maximum
         dimensions of the quadrat")
      }
    }
  }

  #buffGenet ## either a single value (applied to all spp.)or a data.frame with
  # the same number of rows as the number of species in dat that indicates how
  # close together ramets must be to be considered the same genet (in the same
  # units as distances in dat). If multiple values, is subset by spp. before
  # being passed to assign() (and then passed to groupByGenet())
  ## check buffGenet argument
  if(is.null(buffGenet)==TRUE) {
    stop("The 'buffGenet' argument must have a value, since buffGenet values for
    all species are not specified in the 'sppArg' data.frame.")
  } else {
    if (is.na(buffGenet) == FALSE) {
      if(!is.numeric(buffGenet) | ## buffGenet must be numeric
         buffGenet > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buffGenet must
         # not be larger than the dimensions of the quadrat
         buffGenet < 0 ## buffGenet must be greater than or equal to zero
      ) {
        stop("'buffGenet' argument must be numeric and cannot exceed the maximum
         dimensions of the quadrat")
      }
    }
  }

  #clonal ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates whether
  # or not a species is allowed to be clonal. One column contains the species
  # names, he second column contains clonal args. If multiple values, is subset
  # by spp. before being passed to assign()
  ## check clonal argument
  if(is.null(clonal)==TRUE) {
    stop("The 'clonal' argument must have a value, since clonal values for
    all species are not specified in the 'sppArg' data.frame.")
  } else {
    if (is.na(clonal)==FALSE) {
      if(clonal != 1 & clonal != 0 | ## clonal must be either 0 or 1
         !is.numeric(clonal) | ## clonal must be numeric
         length(clonal)!=1){ ## clonal must be a vector of length = 1
        stop("'clonal' argument must be a numeric boolean vector of length = 1")
      }
    }
  }

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
      }
      else if (is.vector(inv)) { ## if there is inv data for only 1 quadrat
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

        ## are there species-specific arg values at all?
        if (exists('sppArgs')==TRUE) {
          ## are there species-specific values for this species?
          if (sum(sppArgs$Species %in% k) > 0) {
            ## if there ARE sp-specific values for this species, check that there
            # are values for each argument

            ## is there a sp-specific value for dorm?
            if (is.null(sppArgs[sppArgs$Species==k,"dorm"])==FALSE) {
              ## there IS a sp-specific value for dorm
              dormK <- sppArgs[sppArgs$Species==k,"dorm"]
            } else {
              ## there is NOT a sp-specific value for dorm
              dormK <- dorm
            }

            ## is there a sp-specific value for buff?
            if (is.null(sppArgs[sppArgs$Species==k,"buff"])==FALSE) {
              ## there IS a sp-specific value for buff
              buffK <- sppArgs[sppArgs$Species==k,"buff"]
            } else {
              ## there is NOT a sp-specific value for buff
              buffK <- buff
            }

            ## is there a sp-specific value for buffGenet?
            if (is.null(sppArgs[sppArgs$Species==k,"buffGenet"])==FALSE) {
              ## there IS a sp-specific value for buffGenet
              buffGenetK <- sppArgs[sppArgs$Species==k,"buffGenet"]
            } else {
              ## there is NOT a sp-specific value for buffGenet
              buffGenetK <- buffGenet
            }

            ## is there a sp-specific value for clonal?
            if (is.null(sppArgs[sppArgs$Species==k,"clonal"])==FALSE) {
              ## there IS a sp-specific value for clonal
              clonalK <- sppArgs[sppArgs$Species==k,"clonal"]
            } else {
              ## there is NOT a sp-specific value for clonal
              clonalK <- clonal
            }

          } else {
            ## if there are NOT sp-specific values for this species, then use
            # default values
            dormk <- dorm
            buffk <- buff
            buffGenetk <- buffGenet
            clonalK <- clonal
          }
        } else {
          ## if there are not species-specific values at all
          dormk <- dorm
          buffk <- buff
          buffGenetk <- buffGenet
          clonalK <- clonal
        }


        ## put this dataset into the 'assign' function
        datOut <- assign(dat = datSp,
                         inv = invQuad,
                         dorm = dormK,
                         buff = buffK,
                         buffGenet = buffGenetK,
                         clonal = clonalK
        )

        ###AES### add an output showing the progress
        ## see if the output d.f exists yet (trackSpOut)
        ## if it does exist, then add datOut for the current spp. to the output
        if (exists("trackSppOut")==TRUE) {
          trackSppOut <- rbind(trackSppOut, datOut)
        }
        ## if not, then put datOut into trackSppOut b/c it is the first spp.
        if (exists("trackSppOut")==FALSE) {
          trackSppOut <- datOut
        }
        #print(paste0("Finished with ",k,", quadrat ",j,", ", i," site"))
      }
      #print(paste0("---Finished with quadrat ",j,", ",i," site"))
    }
    #print(paste0("-----Finished with ",i,' site'))
  }
# output ------------------------------------------------------------------
return(trackSppOut)
}

# Testing -----------------------------------------------------------------

dat <- grasslandData#[grasslandData$Site == "CO"
                     #& grasslandData$Quad %in% c("unun_11","ungz_5a")
                     #& grasslandData$Species == "Bouteloua gracilis",]
inv <- grasslandInventory
dorm <- 1
buff <- 0.05
buffGenet <- 0.005
clonal <- 1
sppArgs <- data.frame('Species' = unique(dat$Species)[1:5],
                      'clonal' = c(1,1,0,0,0),
                      #'buffGenet' = c(.005, .005, 0, 0, 0),
                      'dorm' = c(1,1,1,1,1),
                      'buff' = c(.05,.05,.05,.05,.05))

testOut <- trackSpp(dat, inv, dorm, buff, buffGenet, clonal, sppArgs)


testDat <- st_drop_geometry(dat)
testDat$test <- "old"
testOutputTest <- st_drop_geometry(testOut)
testOutputTest$test <- "new"

testTest <- full_join(testDat,testOutputTest, by = c("Species", "Clone", "Seedling", "Stems", "Basal", "Type", "Site", "Quad", "Year", "sp_code_4", "sp_code_6", "Area"))
testBad <- testTest[is.na(testTest$test.y),]

testBadSmall <- testTest[testTest$Site=="CO" & testTest$Quad == "unun_11" & testTest$Species == "Bouteloua gracilis",]
