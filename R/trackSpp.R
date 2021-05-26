#' Tracks genets through time for multiple species and sites
#'
#' @return
#' @export
#'
#' @examples

trackSpp <- function(dat, inv, dorm, buff, buffGenet, clonal, ...) {

  # argument checks ---------------------------------------------------------
  ## arguments

  #dat ## an sf d.f that contains all of the raw digitized map data
  # (in grasslandData format) for as many sites and quads as you'd like. Subset
  # by spp. and quad. before being passed to assign()

  #inv ## a list of the sampling years for each quadrat included in dat (in the
  # same format as grasslandInventory). SUbset by quad before being passed to
  # assign()

  #dorm ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates the
  # dormancy (in years) that is allowed. If multiple values, is subset by spp.
  # before being passed to assign()

  #buff ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates the
  # buffer distance-- i.e. the distance a genet can move from year to year (in
  # the same units as distances in dat). If multiple values, is subset by spp.
  # before being passed to assign()

  #buffGenet ## either a single value (applied to all spp.)or a data.frame with
  # the same number of rows as the number of species in dat that indicates how
  # close together ramets must be to be considered the same genet (in the same
  # units as distances in dat). If multiple values, is subset by spp. before
  # being passed to assign() (and then passed to groupByGenet())

  #clonal ## either a single value (applied to all spp.) or a data.frame with the
  # same number of rows as the number of species in dat that indicates whether
  # or not a species is allowed to be clonal. One column contains the species
  # names, he second column contains clonal args. If multiple values, is subset
  # by spp. before being passed to assign()

  # work --------------------------------------------------------------------
  ## get the site(s)
  for(i in unique(dat$Site)) { ## i = the name of the site

    ## get the quadrats w/in that site
    for (j in unique(dat[dat$Site==i,]$Quad)) { ## j = the name of the quad
      ## get the quadratInventory data for this quad
      invQuad <- inv[[j]]

      ## get the species w/in that quad
      for(k in unique(dat[dat$Site==i & dat$Quad==j,]$Species)) {
        ## k = the name of the species
        datSp <- dat[dat$Site == i &
                       dat$Quad == j &
                       dat$Species == k,]

        ## put this dataset into the 'assign' function
        datOut <- assign(dat = datSp,
                         inv = invQuad,
                         dorm = ifelse(is.vector(dorm),
                                       yes = dorm,
                                       no = dorm[dorm$Species==k, "dorm"]),
                         buff = ifelse(is.vector(buff),
                                       yes = buff,
                                       no = buff[buff$Species==k,"buff"]),
                         buffGenet = ifelse(is.vector(buffGenet),
                                            yes = buffGenet,
                                            no = buffGenet[buffGenet$Species==k,
                                                           "buffGenet"]),
                         clonal = ifelse(is.vector(clonal),
                                         yes = clonal,
                                         no = clonal[clonal$Species==k,
                                                     "clonal"])
        )
        ## see if the output d.f exists yet (trackSpOut)
        ## if not, then put datOut into trackSppOut b/c it is the first spp.
        if(exists("trackSppOut")==FALSE) {
          trackSppOut <- datOut
        }
        ## if it does exist, then add datOut for the current spp. to the output
        if(exists("trackSppOut")==TRUE){
          trackSppOut <- rbind(trackSppOut, datOut)
        }
      }
    }
  }

# output ------------------------------------------------------------------
return(trackSppOut)
}



# Testing -----------------------------------------------------------------
#
dat <- grasslandData
inv <- grasslandInventory
dorm <- 1
buff <- 0.05
buffGenet <- 0.005
clonal <- data.frame('Species' = unique(dat$Species),
                     "clonal" = c(1,1,0,0,0,0,1,1,1,0,1,1,0,0))

testOut <- trackSpp(dat, inv, dorm, buff, buffGenet, clonal)


testDat <- st_drop_geometry(dat)
testDat$test <- "old"
testOutputTest <- st_drop_geometry(testOut)
testOutputTest$test <- "new"

testTest <- full_join(testDat,testOutputTest, by = c("Species", "Clone", "Seedling", "Stems", "Basal", "Type", "Site", "Quad", "Year", "sp_code_4", "sp_code_6", "Area"))
testBad <- testTest[is.na(testTest$test.y),]

testBadSmall <- testTest[testTest$Site=="CO" & testTest$Quad == "unun_11" & testTest$Species == "Bouteloua gracilis",]

 ###AES### for some reason is missing quite a few obs? need to figure out why###

