getLambda <- function(dat,
                      inv,
                      species = "Species",
                      quad = "Quad",
                      site = "Site",
                      year = "Year",
                      geometry = "geometry",
                      ...) {

# argument checking -------------------------------------------------------
  ## use the 'checkDat' function to make sure that 'dat' is in the
  # correct format
  checked <- checkDat(dat = dat, inv = inv, species = species, site = site,
                      quad = quad, year = year, geometry = geometry,
                      reformatDat = TRUE)
  ## save the checked data to use in the function
  dat <- checked$dat
  usrNames <- checked$userColNames
  ## save the inputed user names (will reassign later to the returned d.f.)
  inv <- checked$inv

# work --------------------------------------------------------------------
  ## get the basal areas for each species in each quadrat in each year
  # (interested in the "absolute_basalArea" column in the return of
  # getBasalAreas())


  datAreas <- getBasalAreas(dat = dat, inv = inv)

  datLambda <- data.frame(NULL)
  ## calculate lambda (area_t+1/area_t)
  ## loop through each site
  for (i in unique(datAreas$Site)) {
    ## loop through each quad
    for (j in unique(datAreas[datAreas$Site == i,"Quad"])) {
      ## loop through each species
      for (k in unique(datAreas[datAreas$Site == i & datAreas$Quad == j,
                                "Species"])) {
        ## get the data for each species
        datSpp <- datAreas[datAreas$Site == i &
                             datAreas$Quad == j &
                             datAreas$Species == k,
                           c("Site", "Quad", "Species",
                             "Year", "absolute_basalArea")]
        ## make sure that years are in order
        datSpp <- datSpp[order(datSpp$Year),]
        ## print a message if there is only one year of data (need at least two
        # time points to calculate lambda)
        if (nrow(datSpp) <= 1) {
          print(paste0("A value of lambda for ", k, " in quadrat ", j,
          " cannot be calculated, since data was only collected for one",
          " year. Data from at least two time points are required to",
          " calculate lambda."))
        }

        ## calculate lambda (basalArea_tplus1 / basalArea_t)
        names(datSpp)[which(names(datSpp) == "Year")] <- "Year_t"
        names(datSpp)[which(names(datSpp) == "absolute_basalArea")] <-
          "absolute_basalArea_t"
        ## get the 'Year_tplus1' value (the next sequential year of sampling,
        # NOT necessarily the next sequential year)
        datSpp$Year_tplus1 <- c(datSpp$Year_t[2:length(datSpp$Year_t)], NA)
        ## get the 'Yabsolute_basalArea_tplus1' value (the next sequential year
        # of sampling, NOT necessarily the next sequential year)
        datSpp$absolute_basalArea_tplus1 <-
          c(datSpp$absolute_basalArea_t[2:length(datSpp$absolute_basalArea_t)],
            NA)
        ## finally calculate lambda itself
        datSpp$lambda <- datSpp$absolute_basalArea_tplus1/
          datSpp$absolute_basalArea_t

        ## save the output
        if (nrow(datLambda) == 0) {
          datLambda <- datSpp
        } else {
          datLambda <- rbind(datLambda, datSpp)
        }
      }
    }
  }
## rename the columns to the user-defined columns
## from above, user-provided names are stored in 'usrNames'
defaultNames <- c("Species", "Site", "Quad")
## reset the names for the columns that we changed to 'default' values
names(datLambda)[match(defaultNames, names(datLambda))] <-
    usrNames[1:3]

# return ------------------------------------------------------------------
return(datLambda)
}

# testing -----------------------------------------------------------------
#
# dat <- grasslandData[grasslandData$Site == "CO" & grasslandData$Year %in% c(1998:2002),]
# names(dat)[1] <- "speciesName"
# inv <- grasslandInventory[unique(dat$Quad)]
# outDat <- trackSpp(dat = dat, inv = inv, dorm = 1, buff = .05,buffGenet = 0.005,clonal = data.frame("Species" = unique(dat$speciesName),"clonal" = c(1,0)), species = "speciesName",aggregateByGenet = TRUE)
# getLambda(dat = outDat, inv = inv, species = "speciesName")
