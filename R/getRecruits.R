getRecruits <- function(dat,
                        species = "Species",
                        quad = "Quad",
                        site = "Site",
                        year = "Year",
                        trackID = "trackID",
                        recuit = "recruit") {

# work --------------------------------------------------------------------
  ## subset 'dat' to only get the individuals that are recruits
  datRecs <- dat[dat$recruit == 1 & is.na(dat$recruit) == FALSE,]

  ## remove the 'geometry' column, since it is not needed in this context
  datRecs <- sf::st_drop_geometry(datRecs)

  ## make sure that each genet has only one row in each year (if not, use
  # aggregateByGenet())
  if (nrow(unique(datRecs[,c("Year","trackID")])) != nrow(datRecs) ) { ## if the
    # number of year/trackID combos is not the same as the number of rows in
    # 'datRecs', then 'datRecs' needs to be aggregated by genet
    datRecs <- aggregateByGenet(dat = datRecs)
  }

  ## count the number of recruits in each year/species/quad/site combo
  datRecruits <- aggregate(x = datRecs[,c("recruit")], by = list(
    Site = datRecs$Site,
    Quad = datRecs$Quad,
    Species = datRecs$Species,
    Year = datRecs$Year
    ),
    FUN = length
  )

  names(datRecruits)[which(names(datRecruits) == "x")] <- "recruits"

# output ------------------------------------------------------------------
  return(datRecruits)
}
