## get example data
dat <- grasslandData[grasslandData$Site == "CO" &
                       grasslandData$Year %in% c(1998:2002), ]
names(dat)[1] <- "speciesName"
inv <- grasslandInventory[unique(dat$Quad)]
outDat <- trackSpp(
  dat = dat,
  inv = inv,
  dorm = 1,
  buff = .05,
  buffGenet = 0.005,
  clonal = data.frame(
    "Species" = unique(dat$speciesName),
    "clonal" = c(TRUE, FALSE)
  ),
  species = "speciesName",
  aggByGenet = TRUE
)
areas <- getBasalAreas(dat = outDat, inv = inv,
              species = "speciesName")

## tests
