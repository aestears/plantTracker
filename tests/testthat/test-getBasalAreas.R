## get example data
dat <- grasslandData[grasslandData$Site == c("AZ") &
                       grasslandData$Species %in% c("Bouteloua rothrockii",
                                                    "Calliandra eriophylla"),]
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
test_that("percent basal area is absolute/quad basal area", {
  expect_equal(object = (areas$absolute_basalArea/areas$quad_basalArea)*100,
               expected = areas$percent_basalArea)
  }
)
