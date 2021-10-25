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

lambdas <- getLambda(dat = dat, inv = inv, method = "area",
                     species = "speciesName")

## tests
test_that("lambda calcs. are correct", {
  expect_equal(object = lambdas$absolute_basalArea_tplus1/lambdas$absolute_basalArea_t,
               expected = lambdas$lambda)
  }
)
