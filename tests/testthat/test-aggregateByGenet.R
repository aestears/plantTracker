## get a test data.frame
dat <- grasslandData[grasslandData$Site == c("AZ") &
                       grasslandData$Species %in% c("Bouteloua rothrockii", "Calliandra eriophylla"),]
names(dat)[1] <- "speciesName"
inv <- grasslandInventory[unique(dat$Quad)]
outDat <- trackSpp(dat = dat,
                   inv = inv,
                   dorm = 1,
                   buff = .05,
                   buffGenet = 0.005,
                   clonal = data.frame("Species" = unique(dat$speciesName),
                                       "clonal" = c(TRUE,FALSE)),
                   species = "speciesName",
                   aggByGenet = FALSE
)
## get a test aggregateByGenet() output
testDat <- aggregateByGenet(dat = outDat,
                            species = 'speciesName')
## test that the resulting d.f has no non-unique genets
test_that("result has no non-unique genets", {
  expect_equal(nrow(testDat), nrow(
    unique(testDat[,c("Site", "Quad", "Year", "speciesName", "trackID")]
    )
  )
  )
}
)
