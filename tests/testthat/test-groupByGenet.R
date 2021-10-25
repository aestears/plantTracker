## get a test data.frame
dat <- grasslandData[grasslandData$Site == c("AZ") &
                       grasslandData$Species %in% c("Bouteloua rothrockii") &
                       grasslandData$Quad == "SG4",]
inv <- grasslandInventory
outDat <- trackSpp(
  dat = dat,
  inv = inv,
  dorm = 1,
  buff = .05,
  buffGenet = 0.005,
  clonal = TRUE, aggByGenet = FALSE
)

genets <- aggregateByGenet(dat = outDat)

#test
test_that("test that there are only as many rows as unique trackID values", {
  expect_equal(object = nrow(genets),
               expected = nrow(unique(genets[,c("trackID", "Year")])))
})
