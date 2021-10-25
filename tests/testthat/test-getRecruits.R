## get a test data.frame
dat <- grasslandData[grasslandData$Site == c("AZ") ,]
inv <- grasslandInventory
outDat <- trackSpp(
  dat = dat,
  inv = inv,
  dorm = 1,
  buff = .05,
  buffGenet = 0.005,
  clonal = TRUE
)

recruits <- getRecruits(dat = outDat, byGenet = TRUE)

#test
test_that("sum of recruits is the actual number of recruits in outDat", {
  expect_equal(object = sum(recruits$recruits),
               expected = sum(outDat$recruit, na.rm = TRUE))
})
