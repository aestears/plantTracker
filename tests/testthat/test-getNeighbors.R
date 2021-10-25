## get a test data.frame
dat <- grasslandData[grasslandData$Site == c("AZ") &
                       grasslandData$Species %in% c("Bouteloua rothrockii") &
                       grasslandData$Quad == "SG4",]
inv <- grasslandInventory[["SG4"]]
outDat <- assign(
  dat = dat,
  inv = inv,
  dorm = 1,
  buff = .05,
  buffGenet = 0.005,
  clonal = TRUE
)

locals <- getNeighbors(dat = outDat, buff = .1, method = "area",
                       compType = "oneSpp")

test_that("output is the same number of rows as the input", {
 expect_equal(object = nrow(locals), expected = nrow(dat))
})

test_that("output has only one more column than the input", {
  expect_equal(object = ncol(locals), expected = ncol(outDat)+1)
})
