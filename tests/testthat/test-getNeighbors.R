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

# using the 'area' method
locals_1 <- getNeighbors(dat = outDat, buff = .1, method = "area",
                       compType = "oneSpp")

# using the 'count' method
locals_2 <- getNeighbors(dat = outDat, buff = .1, method = "count",
                         compType = "oneSpp")

# tests:
test_that("output is the same number of rows as the input", {
 expect_equal(object = nrow(locals_1), expected = nrow(dat))
})
test_that("output is the same number of rows as the input", {
  expect_equal(object = nrow(locals_2), expected = nrow(dat))
})

test_that("output has only one more column than the input (for method = 'count')", {
  expect_equal(object = ncol(locals_2), expected = ncol(outDat)+1)
})
test_that("output has only two more columns than the input (for method = 'area')", {
  expect_equal(object = ncol(locals_1), expected = ncol(outDat)+2)
})
