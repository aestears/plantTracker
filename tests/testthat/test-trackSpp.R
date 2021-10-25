## get a test data.frame
dat <- grasslandData[grasslandData$Site == c("AZ") &
                       grasslandData$Species %in% c("Bouteloua rothrockii") ,]
inv <- grasslandInventory
outDat <- trackSpp(
  dat = dat,
  inv = inv,
  dorm = 1,
  buff = .05,
  buffGenet = 0.005,
  clonal = TRUE
)

## test that the area of each individual is not larger than the area of the
# bounding box
test_that("individual area is not larger than the quadrat area", {
  apply(
    X = outDat,
    MARGIN = 1,
    FUN = function(x)
      expect_lte(
        object = x$basalArea_genet,
        expected =
          sf::st_bbox(outDat)["xmax"] * sf::st_bbox(outDat)["ymax"]
      )
  )
})

## test that the output data has the required column names
test_that("output data has required column names", {
  expect_equal(object = sum(
    c(
      "Species",
      "Site",
      "Quad",
      "Year",
      "trackID",
      "age",
      "size_tplus1",
      "recruit",
      "survives_tplus1",
      "basalArea_genet",
      "geometry",
      "nearEdge"
    ) %in% names(outDat)
  ),
  expected = 12)
})
