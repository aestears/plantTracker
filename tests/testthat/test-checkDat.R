## get a test data.frame
dat <- grasslandData[grasslandData$Site == c("AZ") &
                       grasslandData$Species %in% c("Bouteloua rothrockii",
                                                    "Calliandra eriophylla"),]
inv <- grasslandInventory[unique(dat$Quad)]
names(dat)[c(1, 9)] <- c("speciesName", "date")

## get output of checkDat
outDat <- checkDat(
  dat = dat,
  inv = inv,
  species = "speciesName",
  year = "date",
  reformatDat = TRUE
)

## check that the function output is a list
test_that("output data is a list", {
  expect_true(object = is.list(outDat))
})

## check that the data.frame part of the list is an sf data.frame with the
# required data.types (POLYGON or MULTIPOLYGON)
test_that("output data.frame is an sf data.frame with required data types", {
  expect_equal(object = sum(sf::st_is(
    x = outDat$dat$geometry,
    type = c("POLYGON", "MULTIPOLYGON")
  )),
  expected = nrow(outDat$dat))
})

## test that the output 'oldNames' part of list is the same as names of the
# input 'dat' data.frame
test_that(
  "names of input data.frame are the same as the names in the
          'userColNames' element of the output list", {
    expect_equal(object = sum(outDat$userColNames %in% names(dat)),
                 expect = 5)
  })
