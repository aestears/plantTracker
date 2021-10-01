## get example data
dat <- grasslandData[grasslandData$Site == "CO" &
                       grasslandData$Quad == "unun_11" &
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
maps <- drawQuadMap(
  dat = outDat,
  type = "bySpecies",
  addBuffer = FALSE,
  species = "speciesName"
)

## test that output is a list with 'rect' and 'text' elements
test_that("output is a plot", {
  expect_true(object = is.list(maps) &
                (sum(c("rect","text") %in% names(maps)) == 2)
                )
}
)
