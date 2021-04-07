test_that("adding() adds two numeric variables together", {
  x <- 2
  y <- 3
  z <- 3+2
  expect_identical(adding(2,3), z)
})

