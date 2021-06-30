#' aggregateByGenet
#'
#' @param dat
#' @param site
#' @param quad
#' @param species
#' @param year
#' @param trackID
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
aggregateByGenet <-  function(dat,
                              site = "Site",
                              quad = "Quad",
                              species = "Species",
                              year = "Year",
                              trackID = "trackID",
                              ...) {

# argument checking -------------------------------------------------------

#dat ## an sf d.f that contains data returned from the 'trackSpp' function.
# If 'dat' already only has one row for each unique trackID in each unique year
# (i.e. there are no vegetative individuals), then the output of
# aggregateByGenet() will have the same number of rows as 'dat'. If 'dat' in
# some cases has multiple rows (each of which is a ramet) for the same trackID
# in the same year, then the output of aggregateByGenet() will have fewer rows
# that the input 'dat' d.f. 'dat' MUST have columns called 'basalArea_genet',
# 'age', 'recruit', 'survives_tplus1', and 'size_tplus1'. The trackSpp()
# function adds these columns, so if you have made no changes to the d.f that
# was output from trackSpp(), then your d.f should have these columns.

  ## check that 'dat' is an s.f data.frame
  if (sum(!st_is(dat, type = c("POLYGON", "MULTIPOLYGON"))) != 0) {
print("The 'dat' argument must be an sf data.frame with only 'POLYGON' or 'MULTIPOLYGON' geometries.")
  } else {
    ## check that the 'dat' d.f. has the required columns
    if (sum(c('basalArea_genet', 'age', 'recruit','survives_tplus1',
              'size_tplus1') %in% names(dat)) != 5) {
      print("The 'dat' argument must contain columns that are called 'basalArea_genet', 'age', 'recruit', 'survives_tplus1', and 'size_tplus1'.")
    }
  }

  ## check that the exact required columns in 'dat' ('basalArea_genet','age',
  # 'recruit', 'survives_tplus1', and 'size_tplus1') have the correct data types
  ## 'basalArea_genet'
  if (is.numeric(dat$basalArea_genet)==FALSE |
      sum(!dat$basalArea_genet > 0) != 0 ) {
    print("The 'basalArea_genet' column in 'dat' must be a numeric vector with only positive values.")
  } else {
    ## make sure that the values of 'basalArea_genet' aren't extremely large
    if (sum(dat$basalArea_genet > st_bbox(dat)$xmax * st_bbox(dat)$ymax) > 0)
      {
      print("At least one of the values in the 'basalArea_genet' column is extremely large--larger than the area of the entire quadrat! Double-check the accuracy of the input data.")
    }
  }
  ## 'age'
  if (is.numeric(dat$age) == FALSE | ## age must be numeric
      sum(dat$age < 0, na.rm = TRUE) != 0 ## age values must be positive
      ) {
    print("The 'age' column in 'dat' must be a numeric vector with positive values.")
  }
  ## 'recruit'
  if (is.numeric(dat$recruit) == FALSE | ## recruit must be numeric
      sum((dat$recruit == 0 | dat$recruit == 1 | is.na(dat$recruit))) !=
      nrow(dat)
      ## recruit values must be 1/0/NA
  ) {
    print("The 'recruit' column in 'dat' must be a numeric vector with values of only 1, 0, or NA.")
  }
  ## 'survives_tplus1'
  if (is.numeric(dat$survives_tplus1) == FALSE | ## survives_tplus1 must be numeric
      sum((dat$survives_tplus1 == 0 | dat$survives_tplus1 == 1 | is.na(dat$survives_tplus1))) !=
      nrow(dat)
      ## survives_tplus1 values must be 1/0/NA
  ) {
    print("The 'survives_tplus1' column in 'dat' must be a numeric vector with values of only 1, 0, or NA.")
  }
  ## 'size_tplus1'
  if (is.numeric(dat$size_tplus1)==FALSE |
      sum(!dat$size_tplus1 > 0, na.rm = TRUE) != 0
      ) {
    print("The 'size_tplus1' column in 'dat' must be a numeric vector with only positive values (or NA).")
  } else {
    ## make sure that the values of 'size_tplus1' aren't extremely large
    if (sum(dat$size_tplus1 > st_bbox(dat)$xmax * st_bbox(dat)$ymax,
            na.rm = TRUE) > 0)
    {
      print("At least one of the values in the 'size_tplus1' column is extremely large--larger than the area of the entire quadrat! Double-check the accuracy of the input data.")
    }
  }

#site ## the name of the column in 'dat' that contains the values for site

#quad ## the name of the column in 'dat' that contains the data for quadrat

#species ## the name of the column in 'dat' that contains the data for species

#year ## the name of the column in 'dat' that contains the data for Year

#trackID ## the name of the column in 'dat' that contains the data for trackID
  ## put column name args. into a list for basic checks
  newNames <- list("species" = species, "site" = site, "quad" = quad,
                   "year" = year, "trackID" = trackID)
  ## check that each arg. is a character vector
  if (sum(sapply(newNames, is.character)) != 5) { ## if there is one or more
    # elements of the newNames list that is not a character vector
    ## find which elements are not character vectors
    badArgs <- paste("'",names( which(sapply(newNames, is.character) == FALSE)),
                     "'", collapse = ", and ")

    stop(paste0("The argument(s) ", badArgs, " must each contain a single character
string that gives the name(s) of the column(s) in 'dat' that contain the data
for ", badArgs))

  } else { ## if each of the elements of 'newNames' is a character vector
    ## make sure that each of the elements of newNames is present as a column
    # name in 'dat'
    if (sum(unlist(newNames) %in% names(dat)) != 5) { ## if the column names of
      # 'dat' do NOT match the values provided in 'newNames'
      badBadArgs <- paste("'",names(newNames)[which(!unlist(newNames) %in%
                                                      names(dat))],"'",
                          collapse = ", and ")
      stop(paste0("The argument(s) ", badBadArgs, " contain values that are not column
names in 'dat'. These arguments must be character vectors that give the name(s)
of the column(s) in 'dat' that contain the data for ", badBadArgs, ". Check for
spelling errors, or make sure that you have included values for these arguments that give the name of the columns in 'dat' that contain these data types." ))
    }
  }

# work --------------------------------------------------------------------

## make sure the names of 'dat' are what the function expects
  ## re-assign the names of dat to the default column names
  ## make a vector of 'user names'
  usrNames <- unlist(newNames)
  ## make a vector of 'default names'
  defaultNames <- c("Species", "Site", "Quad", "Year",
                    "trackID")


  ## replace the user-provided names in 'dat' with the default names
  names(dat)[match(usrNames, names(dat))] <- defaultNames

  ## aggregate the 'dat' argument by trackID
  datOut <- aggregate(x = dat[,c('basalArea_genet', 'age', 'recruit',
                               'survives_tplus1', 'size_tplus1')],
            by = list("Site" = dat$Site,
                      "Quad" = dat$Quad,
                      "Species" = dat$Species,
                      "trackID" = dat$trackID,
                      "Year" = dat$Year),
            do_union = TRUE,
            FUN = mean)

  ## change the column name sback to what were present in 'dat'
  ## reset the names for the columns that we changed to 'default' values
  names(datOut)[match(defaultNames, names(datOut))] <-
    usrNames
# output ------------------------------------------------------------------
return(datOut)
  }

