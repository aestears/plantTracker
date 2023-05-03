#' Tracks genets through time for multiple species and sites
#'
#' @description This function tracks individual organisms in mapped quadrats
#' through time to generate a demographic dataset that includes survival
#' and growth for each individual.
#'
#' @details
#' This is a wrapper function that applies \code{\link{assign}} across multiple
#' species, quadrats, and sites. For each species and quadrat, [trackSpp()]
#' loads a spatially referenced data.frame ('dat'), and then uses the
#' \code{\link{groupByGenet}} function to assign genetIDs to polygons (if
#' 'clonal' = TRUE) such that polygons that form the same genet have the same
#' genetID. A buffer of a distance defined by 'buff' is applied around each
#' genet polygon. Then, the spatial data for each genet from the current year
#' (year `t`) is compared to individuals in the next year (year `t+1`). Then
#' [trackSpp()] calculates the amount of overlapping area between polygons
#' of each year `t` genet and polygons of each year `t+1` genet (using
#' \code{\link[sf]{st_intersection}}). If there is unambiguous overlap between a
#' 'parent' genet from year `t` and a 'child' genet from year `t+1`, then that
#' 'child' gets the same identifying trackID as the parent. If there is a 'tie,'
#' where more than one parent overlaps the same child or more than one child
#' overlaps the same parent, the parent-child pair with the greatest amount of
#' overlap receives the same trackID. Polygons in year `t+1` that do not have a
#' parent are given new trackIDs and are identified as new recruits. If dormancy
#' is not allowed, then polygons in year `t` that do not have child polygons get
#' a '0' in the 'survival' column. If dormancy is allowed, parent polygons
#' without child polygons are stored as 'ghosts' and are then compared to data
#' from year `t+1+i` to find potential child polygons, where `i`='dorm'
#' argument. For a more detailed description of the [trackSpp()] function, see
#' the vignette:
#' \code{vignette("Using_the_plantTracker_trackSpp_function",
#' package = "plantTracker")}
#'
#' @param dat An sf data.frame of the same format as
#' \code{\link{grasslandData}}. It must have columns that contain ...
#' #' * a unique identification for each research site in character format
#' with no NAs (the default column name is "Site")
#' * species name in character format with no NAs (the default column
#' name is "Species")
#' * unique quadrat identifier in character format with no NAs (the default
#' column name is "Quad")
#' *  year of data collection in integer format with no NAs (the
#' default column name is "Year")
#' * an s.f 'geometry' column that contains a polygon or multipolygon data type
#' for each individual observation (the default column name is "geometry")
#' and an s.f 'geometry' column
#'
#' This function will add columns called
#' "basalArea_ramet", "trackID", "age", "size_tplus1", "recruit," "nearEdge,"
#' and "survives_tplus1", so 'dat' should not contain columns with these names.
#' @param inv A named list of the same format as
#' \code{\link{grasslandInventory}}. The name of each element of the list is a
#' quadrat name in 'dat', and the contents of that list element is a numeric
#' vector of all of the years in which that quadrat (or other unique spatial
#' area) was sampled. Make sure this is the years the quadrat was actually
#' sampled, not just the years that have data in the 'dat' argument! This
#' argument allows the function to differentiate between years when the quadrat
#' wasn't sampled and years when there just weren't any individuals of a species
#' present in that quadrat.
#' @param dorm A numeric vector of length 1, indicating the number of years an
#' individual of these species is allowed to go dormant, i.e. be absent from the
#' map but be considered the same individual when it reappears. This must be an
#' integer greater than or equal to 0. *OR* dorm can be a data.frame with
#' the columns "Species" and "dorm". This data.frame must have a row for each
#' unique species present in 'dat', with species name as a character string in
#' the "Species" column, and a numeric value greater than or equal to 0 in the
#' 'dorm' column that indicates the number of years each species is allowed to
#' go dormant.
#' @param buff A numeric vector of length 1 that is greater than or equal to
#' zero, indicating how far (in the same units as spatial values in 'dat') a
#' polygon can move from year `i` to year `i+1` and still be
#' considered the same individual. *OR* buff can be a data.frame with
#' the columns "Species" and "buff". This data.frame must have a row for each
#' unique species present in 'dat', with species name as a character string in
#' the "Species" column, and a numeric value in the 'buff' column specifying the
#' 'buff' argument for each species.
#' @param buffGenet A numeric vector of length 1 that is greater than or equal
#' to zero, indicating how close (in the same units as spatial values in 'dat')
#' polygons must be to one another in the same year to be grouped as a genet
#' (if 'clonal' argument = TRUE). *OR* buffGenet can be a data.frame with
#' the columns "Species" and "buffGenet". This data.frame must have a row for
#' each unique species present in 'dat', with species name as a character string
#' in the "Species" column, and a numeric value greater than or equal to 0 in
#' the 'buffGenet' specifying the 'buffGenet' argument for each species. This
#' argument is passed to the \code{\link{groupByGenet}} function, which is used
#' inside the \code{\link{assign}} function.
#' @param clonal A logical vector of length 1, indicating whether a
#' species is allowed to be clonal or not (i.e. if multiple polygons (ramets)
#' can be grouped as one individual (genet)). If clonal = TRUE, the species is
#' allowed to be clonal, and if clonal = FALSE, the species is not allowed to
#' be clonal. *OR* clonal can be a data.frame with the columns "Species" and
#' "clonal". This data.frame must have a row for each unique species present in
#' 'dat', with species name as a character string in the "Species" column, and a
#' logical value in the 'clonal' specifying the 'clonal' argument for
#' each species.
#' @param species An optional character string argument. Indicates
#' the name of the column in 'dat' that contains species name data. It is
#' unnecessary to provide a value for this argument if the column name is
#' "Species" (default value is 'Species').
#' @param site An optional character string argument. Indicates
#' the name of the column in 'dat' that contains site name data. It is
#' unnecessary to provide a value for this argument if the column name is
#' "Site" (default value is 'Site').
#' @param quad An optional character string argument. Indicates
#' the name of the column in 'dat' that contains quadrat name data. It is
#' unnecessary to provide a value for this argument if the column name is
#' "Quad" (default is 'Quad').
#' @param year An optional character string argument. Indicates
#' the name of the column in 'dat' that contains data for year of sampling. It
#' is unnecessary to provide a value for this argument if the column name is
#' "Year" (default is 'Year').
#' @param geometry An optional character string argument. Indicates
#' the name of the column in 'dat' that contains sf geometry data. It is
#' unnecessary to provide a value for this argument if the column name is
#' "geometry" (default is 'geometry').
#' @param aggByGenet A logical argument that determines whether the output
#' of [trackSpp()] will be aggregated by genet. If the value is TRUE
#' (the default), then each unique trackID (or genet) in each year will be
#' represented by only one row in the output data.frame. This prepares the
#' dataset for most demographic analyses. If the value is FALSE, then each
#' unique trackID in each year may be represented by multiple rows in the data
#' (each ramet gets a row). Note that if the value is TRUE, then some columns
#' present in the input data.frame 'dat' will be dropped. If you do not wish
#' this to happen, then you can aggregate the data.frame to genet by hand.
#' @param printMessages A logical argument that determines whether this
#' function returns messages about genet aggregation, as well as messages
#' indicating which year is the last year of sampling in each quadrat and which
#' year(s) come before a gap in sampling that exceeds the 'dorm' argument (and
#' thus which years of data have an 'NA' for "survives_tplus1" and
#' "size_tplus1"). If printMessages = TRUE (the default), then messages are
#' printed. If printMessages = FALSE, messages are not printed.
#' @param flagSuspects A logical argument of length 1, indicating whether
#' observations that are 'suspect' will be flagged. The default is
#' `flagSuspects = FALSE`. If `flagSuspects = TRUE`, then a column called
#' 'Suspect' is added to the output data.frame. Any suspect observations get a
#' 'TRUE' in the 'Suspect' column, while non-suspect observations receive a
#' 'FALSE'. There are two ways that an observation can be classified as
#' 'suspect'. First, if two consecutive observations have the same trackID, but
#' the basal area of the observation in year `t+1` is less that a certain
#' percentage (defined by the `shrink` arg.) of the basal area of the
#' observation in year `t`, it is possible that the observation in year `t+1` is
#' a new recruit and not the same individual. The second way an observation can
#' be classified as 'suspect' is if it is very small before going dormant. It is
#' unlikely that a very small individual will survive dormancy, so it is
#' possible that the function has mistakenly given a survival value of '1' to
#' this individual. A 'very small individual' is any observation with an area
#' below a certain percentile (specified by 'dormSize') of the size{} distribution f
#' or this species, which is generated using all of the size data for this
#' species in 'dat'.
#' @param shrink A single numeric value. This value is only used when
#' `flagSuspects = TRUE`. When two consecutive observations have the same
#' trackID, and the ratio of size `t+1` to size `t` is smaller than the value of
#' `shrink`, the observation in year `t` gets a 'TRUE' in the 'Suspect' column.
#' For example, `shrink = 0.2`, and an individual that the tracking function has
#' identified as 'BOUGRA_1992_5' has an area of 9 cm^2 in year `t` and an area
#' of 1.35 cm^2 in year `t+1`. The ratio of size `t+1` to size `t` is
#' 1.35/9 = 0.15, which is smaller than the cutoff specified by `shrink`, so the
#' observation of BOUGRA_1992_5' in year `t` gets a 'TRUE' in the 'Suspect'
#' column. The default value is `shrink = 0.10`.
#' @param dormSize A single numeric value. This value is only used when
#' `flagSuspects = TRUE` and `dorm` is greater than or equal to 1. An individual
#' is flagged as 'suspect'
#' if it 'goes dormant' and has a size that is less than or equal to the
#' percentile of the size distribution for this species that is designated by
#' `dormSize`. For example `dormSize = 0.05`, and an individual has a basal area
#' of 0.5 cm^2. The 5th percentile of the distribution of size for this species,
#' which is made using the mean and standard deviation of all observations in
#' 'dat' for the species in question, is 0.6 cm^2. This individual does not have
#' any overlaps in the next year (year `t+1`), but does have an overlap in year
#' `t+2`. However, because the basal area of this observation is smaller than the
#' 5th percentile of size for this species, the observation in year t will get a
#' 'TRUE' in the 'Suspect' column. It is possible that the tracking function has
#' mistakenly assigned a '1' for survival in year `t`, because it is unlikely
#' that this individual is large enough to survive dormancy. The default value
#' is `dormSize = .05`.
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return An sf data.frame with the same columns as 'dat,' but with the
#' following additional columns:
#'
#' \item{trackID}{A unique value for each individual genet, consisting of the
#' 6-letter species code, the year in which this individual was recruited, and a
#' unique index number, all separated by a "_".}
#' \item{age}{An integer indicating the age of this individual in year `t`.
#' Values of NA indicate that an accurate age cannot be calculated because this
#' individual was observed either in the first year of sampling or in a year
#' following a gap in sampling, so the exact year of recruitment is not known.}
#' \item{size_tplus1}{The  size of this **genet** in year `t+1`, in the same
#' units as the 'area' column in 'dat'.}
#' \item{recruit}{A Boolean integer indicating whether this individual is a new
#' recruit in year `t` (1), or existed in a previous year (0). Values of NA
#' indicate that this individual was observed either in the first year of
#' sampling or in a year following a gap in sampling, so it is not possible to
#' accurately determine whether or not it is a new recruit in year `t`.}
#' \item{survives_tplus1}{A Boolean integer indicating whether this individual
#' survived (1), or died (0) in year `t+1`.}
#' \item{basalArea_genet}{The size of this entire genet in year `t`, in the same
#' units as the 'area' column in 'dat.' If the 'clonal' argument = FALSE, then
#' this number will be identical to the 'basalArea_ramet' column. }
#' \item{basalArea_ramet}{This is only included if 'aggByGenet' = FALSE.
#' This is the size of this ramet in year `t`, in the same units as the 'area'
#' column in 'dat'. If the 'clonal' argument = FALSE , then this number will be
#' identical to the 'basalArea_genet' column.}
#' \item{nearEdge}{A logical value indicating whether this individual is within
#' a buffer (specified by the 'buff' argument) from the edge of the quadrat.}
#'
#' @seealso [assign()], which is used inside the [trackSpp()] function.
#' [trackSpp()] applies the [assign()] function over multiple species and
#' quadrats, and uses the [aggregateByGenet()] function to aggregate the
#' demographic results by genet (if 'aggByGenet' = TRUE). The [assign()]
#' function uses the [groupByGenet()] function to group ramets into genets
#' (if 'clonal' argument = TRUE).
#'
#' @examples
#' dat <- grasslandData[grasslandData$Site == c("AZ") &
#'  grasslandData$Species %in% c("Bouteloua rothrockii",
#'   "Calliandra eriophylla"),]
#' names(dat)[1] <- "speciesName"
#' inv <- grasslandInventory[unique(dat$Quad)]
#' outDat <- trackSpp(dat = dat,
#'  inv = inv,
#'  dorm = 1,
#'  buff = .05,
#'  buffGenet = 0.005,
#'  clonal = data.frame("Species" = unique(dat$speciesName),
#'  "clonal" = c(TRUE,FALSE)),
#'  species = "speciesName",
#'  aggByGenet = TRUE
#'  )
#'
#' @export
#' @import sf
#' @importFrom units drop_units

trackSpp <- function(dat, inv, dorm , buff , buffGenet , clonal,
                     species = "Species",
                     site = "Site",
                     quad = "Quad",
                     year = "Year",
                     geometry = "geometry",
                     aggByGenet = TRUE,
                     printMessages = TRUE,
                     flagSuspects = FALSE,
                     shrink = .1,
                     dormSize = .05,
                     ...) {
  # argument checks ---------------------------------------------------------
  ## arguments

  ## check the 'dat' and 'inv' arguments using the 'checkDat' function
  checkData <- checkDat(dat = dat, inv = inv,
                        species = species,
                        site = site,
                        quad = quad,
                        year = year,
                        geometry = geometry,
                        reformatDat = TRUE)

  dat <- checkData$dat
  inv <- checkData$inv
  usrNames <- checkData$userColNames

  # dorm
  if (missing(dorm)) {
    stop("The 'dorm' argument must have a value.")
  } else if (is.numeric(dorm) & length(dorm == 1)) { ## is the value of dorm a
    # single numeric integer?
    if (dorm < 0 | ## dorm must be greater than or equal to 0
        round(dorm) != dorm | ## dorm must be a whole number
        length(dorm)!=1) { ## dorm must be a vector of length 1
      stop("If 'dorm' is not specified for every species, it must be a single
      numeric value that is a whole number greater than or equal to 0")
    }
    ## make the 'dorm' arg. into a data.frame to make things easier
    dorm <- data.frame("Species" = unique(dat$Species), "dorm" = dorm)
  } else if (is.data.frame(dorm)) {
    if (sum(!names(dorm) %in% c("Species", "dorm")) == 0) {
      if(sum(!unique(dat$Species) %in% dorm$Species) > 0 | ## dorm must have
         # data for all species
         sum(is.na(dat$dorm)) > 0 | ## can't have NA values in dorm
         !is.numeric(dorm$dorm) | ## can't have non-numeric values for
         # dorm$dorm
         sum(dorm$dorm < 0) > 0 | ## can't be less than 0
         round(dorm$dorm) != dorm$dorm ## must be whole numbers
      ) {
        stop("If the 'dorm' argument is specified by species, it must be a
        data.frame that includes a 'Species' column with a row for every species
        in 'dat', and a 'dorm' column that contains positive, whole number
        values for each species with no NAs.")
      }
    } else {
      stop("If the 'dorm' argument is specifed by species, the column names must
      be 'Species' and 'dorm'")
    }
  } else {
    stop("The 'dorm' argument must be either a single numeric value that is a
    whole number greater than or equal to 0, OR a data.frame that has a
    'Species' column with values for each species in 'dat', and a 'dorm' column
    with numeric, positive whole number values for each species.")
  }

  #buff
  if(missing(buff)) {
    stop("The 'buff' argument must have a value.")
  } else {
    if (is.numeric(buff) & length(buff) == 1) { ## is the value of buff a single
      # numeric value?
      if (buff < 0 | ## buff must be greater than or equal to 0
          buff > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buff must
          # not be larger than the dimensions of the quadrat
          length(buff)!=1) { ## buff must be a vector of length 1
        stop("If 'buff' is not specified for every species, it must be a single
        numeric value that is greater than or equal to 0")
      }
      ## make the 'buff' arg. into a data.frame to make things easier
      buff <- data.frame("Species" = unique(dat$Species), "buff" = buff)
    } else if (is.data.frame(buff)) {
      if (sum(!names(buff) %in% c("Species", "buff")) == 0) {
        if(sum(!unique(dat$Species) %in% buff$Species) > 0 | ## buff must have
           # data for all species
           sum(is.na(dat$buff)) > 0 | ## can't have NA values in buff
           !is.numeric(buff$buff) | ## can't have non-numeric values for
           # buff$buff
           sum(buff$buff < 0) > 0 | ## can't be less than 0
           round(buff$buff) != buff$buff ## must be whole numbers
        ) {
          stop("If the 'buff' argument is specified by species, it must be a
          data.frame that includes a 'Species' column with a row for every
          species in 'dat', and a 'buff' column that contains positive, numeric
          values for each species with no NAs.")
        }
      } else {
        stop("If the 'buff' argument is specifed by species, the column names
        must be 'Species' and 'buff'")
      }
    } else {
      stop("The 'buff' argument must be either a single numeric value that is
      greater than or equal to 0, OR a data.frame that has a 'Species' column
      with values for each species in 'dat', and a 'buff' column with numeric
      values for each species.")
    }
  }

  #clonal
  ## check clonal argument
  if (missing(clonal)) {
    stop("The 'clonal' argument must have a value.")
  } else {
    if (is.logical(clonal) & length(clonal == 1)) { ## is the value of clonal a
      # single logical vector?
      if (clonal != TRUE & clonal != FALSE | ## clonal must be either 0 or 1
          !is.logical(clonal) | ## clonal must be logical
          length(clonal)!=1){ ## clonal must be a vector of length = 1
        stop("If 'clonal' is not specified for every species, it must be a
        single logical value that is either FALSE or TRUE.")
      }
      ## make the 'clonal' arg. into a data.frame to make things easier
      clonal <- data.frame("Species" = unique(dat$Species), "clonal" = clonal)
    } else if (is.data.frame(clonal)) {
      if (sum(!names(clonal) %in% c("Species", "clonal")) == 0) {
        if (sum(!unique(dat$Species) %in% clonal$Species) > 0 | ## clonal
           # must have data for all species present in 'dat'
           length(unique(clonal$Species)) != nrow(clonal) |## clonal cannot have
           # more than one value for each species
           sum(is.na(clonal$clonal)) > 0 | ## can't have NA values in clonal
           (!is.numeric(clonal$clonal) & !is.logical(clonal$clonal)) | ## can't
           # have non-numeric values for clonal$clonal
           sum((clonal$clonal != TRUE & clonal$clonal != FALSE)) > 0
           ## clonal values must be either 0 or 1
        ) {
          stop("If the 'clonal' argument is specified by species, it must be a
          data.frame that includes a 'Species' column with a row for every
          species in 'dat', and a 'clonal' column that contains logical values
          of either FALSE or TRUE for each species with no NAs. There cannot be
          multiple rows for the same species.")
        }
      } else {
        stop("If the 'clonal' argument is specifed by species, the column names
        must be 'Species' and 'clonal'")
      }
    } else {
      stop("The 'clonal' argument must be either a single logical value that is
either TRUE or FALSE, OR a data.frame that has a 'Species' column with
values for each species in 'dat', and a 'clonal' column that contains logical
values of either FALSE or TRUE for each species with no NAs.")
    }
  }

  #buffGenet ## either a single value (applied to all spp.) or a data.frame with
  # the same number of rows as the number of species in dat that indicates how
  # close together ramets must be to be considered the same genet (in the same
  # units as distances in dat). If multiple values, is subset by spp. before
  # being passed to assign() (and then passed to groupByGenet())
  ## check buffGenet argument
  ## is the clonal d.f contain at least one value that == TRUE? If so, then you
  # need a buffGenet argument.
  # does the clonal arg. have have any arguments that are 1?
  if (sum(clonal$clonal) > 0) {
    if (missing(buffGenet)) {
      stop("The 'buffGenet' argument must have a value.")
    } else {
      if (is.numeric(buffGenet) & length(buffGenet == 1)) { ## is the value of
        # buffGenet a single numeric?
        if (buffGenet < 0 | ## buffGenet must be greater than or equal to 0
            buffGenet > max(st_bbox(dat)[c("xmax", "ymax")]) | ## buffGenet
            # must not be larger than the dimensions of the quadrat
            length(buffGenet)!=1) { ## buffGenet must be a vector of length 1
          stop("If 'buffGenet' is not specified for every species, it must be
            a single numeric value that is greater than or equal to 0")
        }
        ## make the 'buffGenet' arg. into a data.frame to make things easier
        buffGenet <- data.frame("Species" = unique(dat$Species),
                                "buffGenet" = buffGenet)
      } else if (is.data.frame(buffGenet)) {
        if (sum(!names(buffGenet) %in% c("Species", "buffGenet")) == 0) {
          ## if the buffGenet d.f contains the correct column names
          if (sum(!unique(dat$Species) %in% buffGenet$Species) > 0 |
              ##if buffGenet does not have data for all species
              sum(is.na(dat$buffGenet)) > 0 | ## can't have NA values in
              # buffGenet
              !is.numeric(buffGenet$buffGenet) | ## can't have non-numeric
              # values for buffGenet$buffGenet
              sum(buffGenet$buffGenet < 0) > 0 ## can't be less than 0
          ) {
            stop("If the 'buffGenet' argument is specified by species, it must
              be a data.frame that includes a 'Species' column with a row for
              every species in 'dat', and a 'buffGenet' column that contains
              positive, numeric values for each species with no NAs.")
          }
        } else {
          stop("If the 'buffGenet' argument is specifed by species, the column
            names must be 'Species' and 'buffGenet'")
        }
      } else {
        stop("The 'buffGenet' argument must be either a single numeric value
          that is greater than or equal to 0, OR a data.frame that has a
          'Species' column with values for each species in 'dat', and a
          'buffGenet' column with numeric values for each species.")
      }
    }
    #check aggByGenet
    if (!is.logical(aggByGenet)) {
      stop("The 'aggByGenet' argument must be a logical value. TRUE
        means that every row in the output of trackSpp() represents a unique
        genetic individual (genet) in a given year. FALSE means that every row
        in the output of trackSpp() represents a unique stem (ramet) in a given
        year.")
    }
  } else {
    buffGenet <- data.frame("Species" = unique(dat$Species),
                            "buffGenet" = NA)
    aggByGenet <- FALSE
  }

  #check printMessages
  if (!is.logical(printMessages)) {
    stop("The 'printMessages' argument must be a logical value.")
  }

  ## need to check the 'flagSuspects', 'shrink', and 'dormSize' arguments
  ## check 'flagSuspects' argument
  if (!is.logical(flagSuspects) | ## flagSuspects must be logical
      length(flagSuspects)!=1){ ## flagSuspects must be a vector of length = 1
    stop("'flagSuspects' must be a single logical value that is either
           FALSE or TRUE.")
  }

  ## check 'shrink' argument
  if (is.numeric(shrink) == FALSE | ## is shrink a numeric argument?
      length(shrink) > 1 ## is shrink longer than one value?
  ) {
    stop("'shrink' must be a single numeric value")
  }

  ## check 'dormSize' argument
  if (is.numeric(dormSize) == FALSE | ## is dormSize a numeric argument?
      length(dormSize) > 1 ## is dormSize longer than one value?
  ) {
    stop("'dormSize' must be a single numeric value")
  }

  # work --------------------------------------------------------------------
  ## first, break the 'dat' d.f into two pieces, one with columns we need, and
  # another with columns that are 'extra'--will rejoin at the end of the
  # function
  ## assign an arbitary index number so we can re-join the data at the end
  dat$indexStore <- c(1:nrow(dat))
  ## put the 'extra' data into a separate d.f to 'store'
  datStore <- dat[, !names(dat) %in% c("Site", "Quad", "Year", "Species",
                                       "geometry")]
  datStore <- sf::st_drop_geometry(datStore)
  ## put the data we actually need into the 'dat' object
  dat <- dat[,names(dat) %in% c("Site", "Quad", "Year", "Species",
                                "geometry", "indexStore")]

  ## get the basal area for each observation
  dat$basalArea_ramet <- sf::st_area(dat)
  units(dat$basalArea_ramet) <- NULL

  ## get the site(s)
  for(i in unique(dat$Site)) { ## i = the name of the site
    if(printMessages==TRUE){
      message(paste0("Site: ",i, "\n"))
    }
    ## get the quadrats w/in that site
    for (j in unique(dat[dat$Site==i,]$Quad)) { ## j = the name of the quad
      if(printMessages==TRUE) {
        message(paste0("-- Quadrat: ",j, "\n"))
      }
      ## get the quadratInventory data for this quad
      if (is.list(inv)==TRUE) { ## if there is inv data for >1 quadrat
        invQuad <- inv[[j]]
      } else if (is.vector(inv)) { ## if there is inv data for only 1 quadrat
        invQuad <- inv
      }
      ## get the boundary box for this quadrat (for calculating nearEdge
      # later on)
      ## make a boundary box that is within the 'buff' argument of the quad
      ## first, define the 'buff' zone (the average of 'buff' for each species)
      buffAvg <- mean(buff$buff )
      buffEdgeOutside <- sf::st_as_sfc(sf::st_bbox(dat[dat$Site==i,]))
      buffEdgeInside <- sf::st_as_sfc(sf::st_bbox(dat[dat$Site==i,]) +
                                        c(buffAvg, buffAvg,-buffAvg, -buffAvg))
      buffEdge <-  sf::st_difference(buffEdgeOutside, buffEdgeInside)
      ## get the species w/in that quad
      for (k in unique(dat[dat$Site==i & dat$Quad==j,]$Species)) {
        ## k = the name of the species
        ## get the data for this site/quad/species
        datSp <- dat[dat$Site == i &
                       dat$Quad == j &
                       dat$Species == k,]

        ## get dorm value
        dormK <- dorm[dorm$Species==k,"dorm"]

        ## get clonal value
        clonalK <- clonal[clonal$Species==k,"clonal"]

        ## get buff value
        buffK <- buff[buff$Species==k,"buff"]

        ## get buffGenet value
        if (sum(is.na(buffGenet$buffGenet)) == 0) {
          buffGenetK <- buffGenet[buffGenet$Species==k,"buffGenet"]
        } else {
          buffGenetK <- NA
        }

        ## put this dataset into the 'assign' function
        datOut <- assign(dat = datSp,
                         inv = invQuad,
                         dorm = dormK,
                         buff = buffK,
                         buffGenet = buffGenetK,
                         clonal = clonalK,
                         flagSuspects = flagSuspects,
                         shrink = shrink,
                         dormSize = dormSize,
                         inheritsFromTrackSpp = TRUE,
                         nearEdgeBox = buffEdge
        )
        ## see if the output d.f exists yet (trackSpOut)
        ## if it does exist, then add datOut for the current spp. to the output
        if (i == unique(dat$Site)[1] & ## if this is the first site
            ## if this is the first quad in the first site
            j == unique(dat[dat$Site==i,]$Quad)[1] &
            ## if this is the first species in the first quad in the first site
            k == unique(dat[dat$Site==i & dat$Quad==j,]$Species)[1]
            ) {
          trackSppOut <- datOut
        } else {
          ## if not, rbind the datOut to trackSppOut
          trackSppOut <- rbind(trackSppOut, datOut)
        }
        ## print the name of the species that was just finished
        if (k == unique(dat[dat$Site==i & dat$Quad==j,]$Species)[1]) {
          if(printMessages==TRUE) {
            message(paste0("---- Species: ",k))
          }
        } else if (k == tail(unique(dat[dat$Site==i & dat$Quad==j,]$Species),
                             n = 1)) {
          if(printMessages==TRUE) {
            message(paste0("; ",k, "\n"))
          }
        } else {
          if(printMessages==TRUE) {
            message(paste0("; ",k))
          }
        }
        if (printMessages == TRUE) {
          ## find years that exceed the 'dorm' gap
          invComp <- data.frame(inv = c(NA, invQuad), invNext = c(invQuad, NA))
          invComp$diff <- invComp$invNext - invComp$inv
          gapYears <- invComp[invComp$diff > (dormK + 1) &
                                is.na(invComp$diff) == FALSE,"inv"]
          if (length(gapYears) > 0) {
            message(paste0("Also Note: Individuals of the species ",
                         unique(datSp$Species)," in year(s) ", gapYears,
                         " have a value of 'NA' in the 'survives_tplus1' and",
                         " 'size_tplus1' columns because ", gapYears," is the"
                         ," last year of sampling in this quadrat before a gap",
                         " that exceeds the 'dorm' argument for this species."))
          }
        }
      }
      ## notify user of last year of sampling (or last year of sampling before a
      # gap)
      if (printMessages == TRUE) {
        message(paste0("Note: Individuals in year ", max(invQuad)," have a value "
                 ,"of 'NA' in the 'survives_tplus1' and 'size_tplus1' columns ",
                   "because ",max(invQuad), " is the last year of sampling in ",
                     "this quadrat."))
      }
    }
  }

  ## rejoin the trackSppOut d.f with the 'extra' data stored in 'datStore'
  trackSppOut <- merge(trackSppOut, datStore, by = "indexStore")

  ## remove the 'indexStore' value
  trackSppOut <- trackSppOut[,names(trackSppOut) != "indexStore"]

  ## aggregate the output by trackID (if aggByGenet == TRUE)
  if (aggByGenet == TRUE) {
    ## aggregate demographic data by trackID/Quad/Year/Site/Species
    trackSppOut <- aggregateByGenet(dat = trackSppOut)

    if (printMessages == TRUE) {
      message(paste0("Note: The output data.frame from this function is shorter",
                   " than your input data.frame because demographic data has",
                   " been aggregated by genet. Because of this, some columns",
                   " that were present in your input data.frame may no longer",
                   " be present. If you don't want the output to be aggregated",
                   " by genet, include the argument 'aggByGenet == FALSE' in",
                   " your call to trackSpp()."))
    }
  }

  ## re-name the appropriate columns in 'trackSppOut' data.frame with the
  # user-provided names of 'dat'
  ## from above, user-provided names are stored in 'usrNames'
  ## make a vector of default column names
  defaultNames <- c("Species", "Site", "Quad", "Year",  "geometry")

  ## reset the names for the columns that we changed to 'default' values
  names(trackSppOut)[match(defaultNames, names(trackSppOut))] <-
    usrNames

  ## remove the '_USER' from the 'extra' column names
  names(trackSppOut) <- gsub(names(trackSppOut),
                             pattern = "_USER", replacement = "")

  # output ------------------------------------------------------------------
  return(trackSppOut)
}

# Testing -----------------------------------------------------------------
# dat <- grasslandData[grasslandData$Site == "CO"
#                      & grasslandData$Quad %in% c("unun_11","ungz_5a"),]
#                      #& grasslandData$Species == "Bouteloua gracilis",]
# names(dat)[1]<- "Species_Name"
# names(dat)[8] <- "location"
# #dat <- dat[dat$location != "ungz_5a",]
# #dat <- dat[,c(1:6,8:13)]
# inv <- grasslandInventory
# #inv <- inv[1:5]
# dorm <- 1
# buff <- .05
# buffGenet <- 0.005
# clonal <- data.frame(Species = unique(dat$Species),
#                      clonal = c(TRUE))
#
# testOut <- trackSpp(dat = dat, inv = inv, dorm = dorm, buff = buff,
#                     buffGenet = buffGenet,
#                     clonal = clonal , species = "Species_Name",
#                     quad = "location", printMessages = FALSE,
#                     flagSuspects = TRUE)


### AES make an example in the documentation that specifies all args as numeric,
# and another example where they specify all four arguments as data.frames

### be careful about warnings: because people freak out about them

###AES make fake errors in data to make sure that the argument checks are
# working correctly  (i.e. a random NA in a column, character for dorm, etc.)

###AES maybe also try to practice on larger subset of data.frame
## find years that exceed the 'dorm' gap
# invComp <- data.frame(inv = c(NA, invQuad), invNext = c(invQuad, NA))
# invComp$diff <- invComp$invNext - invComp$inv
# gapYears <- invComp[invComp$diff > (dorm + 1) &
#                       is.na(invComp$diff) == FALSE,"inv"]
# if (length(gapYears) > 0) {
#   print(paste0("Also Note: Individuals in year(s) ", gapYears," have a",
#                " value of 'NA' in the 'survives_tplus1' and",
#                " 'size_tplus1' columns because ", gapYears," is the last"
#                , " year of sampling in this quadrat before a gap that
#                exceeds the 'dorm' argument."))
# }
