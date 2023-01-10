#' Calculates lambda, the population growth rate, for each species in a quadrat
#' based on changes in basal cover.
#'
#' @description This function calculates the population growth rate (lambda) for
#' every species in a quadrat. This value is the ratio of basal area or number
#' of individuals in the next year to basal area or number of individuals in the
#' current year (basal area in year t+1/ basal area in year t). A lambda value
#' greater than 1 indicates a population is growing, while a value less than 1
#' indicates population decline. Lambda is 'infinity' when the basal area or
#' number of individuals in year t is 0, and is NA when basal area or number of
#' individuals in year t is zero (i.e. when there are no plants present in year
#' t). Note that a lambda value is calculated between of the years when a
#' quadrat was sampled, even if there is a gap in sampling. For example, a
#' quadrat is sampled in 1998, 1999, 2001, and 2002 (but skipped in 2000). A
#' lambda value will be calculated for 1998-1999 and 2001-2002, which is a
#' transition from year `t` to year `t+1`. However, a lambda value is calculated
#' in the same manner for 1999-2001, which is actually a transition from year
#' `t` to year `t+2`. You can remove these values by subsetting the data.frame
#' returned by `getLambda()` for rows when "Year_tplus1"- "Year_t" is equal
#' to 1.
#'
#' @param dat An sf data.frame in which each row represents a unique polygon
#' (either a genet or a ramet) in a unique site/quadrat/year combination. A
#' data.frame returned by \code{\link{trackSpp}} can be put directly into this
#' function. However, it is not necessary for 'dat' to have demographic data or
#' unique identifiers (i.e. 'trackIDs') assigned. 'dat' must have columns that
#' contain a unique identification for each research site (default name is
#' "Site"), species name (default name is "Species"), quadrat identifier
#' (default name is "Quad"), year of data collection (default name is "Year"),
#' and an s.f 'geometry' column that contains a polygon or multipolygon data
#' type for each individual observation.
#' @param inv The name of each element of the list is a
#' quadrat name in 'dat', and the contents of that list element is a numeric
#' vector of all of the years in which that quadrat (or other unique spatial
#' area) was sampled. Make sure this is the years the quadrat was actually
#' sampled, not just the years that have data in the 'dat' argument! This
#' argument allows the function to differentiate between years when the quadrat
#' wasn't sampled and years when there just weren't any individuals of a species
#' present in that quadrat.
#' @param method A single character argument that determines the method for
#' calculating lambda. The default value is "area", which means that lambda is
#' calculated by comparing total basal area for a given species in year t+1 to
#' year t. If 'method' = "count", then lambda is calculated by comparing total
#' number of individuals for a given species in year t+1 to year t. If each
#' individual in 'dat' is mapped as a point, it is best to use
#' 'method' = "count".
#' @param species An optional character string argument. Indicates
#' the name of the column in 'dat' that contains species name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Species" (default value is 'Species').
#' @param quad An optional character string argument. Indicates
#' the name of the column in 'dat' that contains quadrat name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Quad" (default is 'Quad').
#' @param site An optional character string argument. Indicates
#' the name of the column in 'dat' that contains site name data. It is
#' unnecessary to include a value for this argument if the column name is
#' "Site" (default value is 'Site').
#' @param year An optional character string argument. Indicates
#' the name of the column in 'dat' that contains data for year of sampling. It
#' is unnecessary to include a value for this argument if the column name is
#' "Year" (default is 'Year').
#' @param geometry An optional character string argument. Indicates
#' the name of the column in 'dat' that contains sf geometry data. It is
#' unnecessary to include a value for this argument if the column name is
#' "geometry" (default is 'geometry').
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return This function returns a data.frame with columns containing site,
#' quadrat, and species data, as well as the following columns:
#' \describe{
#' \item{Year_t}{the 'current' year}
#' \item{absolute_basalArea_t}{basal area (in the same units as the spatial
#' elements of 'dat') for this species in this quadrat in year 't'}
#' \item{Year_tplus1}{the 'next' year}
#' \item{absolute_basalArea_tplus1}{basal area (in the same units as the spatial
#' elements of 'dat') for this species in this quadrat in year 't+1'}
#' \item{lambda}{The population growth rate for this species in this quadrat
#' from year t to year t+1.}
#' }
#'
#' @seealso [getBasalAreas()], used internally in this function, which
#' calculates absolute and relative basal areas for each species in each quadrat
#' for each year of sampling.
#'
#' @import sf
#'
#' @export
#'
#' @examples
#' dat <- grasslandData[grasslandData$Site == "AZ" &
#'  grasslandData$Year %in% c(1922:1925),]
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
#' getLambda(dat = outDat, inv = inv, method = "area",
#' species = "speciesName")
getLambda <- function(dat,
                      inv,
                      method = "area",
                      species = "Species",
                      quad = "Quad",
                      site = "Site",
                      year = "Year",
                      geometry = "geometry",
                      ...) {

  # argument checking -------------------------------------------------------
  ## use the 'checkDat' function to make sure that 'dat' is in the
  # correct format
  checked <- checkDat(dat = dat, inv = inv, species = species, site = site,
                      quad = quad, year = year, geometry = geometry,
                      reformatDat = TRUE)
  ## save the checked data to use in the function
  dat <- checked$dat
  usrNames <- checked$userColNames
  ## save the inputed user names (will reassign later to the returned d.f.)
  inv <- checked$inv

  ## check the 'method' argument
  #method
  if (is.na(method) == TRUE) {
    stop("The 'method' argument must have a value. It must be a character vector
    with either the value 'area' (meaning you want to calculate lambda based on
         change in total basal area for each species) or 'count' (meaning you
         want to calculate lambda based on change in number of individuals for
         each species.")
  } else if (is.character(method) & length(method) == 1) { ## must be a
    # character of length one
    ## make sure that the 'method' argument is lowercase
    method <- tolower(method)
    if (method != 'area' & method != 'count') {
      stop("The 'method' argument must have a value of 'area' or 'count'.")
    }
  } else {
    stop("'method' must be a character vector of length one.")
  }

  # work --------------------------------------------------------------------
  if (method == "count") {
    ## get the basal areas for each species in each quadrat in each year

    ## if there is a trackID column, use that to aggregate
    if ("trackID_USER" %in% names(dat)) {
      dat <- sf::st_drop_geometry(dat)
      dat$trackID_USER <- as.factor(dat$trackID_USER)
      datCounts <- aggregate(x = dat[,c("trackID_USER")], by = list(
        "Site" = dat$Site,
        "Species" = dat$Species,
        "Quad" = dat$Quad,
        "Year" = dat$Year
      ),
        FUN = function(x)
          length(unique(x)))
      names(datCounts)[5] <- "Count"
    } else {
      ## otherwise, use each row as a unique individual
      dat <- sf::st_drop_geometry(dat)
      datCounts <- aggregate(x = rownames(dat), by = list(
        "Site" = dat$Site,
        "Species" = dat$Species,
        "Quad" = dat$Quad,
        "Year" = dat$Year
      ),
      FUN = length)
      names(datCounts)[5] <- "Count"
    }

    datLambda <- data.frame(NULL)
    ## calculate lambda (count_t+1/count_t)
    ## loop through each site
    for (i in unique(datCounts$Site)) {
      ## loop through each quad
      for (j in unique(datCounts[datCounts$Site == i,"Quad"])) {
        ## loop through each species
        for (k in unique(datCounts[datCounts$Site == i & datCounts$Quad == j,
                                  "Species"])) {
          ## get the data for each species
          datSpp <- datCounts[datCounts$Site == i &
                               datCounts$Quad == j &
                               datCounts$Species == k,
                             c("Site", "Quad", "Species",
                               "Year", "Count")]
          ## make sure that years are in order
          datSpp <- datSpp[order(datSpp$Year),]
          ## print a message if there is only one year of data (need at least
          # two time points to calculate lambda)
          if (nrow(datSpp) <= 1) {
            message(paste0("A value of lambda for ", k, " in quadrat ", j,
                 " cannot be calculated, since data was only available for one",
                   " year (or this species was not observed in more than one ",
                 " year). Data from at least two time points are required to",
                         " calculate lambda."))
            names(datSpp)[which(names(datSpp) == "Year")] <- "Year_t"
            names(datSpp)[which(names(datSpp) == "Count")] <-
              "Count_t"
            datSpp$Year_tplus1 <- NA
            datSpp$Count_tplus1 <- NA
            datSpp$lambda <- NA

            next
          }

          ## calculate lambda
          names(datSpp)[which(names(datSpp) == "Year")] <- "Year_t"
          names(datSpp)[which(names(datSpp) == "Count")] <-
            "Count_t"
          ## get the 'Year_tplus1' value (the next sequential year of sampling,
          # NOT necessarily the next sequential year)
          datSpp$Year_tplus1 <- c(datSpp$Year_t[2:length(datSpp$Year_t)], NA)
          ## get the 'Yabsolute_basalArea_tplus1' value (the next sequential
          # year of sampling, NOT necessarily the next sequential year)
          datSpp$Count_tplus1 <- c(datSpp$Count_t[2:length(datSpp$Count_t)],NA)
          ## finally calculate lambda itself
          datSpp$lambda <- datSpp$Count_tplus1/datSpp$Count_t

          ## save the output
          if (i == unique(datCounts$Site)[1] &
              j == unique(datCounts[datCounts$Site == i,"Quad"])[1] &
              k == unique(datCounts[datCounts$Site == i & datCounts$Quad == j,
                                    "Species"])[1]) {
            datLambda <- datSpp
          } else {
            datLambda <- rbind(datLambda, datSpp)
          }
        }
      }
    }
  } else if (method == "area") {
    ## get the counts for each species in each quadrat in each year

    datAreas <- getBasalAreas(dat = dat, inv = inv, species = "Species",
                              quad = "Quad", site = "Site", year = "Year",
                              geometry = "geometry")

    datLambda <- data.frame(NULL)
    ## calculate lambda (area_t+1/area_t)
    ## loop through each site
    for (i in unique(datAreas$Site)) {
      ## loop through each quad
      for (j in unique(datAreas[datAreas$Site == i,"Quad"])) {
        ## loop through each species
        for (k in unique(datAreas[datAreas$Site == i & datAreas$Quad == j,
                                  "Species"])) {
          ## get the data for each species
          datSpp <- datAreas[datAreas$Site == i &
                               datAreas$Quad == j &
                               datAreas$Species == k,
                             c("Site", "Quad", "Species",
                               "Year", "absolute_basalArea")]
          ## make sure that years are in order
          datSpp <- datSpp[order(datSpp$Year),]
          ## print a message if there is only one year of data (need at least
          # two time points to calculate lambda)
          if (nrow(datSpp) <= 1) {
            message(paste0("A value of lambda for ", k, " in quadrat ", j,
                 " cannot be calculated, since data was only available for one",
                   " year (or this species was not observed in more than one ",
                   " year). Data from at least two time points are required to",
                         " calculate lambda."))
            names(datSpp)[which(names(datSpp) == "Year")] <- "Year_t"
            names(datSpp)[which(names(datSpp) == "Count")] <-
              "absolute_basalArea_t"
            datSpp$Year_tplus1 <- NA
            datSpp$absolute_basalArea_tplus1 <- NA
            datSpp$lambda <- NA

            next
          }

          ## calculate lambda (basalArea_tplus1 / basalArea_t)
          names(datSpp)[which(names(datSpp) == "Year")] <- "Year_t"
          names(datSpp)[which(names(datSpp) == "absolute_basalArea")] <-
            "absolute_basalArea_t"
          ## get the 'Year_tplus1' value (the next sequential year of sampling,
          # NOT necessarily the next sequential year)
          datSpp$Year_tplus1 <- c(datSpp$Year_t[2:length(datSpp$Year_t)], NA)
          ## get the 'Yabsolute_basalArea_tplus1' value (the next sequential
          # year of sampling, NOT necessarily the next sequential year)
          datSpp$absolute_basalArea_tplus1 <-
            c(datSpp$absolute_basalArea_t[2:length(
              datSpp$absolute_basalArea_t)],NA)
          ## finally calculate lambda itself
          datSpp$lambda <- datSpp$absolute_basalArea_tplus1/
            datSpp$absolute_basalArea_t

          ## save the output
          if (i == unique(datAreas$Site)[1] &
              j == unique(datAreas[datAreas$Site == i,"Quad"])[1] &
              k == unique(datAreas[datAreas$Site == i & datAreas$Quad == j,
                                    "Species"])[1]) {
            datLambda <- datSpp
          } else {
            datLambda <- rbind(datLambda, datSpp)
          }
        }
      }
    }
  }

  ## remove rows for quad/spp/year_t combos that don't have a year_tplus1
  # (have an NA for year_tplus1)
  datLambda <- datLambda[is.na(datLambda$Year_tplus1) == FALSE,]

  ## change the values for lambda when area_t and area_tplus1 are both zero
  # from NaN to NA
  datLambda[is.nan(datLambda$lambda) == TRUE,"lambda"] <- NA

  ## rename the columns to the user-defined columns
  ## from above, user-provided names are stored in 'usrNames'
  defaultNames <- c("Species", "Site", "Quad")
  ## reset the names for the columns that we changed to 'default' values
  names(datLambda)[match(defaultNames, names(datLambda))] <-
    usrNames[1:3]

  ## fix the rownames so they're not all wonky
  rownames(datLambda) <- 1:nrow(datLambda)

  # return ------------------------------------------------------------------
  return(datLambda)
}

# testing -----------------------------------------------------------------
#
# dat <- grasslandData[grasslandData$Site == "CO" & grasslandData$Year %in% c(1998:2002),]
# names(dat)[1] <- "speciesName"
# inv <- grasslandInventory[unique(dat$Quad)]
# outDat <- trackSpp(dat = dat, inv = inv, dorm = 1, buff = .05,buffGenet = 0.005,clonal = data.frame("Species" = unique(dat$speciesName),"clonal" = c(TRUE,FALSE)), species = "speciesName",aggregateByGenet = TRUE)
# getLambda(dat = outDat, inv = inv, species = "speciesName", method = "count")
