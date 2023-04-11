#' Calculates basal area for each species in a quadrat
#'
#' @description This function calculates the total basal area for each species
#' present in a quadrat for each year of data collection. Both the absolute
#' basal area (in the same units of area as the input data.frame) is returned,
#' as well as the percentage of the total occupied basal area in the quadrat
#' that is occupied by a given species ("percent total basal area"). If you'd
#' like to ultimately calculate the population growth rate (lambda) for each
#' species, you can use the \code{\link{getLambda}} function directly, which
#' takes raw species occurrence data (like in 'dat' here) and returns lambda
#' values for each species and quadrat for each transition in the dataset.
#' This function should only be used if the individuals in 'dat' were mapped as
#' polygons that are representative of their actual basal area, i.e. were not
#' mapped as points.
#'
#' @param dat An sf data.frame in which each row represents a unique polygon
#' (either a genet or a ramet) in a unique site/quadrat/year combination. A
#' data.frame returned by \code{\link{trackSpp}} can be put directly into this
#' function. However, it is not necessary for 'dat' to have demographic data or
#' unique identifiers (i.e. 'trackIDs') assigned. 'dat' must have columns that
#' contain...
#' * a unique identification for each research site in character format
#' with no NAs (the default column name is "Site")
#' * species name in character format with no NAs (the default column
#' name is "Species")
#' * unique quadrat identifier in character format with no NAs (the default
#' column name is "Quad")
#' *  year of data collection in integer format with no NAs (the
#' default column name is "Year")
#' * an s.f 'geometry' column that contains a polygon or multipolygon data type
#' for each individual observation (the default column name is "geometry")
#'
#' This function should only be used if
#' the individuals in 'dat' were mapped as polygons that are representative of
#' their actual basal area, i.e. were not mapped as points.
#' @param inv The name of each element of the list is a
#' quadrat name in 'dat', and the contents of that list element is a numeric
#' vector of all of the years in which that quadrat (or other unique spatial
#' area) was sampled. Make sure this is the years the quadrat was actually
#' sampled, not just the years that have data in the 'dat' argument! This
#' argument allows the function to differentiate between years when the quadrat
#' wasn't sampled and years when there just weren't any individuals of a species
#' present in that quadrat.
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
#' @import sf
#' @return [getBasalAreas()] returns a data.frame with the columns "Site",
#' "Quad", "Year" and "Species". No two rows will have the same values for all
#' four of these columns. The column "absolute_basalArea" has the area of the
#' quadrat that is occupied by a species in a given unique site/quadrat/year
#' combination. This is in the same units as the area in 'dat'.
#' "quad_basalArea" gives the combined basal area of all organisms in a given
#' site/quadrat/year. "percent_basalArea" gives the percentage of occupied basal
#' area within a quadrat that is occupied by each species in a given
#' site/quadrat/year. For example, species A has a basal area of 22 cm^2 in
#' quadrat 1 in 2005 ("absolute_basalArea = 22). In 2005, there are 50 cm^2 of
#' quadrat 1 that are occupied by organisms ("quad_basalArea" = 55). 44% of the
#' occupied basal area in quadrat 1 in 2005 is occupied by species A
#' ("percent_basalArea" = 44). There may be an 'NA' in the "percent_basalArea"
#' column if the "quad_basalArea" for that species and year is 0.
#'
#' @seealso [getLambda()], which uses this function to calculate basal areas and
#' ultimately return population growth rates (lambdas) for each species in
#' each quadrat.
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
#' getBasalAreas(dat = outDat, inv = inv,
#' species = "speciesName")
getBasalAreas <- function(dat,
                          inv,
                          species = "Species",
                          quad = "Quad",
                          site = "Site",
                          year = "Year",
                          geometry = "geometry",
                          ...) {
  # argument checks ---------------------------------------------------------
  ## check the 'dat' data.frame and the column names (change if needed)
  checked <- checkDat(dat = dat,
           inv = inv,
           species = species,
           quad = quad,
           site = site,
           year = year,
           geometry = geometry,
           reformatDat = TRUE)
  ## put the checked data into 'dat'
  dat <- checked$dat
  ## put the new names into 'usrNames'
  usrNames <- checked$userColNames
  ## put the checked inv into 'inv'
  inv <- checked$inv

  # work --------------------------------------------------------------------

  ## drop all of the columns except those that you'll need for this fxn
  dat <- dat[,c("Site", "Quad", "Species", "Year", "geometry")]

  ## compare the years of quadrat data against the years the quadrat was sampled
  # so we make sure to assign 0 cover values when appropriate
  ## get only the invs that we need for the quadrats present in 'datArea'
  needInvs <- inv[unique(dat$Quad)]
  ## reformat into a long format data.frame
  needInvsDF <- data.frame(NULL)
  for (i in names(needInvs)) {
    temp <- data.frame("quadInv" = names(needInvs[i]),
                       "yearInv" = needInvs[i])
    names(temp) <- c("quadInv", "yearInv")
    if (nrow(needInvsDF) == 0) {
      needInvsDF <- temp
    } else {
      needInvsDF <- rbind(needInvsDF, temp)
    }
  }

  ## add columns for Site and Species (get from 'dat')
  tempInvDF <- merge(x = unique(
    sf::st_drop_geometry(dat[,c("Site", "Quad", "Species")])),
                y = needInvsDF, by.x = c("Quad"),
                by.y = c("quadInv"), all = TRUE)

  ## sum the area in each year/species/quad/site combo
  datArea <- aggregate(x = sf::st_area(dat), by = list(
    Year = dat$Year,
    Species = dat$Species,
    Quad = dat$Quad,
    Site = dat$Site
  ),
  FUN = sum)

  ## rename the 'x' column to the correct value
  names(datArea)[which(names(datArea) == "x")] <- "absolute_basalArea"

  ## reorder the names of columns
  datArea <- datArea[,c("Site", "Quad", "Species", "Year",
                        "absolute_basalArea")]


  ## join to the datArea data.frame to see if there are any years from the
  # inventory that aren't present in datAreas
  tempArea <- merge(x = datArea, y = tempInvDF,
                by.x = c("Site", "Quad", "Species", "Year"),
                by.y = c("Site","Quad","Species", "yearInv"),
                all = TRUE)
  ## put 0s where they are appropriate (i.e. in years when quad was sampled but
  # no organisms of that species were present)
  tempArea[is.na(tempArea$absolute_basalArea),"absolute_basalArea"] <- 0
  datArea <- tempArea

  ## calculate the percent of total basal area in that year
  # (the basal area of species A / basal area of quadrat occupied by
  # any organisms)
  ## get the total organism area by site/quad/year
  quadBasalArea <- aggregate(x = datArea$absolute_basalArea,
                             by = list(
                               Year = datArea$Year,
                               Quad = datArea$Quad,
                               Site = datArea$Site
                             ),
                             FUN = sum)

  names(quadBasalArea)[which(names(quadBasalArea) == 'x')] <- "quad_basalArea"

  datArea<- merge(x = datArea, y = quadBasalArea,
                  by = c("Site","Quad","Year"))

  ## calculate percentBasalArea
  datArea$percent_basalArea <- (datArea$absolute_basalArea /
                                  datArea$quad_basalArea)*100
  ## change NaN values to NA (happen if the quad_basalArea is 0)
  datArea[is.nan(datArea$percent_basalArea)==TRUE,"percent_basalArea"] <- NA

  ## revert the names of the output data.frame to the names that the user input
  ## re-name the appropriate columns in the output data.frame with the
  # user-provided names of 'dat'
  ## from above, user-provided names are stored in 'usrNames'
  defaultNames <- c("Species", "Site", "Quad", "Year")
  ## reset the names for the columns that we changed to 'default' values
  names(datArea)[match(defaultNames, names(datArea))] <-
    usrNames[1:4]

  # output ------------------------------------------------------------------
  return(datArea)
}
