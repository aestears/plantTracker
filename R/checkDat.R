#' checkDat
#'
#' @return
#' @export
#'
#' @examples

## check that the dat and inv arguments contain the appropriate values/formats?
checkDat <- function (dat, inv, datNames = c("Species", "Site", "Quad", "Year",
                                             "sp_code_6", "geometry"), ...) {

# arguments ---------------------------------------------------------------
# dat ## an sf d.f that contains all of the raw digitized map data
  # (in grasslandData format) for as many sites and quads as you'd like.

#inv ## a list of the sampling years for each quadrat included in dat (in the
# same format as grasslandInventory).

# datNames ## a character vector indicating the names in dat for each of the
# required columns (Species, Site, Quad, Year, sp_code_6, geometry)

# work --------------------------------------------------------------------

# output ------------------------------------------------------------------

}
