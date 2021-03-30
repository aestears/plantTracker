# function for a generalized remote database call to read in sites
require(httr)
require(sf)
require(mapview)
require(jsonlite)


wkt<- POST(
  "https://datacorral.uwyo.edu/RemoteDb/querytablewktgeometry",
  body = list(Token = "cZlOF3K5w1KkUJ8leW2S1htYSKy0fNjx6nkOtqqZ",
              Database = "Stears_project_DB",
              Table = "lusites"
              #, Filters = list(list(Field = "Site", Min = "2017-01-01", Max = "2017-01-31"))
              # , Filters = list(list(Field = "Site", In = list("AZs")))#,
              #list(Field = "AirTemperature_C", Min = 0, Max = 2))
  ),
  encode = "json"
)

out<- content(wkt, as = "parsed")

# extract req$content 
#cont <- wkt$content

#Convert to char
char <- rawToChar(wkt$content)

#Convert to df 
dat <- as.data.frame(
  jsonlite::fromJSON(char, flatten = TRUE)$data)  
#assign correct column names
names(dat) <- jsonlite::fromJSON(char, flatten = TRUE)$columnInfo[,1]

#make site data spatially referenced
dat <- st_as_sf(dat, geometry = st_as_sfc(dat$ogr_geometry), crs = paste("EPSG:", out$spatialReferences, sep = ""))

