#### Function shp2df: converting shapefile to data frame ####
##'@author Jared Beck
##'@description Function to convert shapefile into data frame of spatial
##'@description data points with a desired projection
##'@param file.path file pathway to shapefile
##'@param layer name of shapefile layer
##'@param projection desired spatial projection, defaults to WGS lat/long coordinates
##'@return Data frame with projected spatial coordinates
shp2df <- function(file.path, layer, projection = "+proj=longlat +datum=WGS84") {
  require(rgdal)
  require(sp)
  require(broom)
  require(dplyr)
  
  shp <- readOGR(dsn = file.path, layer = layer)
  shp.trans <- spTransform(shp, CRS(projection))
  shp.df <- tidy(shp.trans) %>%
    mutate(x = long, y = lat) %>%
    select(-long, -lat)
} ## end shp2df function


#### Function transform.latlong.points: converting data frame of lat/long coords to spatial points ####
##'@author Jared Beck
##'@description Function to transform data frame including lat/long coordinates into data frame of spatial
##'@description data points with a desired project
##'@param df data frame including columns with lat/long coordinates
##'@param projection desired projection for lat/long coordinates, defaults to Eckhart WGS84 global projection
##'@return Data frame with lat/long coordinates and projected spatial coordinates
project_latlong = function(df, projection = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"){
  require(sp)
  require(broom)
  require(dplyr)
  
  ## Standardize naming conventions
  names(df)[names(df) %in% c("lat", "latitude")] = "lat"
  names(df)[names(df) %in% c("lon", "long", "longitude")] = "lon"
  
  ## Convert lat/long coordinates in data frame to spatial points object with WGS84 spatial reference system
  sp = SpatialPoints(cbind(df$lon, df$lat), proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  ## Create spatial points data frame
  spdf = SpatialPointsDataFrame(sp, df)
  
  ## Transform lat/long coordinates into x/y for selected projection
  spdf.trans = spTransform(spdf, CRS(projection))
  
  out = data.frame(spdf.trans) %>%
    mutate(x = coords.x1, y = coords.x2) %>%
    dplyr::select(-coords.x1, -coords.x2, -optional)
  
  return(out)
} ## end project_latlong function



#### Find intersection of two vectors ####
##'@author Jared Beck
##'@description Function to compare elements in two vectors
##'@description Utility extension of intersect and setdiff functions
##'@param a first vector of elements
##'@param b second vector of elements
##'@return list with vectors describing elements in both a and b, in a only, and in b only
compare = function(a, b) {
  require(dplyr)
  
  out = list(in_both = intersect(a,b), in_a_only = setdiff(a,b), in_b_only = setdiff(b,a))
  return(out)
} ## end compare function


################ bootstrap.ci function #############################
#' @description function for bootstrapping single parameters
#' @param x, a vector of numbers to bootstrap
#' @param n.boot, the number of bootstrap resamples, 1000 resamples by default
#' @param level, level of confidence interval, 95% confidence interval by default
#' @param type, the type of summary statistic to bootstrap, mean by default
#' @param type, options include "mean", "harmonic.mean", "median", and "proportion"
#' @details for type == "proportion", function takes binary (0,1) classification of success/failure
bootstrap.ci <- function(x, n.boot = 1000, level = .95, type = "mean") {
  # set upper and lower quantiles
  lwr = (1 - level)/2
  upr = 1 - (1 - level)/2
  
  # create blank vector to house values
  boot.dist <- rep(NA, n.boot)
  
  if(type == "mean") {
    value <- mean(x, na.rm = TRUE)
    for(i in 1:n.boot) {
      boot <- sample(x[!x %in% c(NA)], length(x[!x %in% c(NA)]), replace = TRUE)
      boot.dist[i] <- mean(boot, na.rm = TRUE)}
    
  } else if(type == "harmonic.mean") {
    value <- 1/mean(1/x, na.rm = TRUE)
    for(i in 1:n.boot) {
      boot <- sample(x[!x %in% c(NA)], length(x[!x %in% c(NA)]), replace = TRUE)
      boot.dist[i] <- 1/mean(1/boot, na.rm = TRUE)}
    
  } else if(type == "median") {
    value <- median(x, na.rm = TRUE)
    for(i in 1:n.boot) {
      boot <- sample(x[!x %in% c(NA)], length(x[!x %in% c(NA)]), replace = TRUE)
      boot.dist[i] <- median(boot, na.rm = TRUE)}
    
  } else if(type == "proportion") {
    value <- length(x[x %in% c(1)])/length(x[x %in% c(0,1)])
    for(i in 1:n.boot) {
      boot <- sample(x[!x %in% c(NA)], length(x[!x %in% c(NA)]), replace = TRUE)
      boot.dist[i] <- length(boot[boot %in% c(1)])/length(boot[boot %in% c(0,1)])}
  }
  
  return(list(value = value,
              lower = quantile(boot.dist,probs = lwr),
              upper = quantile(boot.dist,probs = upr)))
} ## end bootstrap.ci function


