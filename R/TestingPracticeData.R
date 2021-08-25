# library(tidyverse)
# ## CARBONO example ("doi_10.5061_dryad.51c59zw8n__v2")
# ## get parent file name
# tempWd <- "/Users/Alice/Downloads/doi_10.5061_dryad.51c59zw8n__v2"
#
# ## files we want
# # individualsCARBONOProjectClarkClark19972018.csv
# indivs <- read.csv(paste0(tempWd,"/individualsCARBONOProjectClarkClark19972018.csv")) %>%
#   select(plot_treeid, tree, plot, genspcode, first, second, dist, ang)
# # columns we need: "plot_treeid" (unique tree identifier), "tree" (unique tree identifier w/in plot), "plot" (plot ID), "genspcode" (species code), first" ('left', 'lower' y-coord of the 10 m x 10 m subplot of this tree), "second" ('left', 'lower' x-coord of the 10 m x 10 m subplot of this tree), "dist" (distance of the tree from the 'left,lower' corner of the subplot), 'ang' (the angle from the 'let,lower' corner stake to the tree)
#
# # sizeCARBONOProjectClarkClark19972018.csv
# size <- read.csv(paste0(tempWd,"/sizeCARBONOProjectClarkClark19972018.csv")) %>%
#   select(plot_treeid, dia_year, dia_calc, BAsqm)
# # columns we need: "plot_treeid", "dia_year", "dia_calc" (diameter in mm, missing values are -999), "BAsqm" (basal area in m^2)
#
# # "plotsCARBONOProjectClarkClark19972018.csv"
# plots <- read.csv(paste0(tempWd,"/plotsCARBONOProjectClarkClark19972018.csv")) %>%
#   select(plot, short_angle, long_angle)
# # columns we need: "plot" (plot name), "short_angle" (compass bearing from 0,0 along the 50m side/y-axis), "long_angle" (compass bearing from 0,0 along the 100m side/x-axis)
#
# ## merge datasets into one
# # merge indivs and plots on 'plot'
# temp <- left_join(indivs, plots, by = "plot")
# # calculate the x,y coords for each tree
# temp$want_angle <- NA
# # for plots w/ a long_angle < 90, and a short_ang < 360 but > 270 (P6, P4, P1, L2, A2)
# # eqn: want_angle = 90 = (ang - short_ang)
# temp[temp$plot %in% c("P6", "P4", "P1", "L2", "A2"),"want_angle"] <- 90 - (temp[temp$plot %in% c("P6", "P4", "P1", "L2", "A2"), "ang"] - temp[temp$plot %in% c("P6", "P4", "P1", "L2", "A2"), "short_angle"])
#
# # for plots w/ a long_angle > 90 but < 360
# # eqn: want_angle = long_ang - ang
# temp[temp$plot %in% c("A6", "A5", "A3", "L3", "L1", "A4", "L6", "P3", "P5", "P2", "L5", "L4", "A1"),"want_angle"] <- (temp[temp$plot %in% c("A6", "A5", "A3", "L3", "L1", "A4", "L6", "P3", "P5", "P2", "L5", "L4", "A1"), "long_angle"] - temp[temp$plot %in% c("A6", "A5", "A3", "L3", "L1", "A4", "L6", "P3", "P5", "P2", "L5", "L4", "A1"), "ang"])
#
# # fix wacky want_angles , actually get rid of those trees so we don't have to deal with them.
# temp$bad <- "good"
# temp[temp$want_angle < 0 | temp$want_angle > 90, "bad"] <- "bad"
#
# temp <- temp[temp$bad == "good",]
# temp$bad <- NULL
#
# # put angles in radians
# temp$radians <- temp$want_angle * (pi/180)
# # calculate the x,y coordinates for each tree w/in each subplot
# temp$x_temp <- cos(temp$radians) * temp$dist
# temp$y_temp <- sin(temp$radians) * temp$dist
# # calculate the x,y coordinates for each tree w/in the entire plot
# temp$x <- temp$second + temp$x_temp
# temp$y <- temp$first + temp$y_temp
# # drop unnecessary columns
# temp[,c("want_angle", "radians", "x_temp", "y_temp")] <- NULL
#
#
# ## join the size data to the 'temp' data
# trees <- left_join(temp, size, by = 'plot_treeid')
# # remove columns w/ no diameter data
# trees <- trees[is.na(trees$dia_calc) == FALSE,]
#
# ## convert to an sf object
# ## first make the trees points
# trees_sf <- st_as_sf(trees, coords = c("x", "y"))
# ## add a buffer with a radius that's the size of the longest leaf
# trees_sf <- st_buffer(x = trees_sf, dist = ((trees_sf$dia_calc/2)/1000))
# # add a column for site
# trees_sf$Site <- "laSelva"
#
# # ## test w/ a plot
# # ggplot(data = trees_sf[trees_sf$plot == "L1",]) +
# #   geom_sf(aes(col = as.factor(dia_year))) +
# #   facet_wrap(~plot) +
# #   theme_classic()
#
# ## get survival information
# ## assing an arbitrary index to each row in 'butterfly'
# trees_sf$index <- 1:nrow(trees_sf)
# ## put an NA for survival
# trees_sf$survs_tplus1_ACTUAL <- NA
# rm("datOut")
# for (i in unique(trees_sf$plot)) {
#   ## get data just for one quad
#   temp <- trees_sf[trees_sf$plot == i,]
#   for (j in unique(temp$tree)) {
#     ## get data just for one individual
#     temp_1 <- temp[temp$tree == j, ]
#     ## sort by year
#     temp_1 <- temp_1[order(temp_1$dia_year),]
#     ## make sure that there is data for more than one year (otherwise just data
#     # for one year, which already has an NA for survival)
#     if (nrow(temp_1) > 1) {
#       ## get a vector of the survival values
#       temp_1$survs_tplus1_ACTUAL <- c(rep.int(1,(nrow(temp_1) - 1)),0)
#     } else if (nrow(temp_1 == 1) & sum(!(temp_1$dia_year != 2018)) == 0) {
#       temp_1$survs_tplus1_ACTUAL <- 0
#     }
#     if (exists("datOut") == FALSE) {
#       datOut <- temp_1
#     } else {
#       datOut <- rbind(datOut, temp_1)
#     }
#   }
# }
#
# for (i in 1:length(plots$plot)) {
#   if (i == 1) {
#     inv <- list(1997:2018)
#     names(inv) <- plots$plot[i]
#   } else {
#     inv[[plots$plot[i]]] <- c(1997:2018)
#   }
# }
#
# ## try running the dataset through trackSpp
# test <- trackSpp(dat = trees_sf, inv = inv, dorm = 0, buff = .1, clonal = FALSE,
#                  species = "genspcode", quad = "plot", year = "dia_year")
#
# ## compare the actual to trackSpp data
#
# ## same number of individual trackIDs as unique plot_treeids?
# length(unique(test$trackID))
# length(unique(test$plot_treeid))
