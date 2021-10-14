# library(tidyverse)
# library(sf)
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
# ## test w/ a plot
# ggplot(data = trees_sf[trees_sf$plot == "L1" & trees_sf$dia_year %in% c(1997:2002),]) +
#   geom_sf(aes(col = as.factor(dia_year))) +
#   facet_wrap(~dia_year) +
#   theme_classic() +
#   scale_colour_discrete(guide = "none")
#
# ## get survival information
# ## assign an arbitrary index to each row in 'trees'
# trees_sf$index <- 1:nrow(trees_sf)
# trees_no_sf <- st_drop_geometry(trees_sf)
# ## put an NA for survival
# trees_no_sf$survs_tplus1_ACTUAL <- NA
# rm("datOut")
# for (i in unique(trees_no_sf$plot)) {
#   ## get data just for one quad
#   temp <- trees_no_sf[trees_no_sf$plot == i,]
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
# ## fix survival for 2018 (is an 'NA' instead of a 0)
# datOut[datOut$dia_year == 2018, 'survs_tplus1_ACTUAL'] <- NA
# ## put the geometry data back into the 'trees_no_sf' data.frame
# trees_sf <- left_join(x = trees_sf, y = datOut)
#
# ## get the correct inventory list
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
# test <- trackSpp(dat = trees_sf, inv = inv, dorm = 0, buff = .01, clonal = FALSE,
#                  species = "genspcode", quad =
#                    "plot", year = "dia_year", flagSuspects = TRUE)
#
# ## compare the actual to trackSpp data
#
# ## same number of individual trackIDs as unique plot_treeids?
# ## make trackIDs include plot info
# test$uniqueID <- paste0(test$plot,"_",test$trackID)
# length(unique(test$uniqueID))
# length(unique(test$plot_treeid))
#
# test$survDiff <- test$survives_tplus1 - test$survs_tplus1_ACTUAL
#
# testTest <- test[is.na(test$survDiff) == TRUE |
#                    test$survDiff == 1 |
#                    test$survDiff == (-1)
#                    ,]
#
# testBad <- testTest[is.na(testTest$survs_tplus1_ACTUAL) == TRUE &
#                        is.na(testTest$survives_tplus1) == TRUE &
#                       testTest$dia_year != 2018
#                        ,]
#
# testTest <- testTest[!(testTest$index %in% testBad$index),]
#
# ### values for ms ###
# ## % of surv. assignments that were correct
# (77062 - 6)/77062
# # 99.99%
#
# ## no. of trackID's assigned
# # actual no.
# length(unique(trees_no_sf$plot_treeid))
# # 5212
# # fxn. no.
# length(unique(paste0(test$trackID, "_",test$plot)))
# # 5212
#
# ## no. of recruits/quad/year
# # actual no.
# for (i in 1:length(unique(trees$plot))) {
#   ## get data for just one plot
#   temp <- trees[trees$plot == unique(trees$plot)[i],]
#   ## make sure it is in sequential order
#   temp <- temp[order(temp$dia_year),]
#
#   ## find the year when a tree was first observed
#   ## first, find which are duplicates, then remove them
#   temp <- temp[!duplicated(temp$plot_treeid),]
#   temp$recruit_true <- 1
#
#   if (i == 1) {
#     recruitOut <- temp
#   } else {
#     recruitOut <- rbind(recruitOut, temp)
#   }
# }
# # remove values for 1997, since we can't know if it was really a recruit then
# recruitOut <- recruitOut[recruitOut$dia_year != 1997,]
# recruitOut <- aggregate(x = recruitOut$recruit_true, by =
#                           list("plot" = recruitOut$plot,
#                                "genspcode" = recruitOut$genspcode,
#                                "dia_year" = recruitOut$dia_year),
#                         FUN = sum)
#
# # fxn. no.
# test_recruits <- getRecruits(dat = test, byGenet = TRUE, species = "genspcode", quad =
#               "plot", year = "dia_year")
# # compare
# testRecs <- merge(recruitOut, test_recruits, by = c("plot", "dia_year", "genspcode"))
# testRecs$diff <- testRecs$x - testRecs$recruits
# unique(testRecs$diff)
# # no differences!
#
# # #### COBP test ####
# #
# # ## testing using COBP data
# load("/Users/Alice/Dropbox/Grad School/Research/Oenothera coloradensis project/Processed_Data/spatial_COBP.RData")
#
# ## add necessary columns
# butterfly$Species <- "Oenothera coloradensis"
# names(butterfly)[10] <- "survives_tplus1_actual"
# ## make a quadrat inventory
# cobpInv <- list("C4" = c(2018:2020),
#                 "C5" = c(2018:2020),
#                 "C8" = c(2018:2020),
#                 "D10" = c(2018:2020),
#                 "D11" = c(2018:2020),
#                 "D7" = c(2018:2020),
#                 "S1" = c(2018:2020),
#                 "S2" = c(2018:2020),
#                 "S3" = c(2018:2020),
#                 "S4" = c(2018:2020),
#                 "S5" = c(2018:2020),
#                 "S6" = c(2018:2020),
#                 "S7" = c(2018:2020),
#                 "S8" = c(2018:2020),
#                 "S9" = c(2018:2020),
#                 "U3" = c(2018:2020),
#                 "U4" = c(2018:2020),
#                 "U6" = c(2018:2020)
#               )
#
# test <- trackSpp(dat = butterfly, inv = cobpInv, dorm = 0, buff = 1, clonal = FALSE, site = "Location", quad = "Plot_ID", aggByGenet = FALSE, flagSuspects = TRUE, shrink = .3)
#
# mapview(butterfly[butterfly$Plot_ID == "D10" & butterfly$Year == 2018,], col.regions = "red") + mapview(butterfly[butterfly$Plot_ID == "D10" & butterfly$Year == 2019,], col.regions = "orange") + mapview(butterfly[butterfly$Plot_ID == "D10" & butterfly$Year == 2020,], col.regions = "yellow") + mapview(butterfly[butterfly$Plot_ID == "D10" & butterfly$Year == 2018 & butterfly$ID == 1,])
#
# ggplot(dat = test[test$Plot_ID == "D10" &
#                     test$trackID %in% c("OENCOL_2018_1","OENCOL_2018_2","OENCOL_2018_3","OENCOL_2018_4","OENCOL_2018_5" ),]) +
#   geom_sf(aes(color = Year, fill = trackID), alpha = .5)
#
# plot(butterfly[butterfly$Plot_ID == "D10" & butterfly$Year == 2018,]$geometry,
#      col = "red")
#
# test$surv_diff <- test$survives_tplus1 - test$survives_tplus1_actual
# test$bigNames <- paste(test$Plot_ID, test$trackID, test$Year, sep = "_")
#
# ## find 'bad' assignments
# # first, remove obs. from the last year
# testGood <- test[test$Year != 2020,]
# testGood$surv_diff <- testGood$survives_tplus1 - testGood$survives_tplus1_actual
# testGood_1 <- testGood[testGood$surv_diff == 0,]
# testGood_2 <- testGood[is.na(testGood$survives_tplus1) == TRUE &
#                          is.na(testGood$survives_tplus1_actual) == TRUE,]
# testGood <- rbind(testGood_1, testGood_2)
# testGood[testGood$surv_diff==-1 & is.na(testGood$surv_diff) == FALSE,]
# testGood[testGood$Suspect == TRUE,]
#
# testBad <- test[test$Year != 2020 &
#                   !(test$bigNames %in% testGood$bigNames), ]
#
#
# ## view some 'bad' obs
# ggplot() +
#   #geom_sf(dat = test[test$Plot_ID == "D10" & test$Year == 2018,], alpha = .5) +
#   geom_sf(dat = test[test$Plot_ID == "D10" & test$trackID == "OENCOL_2018_40",], aes(fill = Year), alpha = .5) +
#   geom_sf(dat = test[test$Plot_ID == "D10" & test$ID == 112,],
#           aes(color = Year), fill = NA, alpha = .1, lty = 2)
#
#
# ## remove 'suspect' observations
# ## remove observations from the last year, since they aren't helpful here
# testBest <- test[test$Year != 2020,]
#
# ### values for ms ###
# ## % of surv. assignments that were correct
# nrow(testBest[is.na(testBest$survives_tplus1_actual) == TRUE &
#        is.na(testBest$survives_tplus1) == TRUE,]) # 0 got NA for both real
# # and fake surv.
# sum(testBest$surv_diff == 0, na.rm = TRUE) # 2705
# table(testBest$surv_diff)
# (3587) / nrow(testBest)
# # 99.3%
#
# ## no. of trackID's assigned
# # actual no.
# length(unique(paste0(testBest$ID, "_",testBest$Plot_ID)))
# # 2650
# # fxn. no.
# length(unique(paste0(testBest$trackID, "_",testBest$Plot_ID)))
# # 2641
#  # (99.66% accuracy)
#
# ## no. of recruits/quad/year (use 'test' d.f, which contains the 2020 data)
# # actual no.
# test$names <- paste(test$Plot_ID,test$ID, sep = "_")
# for (i in 1:length(unique(test$Plot_ID))) {
#   ## get data for just one plot
#   temp <- test[test$Plot_ID == unique(test$Plot_ID)[i],]
#   ## make sure it is in sequential order
#   temp <- temp[order(temp$Year),]
#
#   ## find the year when a tree was first observed
#   ## first, find which are duplicates, then remove them
#   temp <- temp[!duplicated(temp$names),]
#   temp$recruit_true <- 1
#
#   if (i == 1) {
#     recruitOut <- temp
#   } else {
#     recruitOut <- rbind(recruitOut, temp)
#   }
# }
# # remove values for 2018, since we can't know if it was really a recruit then
# recruitOut <- recruitOut[recruitOut$Year != 2018,]
# recruitOut <- aggregate(x = recruitOut$recruit_true, by =
#                           list("Plot_ID" = recruitOut$Plot_ID ,
#                                "Species" = recruitOut$Species ,
#                                "Year" = recruitOut$Year),
#                         FUN = sum)
#
# # fxn. no.
# test_recruits <- getRecruits(dat = test, byGenet = TRUE, quad = "Plot_ID", site = "Location")
# # compare
# testRecs <- merge(recruitOut, test_recruits, by = c("Plot_ID", "Year", "Species"))
# testRecs$diff <- testRecs$x - testRecs$recruits
# unique(testRecs$diff)
# testRecs$percRight <- testRecs$recruits / testRecs$x
# mean(testRecs$percRight)
# # 99.8% correct
