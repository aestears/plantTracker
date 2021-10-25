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
#                  species = "genspcode", quad = "plot", year = "dia_year",
#                  flagSuspects = TRUE, shrink = .5)
# ## compare the actual to trackSpp data
#
# ## same number of individual trackIDs as unique plot_treeids?
# ## make trackIDs include plot info
# test$uniqueID <- paste0(test$plot,"_",test$trackID)
# ## remove obs. that are flagged as 'suspect'
# testBest <- test
#
# length(unique(testBest$uniqueID))
# length(unique(testBest$plot_treeid))
#
# testBest$survDiff <- testBest$survives_tplus1 - testBest$survs_tplus1_ACTUAL
#
# testTest <- testBest[is.na(testBest$survDiff) == TRUE |
#                    testBest$survDiff == 1 |
#                    testBest$survDiff == (-1)
#                    ,]
#
# testBad <- testBest[!((is.na(testBest$survs_tplus1_ACTUAL)==TRUE & is.na(testBest$survives_tplus1)==TRUE) |
#                       (testBest$survs_tplus1_ACTUAL == 1 & testBest$survives_tplus1 == 1) |
#                       (testBest$survs_tplus1_ACTUAL == 0 & testBest$survs_tplus1_ACTUAL == 0))
#                     ,]
#
# testTest <- testTest[!(testTest$index %in% testBad$index),]
#
# ### values for ms ###
# ## % of surv. assignments that were correct
# nrow(testBest[(is.na(testBest$survs_tplus1_ACTUAL)==TRUE & is.na(testBest$survives_tplus1)==TRUE) |
#            (testBest$survs_tplus1_ACTUAL == 1 & testBest$survives_tplus1 == 1) |
#            (testBest$survs_tplus1_ACTUAL == 0 & testBest$survs_tplus1_ACTUAL == 0)
#          ,])
# (77059)/nrow(testBest)
# # 99.99% (only 3 are wrong)
# ## why are they wrong? : trees are extremely close together, and the assignment chose the wrong 'child'
# mapview(trees_sf[trees_sf$plot_treeid=="P5325",])
# ggplot() +
#   geom_sf(data = st_buffer(trees_sf[trees_sf$plot == "L5" & trees_sf$genspcode == "WARSCOCC" & trees_sf$dia_year == 2012,], .1), col = "red") +
#   geom_sf_text(label = trees_sf[trees_sf$plot == "L5" & trees_sf$genspcode == "WARSCOCC" & trees_sf$dia_year == 2012,]$plot_treeid) +
#   geom_sf(data = st_buffer(trees_sf[trees_sf$plot == "L5" & trees_sf$genspcode == "WARSCOCC" & trees_sf$dia_year == 2013,], .1), col = "orange") +
#   xlim(c(50,60)) +
#   ylim(c(0,10))
#
#   mapview(trees_sf[trees_sf$plot == "L5" & trees_sf$genspcode == "WARSCOCC" & trees_sf$dia_year == 2002,], col.regions = "orange")
#
# ## no. of trackID's assigned
# # actual no.
# length(unique(trees_no_sf$plot_treeid))
# # 5212
# # fxn. no.
# length(unique(paste0(testBest$trackID, "_",testBest$plot)))
# # 5212
#
# ## no. of recruits/quad/year
# # actual no.
# testBest$recruit_actual <- NA
# ## make sure data are in sequential order
# testBest <- testBest[order(testBest$dia_year),]
# ## find the year when a tree was first observed
# ## first, find which are NOT duplicates; also find the
# # observations that were taken in 1997--the first year
# recruitIndex <- testBest[!duplicated(testBest$plot_treeid) &
#                            testBest$dia_year != 1997,]$index
# ## give the recruits a '1' for actual recruit value
# testBest[testBest$index %in% recruitIndex,"recruit_actual"] <- 1
# ## get the plot_treeids of the recruits
# noRecruitIDs <- unique(testBest[testBest$index %in% recruitIndex,]$plot_treeid)
# ##
# testBest[testBest$plot_treeid %in% noRecruitIDs &
#            !(testBest$index %in% recruitIndex),"recruit_actual"] <- 0
# ## remove those obs. from 1997 (give them an 'NA')
# testBest[testBest$dia_year == 1997, "recruit_actual"] <- NA
# badIDs <- testBest[testBest$dia_year == 1997, ]$plot_treeid
# testBest[testBest$plot_treeid %in% badIDs,"recruit_actual"] <- NA
#
# ## compare assigned 'recruit' to 'recruit_true'
# testBest$recDiff <- testBest$recruit_actual - testBest$recruit
# testBest[testBest$recruit == 0 & testBest$recruit_actual == 0 &
#            is.na(testBest$recruit) == FALSE &
#            is.na(testBest$recruit_actual) == FALSE, "recDiff"] <- 0
# #16120
# testBest[testBest$recruit == 1 & testBest$recruit_actual == 1 &
#            is.na(testBest$recruit) == FALSE &
#            is.na(testBest$recruit_actual) == FALSE, "recDiff"] <- 0
# #1663
# testBest[is.na(testBest$recruit) &
#            is.na(testBest$recruit_actual), "recDiff"] <- 0
# #59254
#
# testBest[is.na(testBest$recDiff), "recDiff"] <- 9999
#
#
# # fxn. no.
# test_recruits <- getRecruits(dat = testBest, byGenet = TRUE, species = "genspcode", quad =
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
# names(butterfly)[9] <- "survives_tplus1_actual"
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
# testGood <- test[test$Suspect == FALSE,]
# testGood$surv_diff <- testGood$survives_tplus1 - testGood$survives_tplus1_actual
# testGood <- testGood[testGood$surv_diff == 0 |
#                         (is.na(testGood$survives_tplus1) == TRUE &
#                            is.na(testGood$survives_tplus1_actual) == TRUE) ,]
#
# testBad <- test[!(test$index %in% testGood$index), ]
#
#
# ## view some 'bad' obs
# ggplotly(ggplot() +
#   #geom_sf(dat = test[test$Plot_ID == "D10" & test$Year == 2018,], alpha = .5) +
#   geom_sf(dat = st_buffer(test[test$Plot_ID == "D7" & test$trackID == "OENCOL_2019_70",], 2),colour = "red", alpha = .5) +
#   geom_sf(dat = test[test$Plot_ID == "D7" & test$ID == 220,],
#           aes(fill = Year), colour = NA, alpha = .9) +
#   theme_classic())
#
# ggplotly(ggplot(dat = test[test$Plot_ID == "C5" & (test$ID %in%  c(74, 174) | test$trackID == "OENCOL_2018_18"),]) +
#            geom_sf(aes(fill = Year), colour = NA, alpha = .9) +
#            geom_sf_text(aes(label = trackID)) +
#            geom_sf_text(aes(label = ID), nudge_y = .5) +
#            theme_classic())
#
# ## remove 'suspect' observations
# ## remove observations from the last year, since they aren't helpful here
# testBest <- test[test$Suspect == FALSE,]
#
# ### values for ms ###
# ## % of surv. assignments that were correct
# nrow(testBest[is.na(testBest$survives_tplus1_actual) == TRUE &
#        is.na(testBest$survives_tplus1) == TRUE,]) # 1634 got NA for both real
# # and fake surv.
# sum(testBest$surv_diff == 0, na.rm = TRUE) # 3591
# table(testBest$surv_diff)
# (1634 + 3591) / nrow(testBest)
# # 99.6%
#
# ## for some, assigns a small plant the same trackID, even when it's likely a new
# # seedling that was close to the parent plant
# ## others, for indivdiuals that are super close together, the function picks the wrong one for the next year
#
# ## no. of trackID's assigned
# # actual no.
# length(unique(paste0(testBest$ID, "_",testBest$Plot_ID)))
# # 3146
# # fxn. no.
# length(unique(paste0(testBest$trackID, "_",testBest$Plot_ID)))
# # 3128
#  # (99.4% accuracy)
#
# ## no. of recruits/quad/year (use 'test' d.f, which contains the 2020 data)
# # actual no.
# ## no. of recruits/quad/year
# # actual no.
# testBest$recruit_actual <- NA
# testBest$medName <- paste0(testBest$Plot_ID,"_",testBest$ID)
# ## make sure data are in sequential order
# testBest <- testBest[order(testBest$Year),]
# ## find the year when a plant was first observed
# ## first, find which are NOT duplicates; also find the
# # observations that were taken in 2018--the first year
# recruitIndex <- testBest[!duplicated(testBest$medName) &
#                            testBest$Year != 2018,]$index
# ## give the recruits a '1' for actual recruit value
# testBest[testBest$index %in% recruitIndex,"recruit_actual"] <- 1
# ## get the plot_treeids of the recruits
# noRecruitIDs <- unique(testBest[testBest$index %in% recruitIndex,]$medName)
# ##
# testBest[testBest$medName %in% noRecruitIDs &
#            !(testBest$index %in% recruitIndex),"recruit_actual"] <- 0
# ## remove those obs. from 2018 (give them an 'NA')
# testBest[testBest$Year == 2018, "recruit_actual"] <- NA
# badIDs <- testBest[testBest$Year == 2018, ]$medName
# testBest[testBest$medName %in% badIDs,"recruit_actual"] <- NA
#
# ## compare assigned 'recruit' to 'recruit_true'
# testBest$recDiff <- NA
# testBest[testBest$recruit == 0 & testBest$recruit_actual == 0 &
#            is.na(testBest$recruit) == FALSE &
#            is.na(testBest$recruit_actual) == FALSE, "recDiff"] <- 0
# #687
# testBest[testBest$recruit == 1 & testBest$recruit_actual == 1 &
#            is.na(testBest$recruit) == FALSE &
#            is.na(testBest$recruit_actual) == FALSE, "recDiff"] <- 0
# #1521
# testBest[is.na(testBest$recruit) &
#            is.na(testBest$recruit_actual), "recDiff"] <- 0
# #3013
#
# testBest[is.na(testBest$recDiff), "recDiff"] <- 9999
#
