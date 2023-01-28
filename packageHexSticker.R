#////////////////////////
# making a Hex Sticker for plantTracker
# Alice Stears
# 27 January 2023
#///////////////////////

# load packages -----------------------------------------------------------
library(hexSticker)
library(plantTracker)
library(tidyverse)
library(svglite)
library(png)
library(sf)

# make map of plot to use in sticker --------------------------------------
# make a ggplot using plantTracker::grasslandData data
oldMai <- par()$mai

# reset the margin
par(mai = c(0,0,0,0))

exampleDat <- grasslandData[grasslandData$Site == "AZ" &
                              grasslandData$Quad == "SG2" ,]
exampleSmall <- exampleDat[exampleDat$Species == "Heteropogon contortus" &
                             exampleDat$Year %in% c(1922:1923),]
exampleSmall$Site <- "AZs"
exampleSmall$Quad <- "SG2"

## trim the dataset to be smaller
## make a bounding box
pl = list(rbind(c(0,0), c(.5,0), c(.5,.7), c(0,.7), c(0,0)))
box <- st_polygon(pl)
exampleSmall <-suppressWarnings(st_intersection(exampleSmall, box))

## get trackIDs
exampleOut <- suppressMessages(plantTracker::trackSpp(dat = exampleSmall, inv = list("SG2" = c(1922:1925)), buff = .05, clonal = TRUE, dorm = 1, buffGenet = .01, aggByGenet = TRUE, printMessages = FALSE))

labels <- data.frame(trackID = unique(exampleOut$trackID),
                     trackID_new = c(1:length(unique(exampleOut$trackID))))
exampleOut$trackID_new <- labels$trackID_new[match( exampleOut$trackID, labels$trackID)]

plotMap <- ggplot(data = exampleOut) +
  geom_sf(aes(color = trackID, fill = trackID), alpha = .9) +
  xlim(c(0,.35)) +
  ylim(c(0.4,.65)) +
  scale_fill_discrete(guide = "none") +
  scale_color_discrete(guide = "none") +
  #labs(title = Year) +
  facet_wrap(~ Year, ncol = 2) +
  theme_void() +
  theme_transparent() +
  theme(panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
    axis.line = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        strip.text = element_blank())

# save to file
svglite(filename = "../plotForPackageHexSticker.svg")
plotMap
dev.off()

## add arrows in InkScape
## read in InkScape file
imgurl <- "../plotForPackageHexSticker_arrows.png"

# make the hex sticker ----------------------------------------------------
s <- sticker(imgurl, package="plantTracker", p_size=20, p_color = "#274e13",
             s_x=1, s_y=.85, s_width=.8,
             h_fill = "#b6d7a8",
             h_color = "#38761d",
             #3e6073,
             filename="inst/figures/imgfile.png")
plot(s)
