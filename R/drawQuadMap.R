

drawQuadMap() <- function (dat,
                           species = "Species",
                           site = "Site",
                           quad = "Quad",
                           year = "Year",
                           geometry = "geometry") {

  ## get the boundary box dimensions
  quadBox <- st_bbox(dat)

  ## get the number of years in the dat
  numYears <- length(unique(dat$Year)) + 1
  ## get the number of rows we need
  numRows <- ceiling(numYears/4)

  ## change the margins of the plots
  par(mar = c(1,1,1,1))
  ## change the size of the plotting area
  par(mfrow = c(4, numRows))

  ## plot the quadrat maps
  for (i in sort(unique(dat$Year))) {
    datTemp <- dat[dat$Year == i,]
    plot(datTemp$geometry, col = as.factor(datTemp$Species),
         main = paste(unique(datTemp$Year)),
         xlim = c(quadBox$xmin, quadBox$xmax),
         ylim = c(quadBox$ymin, quadBox$ymax))
    lines(x = c(quadBox$xmin, quadBox$xmin),
          y = c(quadBox$ymin, quadBox$ymax))
    lines(x = c(quadBox$xmax, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymax))
    lines(x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymin, quadBox$ymin))
    lines(x = c(quadBox$xmin, quadBox$xmax),
          y = c(quadBox$ymax, quadBox$ymax))
  }
  ## make a 'plot' that is a legend
  plot(x = NULL, y = NULL,
       #xlim = c(quadBox$xmin, quadBox$xmax),
       #ylim = c(quadBox$ymin, quadBox$ymax)
       )
  legend(x = quadBox$xmin,
         y = quadBox$ymax,
         legend = unique(datTemp$Species),
         fill = as.factor(unique(datTemp$Species)),
         bty = "n",
         cex = 1,
         x.intersp= .2)

  dev.off()


  }
# test --------------------------------------------------------------------

dat <- dat[dat$Site == "CO" &
           dat$location == "unun_11",]
