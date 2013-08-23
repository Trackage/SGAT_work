
load("D:\\projects\\SGAT\\Movement_fit.Rdata")
library(SGAT)

prj <-  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
##pZ <- Pimage(fit$model$twilight, grid = raster(extent(129, 142, -51, -27), 300, 150, crs = prj))
pZ <- Pimage(fit$model$twilight, grid = raster(extent(110, 170, -70, -5), 300, 150, crs = prj))
pZ <- chain.bin(fit$z, pZ)


## x is a Pimage
midProj <- function(x, proj = "laea") {
    ## mid point of longlat data
    mid <- apply(bbox(x[1]), 1, mean)
    sprintf("+proj=%s +lon_0=%.1f +lat_0=%.1f +datum=WGS84", proj, mid[1], mid[2])

}


iproj <- "laea"

## reproject the z's

newproj <- midProj(pZ, iproj)
z <- fit$z
library(rgdal)
for (i in seq_len(dim(z)[3])) {
    z[,,i] <- project(z[,,i], newproj)
}

## use the extents and general resolution of the original
projZ <- Pimage(fit$model$twilight, grid = projectRaster(pZ[1], crs = newproj))
projZ <- chain.bin(z, projZ)

plot(projZ[], axes = FALSE, box = FALSE)








## create graticule, just for fun
gratinput <- as(extent(pZ[]), "SpatialPolygons")
projection(gratinput) <- projection(pZ[])
grat <- gridlines(gratinput)
gratlab <- gridat(gratinput, side = "EN", offset = 0)
projection(grat) <- projection(pZ[1])
projection(gratlab) <- projection(pZ[1])
grat <- spTransform(grat, CRS(newproj))
gratlab <- spTransform(gratlab, CRS(newproj))

## and world
library(rgeos)
library(maptools)
data(wrld_simpl)
wrld_simpl <- gIntersection(wrld_simpl, as(extent(-175, 175, -80, 0), "SpatialPolygons"))
## spTransform(gIntersection(wrld_simpl, as(extent(pZ[1]), "SpatialPolygons")), CRS(newproj))
wrld <- suppressWarnings(
    spTransform(
        wrld_simpl[as.vector(gIntersects(wrld_simpl, as(extent(pZ[1]), "SpatialPolygons"), byid = TRUE)), ],
                CRS(newproj))
    )


plot(wrld, add = TRUE, col = rgb(0.75, 0.75, 0.75, 0.5), usePolypath = FALSE) ## FALSE or no transparency
    plot(grat, add = TRUE, lty = 2)
    text(coordinates(gratlab), labels=parse(text=gratlab$labels), pos=gratlab$pos, offset=gratlab$offset)
title(newproj)


