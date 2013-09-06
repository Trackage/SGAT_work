
load("D:\\projects\\SGAT\\Movement_fit.Rdata")
library(SGAT)

## Pimage needs no input but the times
pZ <- Pimage(fit$model$twilight)


## can take a grid (anything that raster() accepts or a raster)
prj <-  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
##pZ <- Pimage(fit$model$twilight, grid = raster(extent(129, 142, -51, -27), 300, 150, crs = prj))
pZ <- Pimage(fit$model$twilight, grid = raster(extent(110, 170, -70, -5), 300, 150, crs = prj))


## 0) nothing supplied
pZ <- chain.bin(fit)

## 1) pimg is supplied, and not empty
pZ <- chain.bin(fit, pimg = pZ)
## pimg is supplied, empty
pZ <- chain.bin(fit, pimg = Pimage(fit$model$twilight, grid = raster(extent(110, 170, -70, -5), 300, 150, crs = prj)))

## 2) no pimg, grid supplied
pZ <- chain.bin(fit, grid = list(x = 110:170, y = -70:5, z = matrix(0, 61, 76)))
## raster, default
pZ <- chain.bin(fit, grid = raster())
## raster, custom
pZ <- chain.bin(fit, grid = raster(xmn = 100, xmx = 180, ymn = -70, ymx = -10, nrows = 100, ncols = 120))


mid <- c(140, -40)
proj <- "laea"
newproj <-  sprintf("+proj=%s +lon_0=%.1f +lat_0=%.1f +datum=WGS84", proj, mid[1], mid[2])


## grid supplied, and is projected
pZ <- chain.bin(fit, grid = projectExtent(raster(xmn = 100, xmx = 180, ymn = -70, ymx = -10, nrows = 100, ncols = 120), newproj))

## only proj is supplied
pZ <- chain.bin(fit, proj = newproj)
