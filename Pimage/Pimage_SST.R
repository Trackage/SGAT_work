load("D:\\projects\\SGAT\\Movement_fit.Rdata")
fit$model$time <- fit$model$twilight

library(SGAT)
library(raadtools)
x <- model.bin(fit)
sst <- crop(readsst(), extent(x[]), snap = "out")

x <- model.bin(fit, grid = sst)


cl2poly <- function(x) {
    x <- x[[1]]
    x <- cbind(x$x, x$y)
    SpatialPolygons(list(Polygons(list(Polygon(rbind(x, x[1,]))), "1")))
}
i <- 1
##x0 <- SGAT:::as.local.pimg(x[[i]])
##cl <- cl2poly(contourLines(x0, levels = quantile(x0$z[x0$z > 0], 0.5)))
##cn <- extract(sst, cl, cellnumbers = TRUE)

xsst <- x
for (i in seq_along(xsst)) {
    x0 <- SGAT:::as.local.pimg(x[[i]])
    xsst[[i]]$image <- t(as.matrix(resample(sst, raster(x0))))[,ncol(x0$z):1 ]
}




x <- model.bin(fit)
curr <- readcurr(magonly = TRUE)

a <- crop(curr, extent(projectExtent(x[], projection(curr))), snap = "out")



fit$model$time <- fit$model$time - 365 * 24 * 3600
curr <- readcurr(fit$model$time, magonly = TRUE, xylim = extent(a))


x <- model.bin(fit, a)
xcurr <- x
for (i in seq_along(xsst)) {
    x0 <- SGAT:::as.local.pimg(x[[i]])
    xcurr[[i]]$image <- t(as.matrix(resample(curr, raster(x0))))[,ncol(x0$z):1 ]
}

