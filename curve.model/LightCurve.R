### R code from vignette source 'C:/Users/michae_sum/Documents/GIT/SGAT/vignettes/LightCurve.Rnw'

###################################################
### code chunk number 1: first
###################################################
set.seed(42)
library(SGAT)
data(ElephantSeal2)
data(ElephantSeal2calib)

calibration <- with(ElephantSeal2calib, approxfun(zenith, light, rule = 2))

library(raster)
library(maptools)
data(wrld_simpl)

poly <- wrld_simpl

lonlim <- c(110, 190)
latlim <- c(-70, -50)





land_dist.mask <- function(xlim, ylim, n=4) {
    data(wrld_simpl)
    r <- raster(nrows = n * diff(ylim), ncols = n * diff(xlim),
              xmn = xlim[1L], xmx = xlim[2L],
              ymn= ylim[1L], ymx = ylim[2L],
              crs = projection(wrld_simpl))

    r <- rasterize(wrld_simpl, r, field = "FIPS")

    ## clip as lines, so we don't get the sealed edges from the intersection
    outline <- gIntersection(as(wrld_simpl, "SpatialLines"), as(extent(r), "SpatialPolygons"))
    wxy <- coordinates(as(outline, "SpatialPoints"))

    ## NA is ocean, the rest need distance values
    ## prepare for distance to coast
    rxy <- as.data.frame(r, xy = TRUE)

    for (i in seq_len(nrow(rxy))[!is.na(rxy$layer)]) {
         dsts <- spDistsN1(wxy, as.matrix(rxy[i,c("x", "y"), drop = FALSE]), longlat = TRUE)
         rxy$layer[i] <- min(dsts)

     }
    r <- setValues(r, rxy$layer)
    r <- as.matrix(r)[nrow(r):1, ]

    xbin <- seq(xlim[1],xlim[2],length=ncol(r)+1)
    ybin <- seq(ylim[1],ylim[2],length=nrow(r)+1)
    function(p) {
        r[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
    }
}

dist.sea <- land_dist.mask(xlim=lonlim,ylim=latlim,n=4)

log.prior <- function(x) {
    e <- dist.sea(x[,1:2])
    ifelse(is.na(e), 0, -e)
}



###################################################
### code chunk number 6: model
###################################################
model <- curve.model(ElephantSeal2$time, ElephantSeal2$light, ElephantSeal2$segment, calibration,
                     logp.x = log.prior,
                     logp.z = log.prior,
                     alpha = c(7, 10), beta = c(8, 3.5)/10)


###################################################
### code chunk number 7: init
###################################################
## find starting points
nx <- 30L
ny <- 25L
grid <- list(x = seq(lonlim[1L], lonlim[2L], length = nx),
             y = seq(latlim[1L], latlim[2L], length = ny),
             z = array(0, c(nx, ny, length(model$time))))
for (i in seq_along(grid$x)) {
        for (j in seq_along(grid$y)) {
            grid$z[i,j,] <- model$logpx(cbind(rep(grid$x[i], length(model$time)), grid$y[j], 0))
        }
    }

x0 <- cbind(as.matrix(expand.grid(grid$x, grid$y)[apply(grid$z, 3, which.max), ]), 0)

## the first and last locations are known
x0[1L, c(1L, 2L) ] <- c(158.950, -54.5)
x0[nrow(x0), c(1L, 2L)] <- x0[1L, c(1L, 2L)]

fixedx <- seq_len(nrow(x0)) %in% c(1, nrow(x0))



###################################################
### code chunk number 8: proposals
###################################################
x.proposal <- mvnorm(S=diag(c(0.005,0.005, 0.05)),n=nrow(x0))
z.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nrow(x0)-1)


###################################################
### code chunk number 9: run1
###################################################
fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=20,thin=20, x0 = x0, z0 = (x0[-nrow(x0),1:2]+ x0[-1,1:2])/2)

plot(x0)
##segments(x0[,1], x0[,2], chain.last(fit$x)[,1], chain.last(fit$x)[,2])

lastx <- chain.last(fit$x)
while((!all(fixedx | log.prior(chain.last(fit$x)) == 0)) | !all(log.prior(chain.last(fit$z)) == 0)) {
fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=40,thin=20, x0 = chain.last(fit$x), z0 = chain.last(fit$z))
##plot(rdist)
##points(lastx)
segments(lastx[,1], lastx[,2], chain.last(fit$x)[,1], chain.last(fit$x)[,2])
lastx <- chain.last(fit$x)


}



x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)




opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
k <- sample(nrow(x0),20)
matplot(t(cbind(fit0$x[k,1,], fit$x[k,1,])),type="l",lty=1,col="dodgerblue",ylab="Lon")
matplot(t(cbind(fit0$x[k,1,], fit$x[k,1,])),type="l",lty=1,col="firebrick",ylab="Lat")
par(opar)

x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)

## Tune proposals based on previous run
x.proposal <- mvnorm(chain.cov(fit$x,discard=100),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z,discard=100),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)
## Tune proposals based on previous run
x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)

opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
matplot(scale(t(fit$x[150,,]),scale=F),type="l",lty=1,col=c(2,4),
        xlab="",ylab=expression(x[150]))
matplot(scale(t(fit$z[150,,]),scale=F),type="l",lty=1,col=c(2,4),
        xlab="",ylab=expression(z[150]))
par(opar)

## Draw final sample
x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=2000,thin=50)






library(raadtools)
d0 <- readsst(xylim = extent())
p <- model.bin(fit, grid = d0)

ps <- p
mkxy <- function(x) expand.grid(x = x$x, y = x$y)
for (i in seq_along(p)) {
    cn <- cellFromXY(d0, mkxy(SGAT:::as.local.pimg(ps[[i]])))
    ok <- ps[[i]]$image > 0
    ps[[i]]$image[ok] <- extract(d0, cn)[ok]
}



##' Extract gridded data from a \code{\link{Pimage}}
##'
##' Read data via function in raadtools queried with a Pimage object.
##' @title extract.pimg
##' @param x a raadtools read function
##' @param y a Pimage object, must match an object returned by x
##' @param ... arguments passed to function x
##' @return Pimage
extract.pimg <- function(x, y, ...) {
    times <- attr(y, "times")

    ## assume daily
    ct <- cut(times, "1 day")

    for (i in seq_len(nlevels(ct)) {
        thistime <- times[ct == levels(ct)[i]][1]
        d <- x(thistime, ...)
    }

}
