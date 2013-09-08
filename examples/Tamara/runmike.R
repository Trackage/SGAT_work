
library(tripEstimation)
library(raster)
library(maptools)

setwd("E:/DATA/people/TamaraB/4782_light")
load("4782a_t1.1.Rdata")
load("4782a_calib.Rdata")
## recent versions of R need this fix for saved calibrations
##plot(calib(-10:10))
####Error in calib(-10:10) : object 'C_R_approxfun' not found
##FIX:
 body(calib) <- parse(text = ".approxfun(x, y, v, method, yleft, yright, f)")

lon.min <- -180
lon.max <- -110
lat.min <- 0
lat.max <- 65

start <- c(-122.338015, 37.109136)
end <- start

fixed.release <- TRUE
fixed.recapture <- TRUE

## how many twilights plus fixed points?
m <- max(d1$segment) + fixed.release + fixed.recapture


## "topo" is a basic "xyz" list as used by ?image, with a topo$z matrix that
## encodes spatially TRUE/FALSE for land/sea (plus a bit of jigger)
data(wrld_simpl)
mask <- rasterize(wrld_simpl, raster(xmn = lon.min, xmx = lon.max, ymn = lat.min, ymx = lat.max, nrows = 120, ncols = 80))

mask <- is.na(mask)
mask <- as.image.SpatialGridDataFrame(as(mask, "SpatialGridDataFrame"))

lookup <- mkLookup(mask, FALSE)


proposal.x <- norm.proposal(m, 3, c(0.5, 0.5, 0.05))
proposal.z <- norm.proposal(m-1, 2, c(1, 1))

speed.mu <- 2
speed.sigma <- 2

elev <- elevation(start[1], start[2], solar(d1$gmt))

d.calib <- d1[d1$segment %in% 1:2, ]
elev <- elev[d1$segment %in% 1:2]
b <- approxfun(elev, d.calib$light, rule = 2)

## create the model
d.model <- solar.model(d1$segment,d1$gmt,d1$light,
                       proposal.x = proposal.x$proposal,
                       proposal.z = proposal.z$proposal,
                       fix.release = fixed.release,
                       fix.recapture = fixed.recapture,
                       mask.x = lookup,
                       mask.z = lookup,
                       calibration=b,
                       behav = "speed",
                       behav.dist = "log",
                       light.sigma=3,
                       k.sigma = 4,
                       behav.mean = speed.mu,
                       behav.sd = speed.sigma,
                       ekstrom = c(-8, 5, 128))

## find some starting points that satisfy the land mask
xx <- seq(lon.min, lon.max, length = 30)
yy <- seq(lat.min, lat.max, length = 30)
plg <- position.logp(d.model, xx, yy, xrest = 0,
start = start, end = end, prob = 0.9)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# initialization

X <- plg$X

xfile <- "initx.bin"
zfile <- "initz.bin"

ch <- metropolis0(d.model,iters=10,thin=1,start.x = X,
                  start.z = (X[-m,1:2]+X[-1,1:2])/2)

chain.write(xfile, ch$x)
chain.write(zfile, ch$z)

while(!all(ch$mask.x) | !all(ch$mask.z)) {
    ch <- metropolis0(d.model, iters=20,
    				  thin=20, start.x=ch$last.x,
				      start.z=ch$last.z)
    chain.write(xfile, ch$x, append = TRUE)
    chain.write(zfile, ch$z, append = TRUE)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Estimation


xfile <- "burnx.bin"
zfile <- "burnz.bin"

ch <- metropolis(d.model,iters=100,thin=20,start.x = ch$last.x,
                  start.z = ch$last.z)

chain.write(xfile, ch$x)
chain.write(zfile, ch$z)


## examine burn-in (do this after a few runs, before committing
## to the X.bin/Z.bin below
opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
k <- sample(m,20)
matplot(t(ch$x[k,1,]),type="l",lty=1,col="dodgerblue",ylab="Lon")
matplot(t(ch$x[k,2,]),type="l",lty=1,col="firebrick",ylab="Lat")
par(opar)


## repeat this loop a few times after checking a couple
 for (k in 1:5) {

for (i in 1:3) {

    ch <- metropolis(d.model,iters=200,thin=10,start.x = ch$last.x,
                      start.z = ch$last.z)

    chain.write(xfile, ch$x, append = TRUE)
    chain.write(zfile, ch$z, append = TRUE)

}

## examine burn-in (do this after a few runs, before committing
## to the X.bin/Z.bin below
opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
k <- sample(m,20)
matplot(t(ch$x[k,1,]),type="l",lty=1,col="dodgerblue",ylab="Lon")
matplot(t(ch$x[k,2,]),type="l",lty=1,col="firebrick",ylab="Lat")
par(opar)

}

for (i in 1:3) {

    ch <- metropolis(d.model,iters=200,thin=10,start.x = ch$last.x,
                      start.z = ch$last.z)

    chain.write(xfile, ch$x, append = TRUE)
    chain.write(zfile, ch$z, append = TRUE)
    proposal.x$tune(ch$x, 0.3)
    proposal.z$tune(ch$z, 0.3)
}


## START SAVING
xfile <- "X.bin"
zfile <- "Z.bin"
for (i in 1:20) {
    ch <- metropolis(d.model,iters=2000,thin=10,start.x = ch$last.x,
                     start.z = ch$last.z)
    xm <- apply(ch$z, 1:2, mean)
    plot(xm)
    plot(wrld_simpl, add = TRUE)
    points(xm)
    proposal.x$tune(ch$x, 0.3)
    proposal.z$tune(ch$z, 0.3)
    chain.write(xfile, ch$x, append = i > 1)
    chain.write(zfile, ch$z, append = i > 1)
}

