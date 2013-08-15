library(SGAT)

if ( setInternet2(NA)) {
    source("C:\\Users\\michae_sum\\Documents\\GitHub\\SGAT\\R\\Curve.R")
    load("D:\\data\\TAGDATA\\mcmcGeo\\b362_99\\readyRun.Rdata")
    load("D:\\data\\TAGDATA\\mcmcGeo\\b362_99\\old\\metropolisMarch06\\x.Rdata")
} else {
    source("C:\\Users\\mdsumner\\Documents\\GitHub\\SGAT\\R\\Curve.R")
    load("E:\\backup\\Elements\\Azimuth\\DATA\\mcmcGeo\\b362_99\\readyRun.Rdata")
load("E:\\backup\\Elements\\Azimuth\\DATA\\mcmcGeo\\b362_99\\old\\metropolisMarch06\\x.Rdata")
}


## don't ask
body(calib) <- parse(text = ".approxfun(x, y, v, method, yleft, yright, f)")

d1 <- d[!is.na(d$segment), ]
d1 <- d1[!is.na(d1$light), ]
d1$segment <- unclass(factor(d1$segment))

library(raster)
library(maptools)
data(wrld_simpl)
land.mask <- function(xlim,ylim,n=4,land=TRUE) {
  r <- raster(nrows=n*diff(ylim),ncols=n*diff(xlim),
              xmn=xlim[1],xmx=xlim[2],
              ymn=ylim[1],ymx=ylim[2],
              crs=proj4string(wrld_simpl))
  r <- rasterize(wrld_simpl,r)
  r <- as.matrix(is.na(r))[nrow(r):1,]
  if(land) r <- !r
  xbin <- seq(xlim[1],xlim[2],length=ncol(r)+1)
  ybin <- seq(ylim[1],ylim[2],length=nrow(r)+1)

  function(p) {
    r[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
  }
}

is.sea <- land.mask(xlim = c(lon.min, lon.max),
                    ylim= c(lat.min, lat.max),
                    n = 4,land = FALSE)

log.prior <- function(p)  {
  f <- is.sea(p)
  ifelse(f|is.na(f),0,-1000)
}



x0 <- x[,1:2]

fixedx <- seq_len(nrow(x0)) %in% c(1, nrow(x0))

z0 <-  (x0[-nrow(x0),1:2] + x0[-1,1:2])/2
model <- curve.model(d1$gmt, d1$light, d1$segment, calib,
                     twilight.model = "LogNormal",
                     alpha = c(7, 10), beta = c(150, 80),
                     logp.x = log.prior,
                     logp.z = log.prior,
                     x0 = x0,z0 = z0,fixedx = fixedx)

x.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nrow(x0))
z.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nrow(x0)-1)

fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=100,thin=20)

## testthat
## fixedx[i/n] is really a logical vector of length n, no NAs (and friends)
## dt is all valid, no funky
## all arguments are right length, sign and mode
## calib returns valid values that are increasing over -10:10


all(fixedx | is.sea(chain.last(fit$x)))
all(is.sea(chain.last(fit$z)))

