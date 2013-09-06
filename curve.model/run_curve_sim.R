library(SGAT)
calib.data <- read.csv("C:\\Users\\michae_sum\\Documents\\GitHub\\SGAT_work\\curve.model\\calib.data.csv")

calibfun <- approxfun(calib.data$zenith, calib.data$light, rule = 2)


## control points c1 and c2.
bezier <- function(t1,p1,c1,t2,p2,c2) {
  function(t) {
    t <- (unclass(t)-unclass(t1))/(unclass(t2)-unclass(t1))
    outer((1-t)^3,p1)+
      outer(3*(1-t)^2*t,c1)+
        outer(3*(1-t)*t^2,c2)+
          outer(t^3,p2)
  }
}
## Endpoints
p0 <- c(130,-30)
p1 <- c(140,-50)
p2 <- c(130,-40)
## Control points
c0 <- p0+c(6,-4)
c1 <- p1+c(-4,8)
c2 <- p1+c(-6,2)
c3 <- p2+c(2,-8)
f1 <- bezier(30,p0,c0,60,p1,c1)
f2 <- bezier(90,p1,c2,120,p2,c3)
days <- c(1:30,31:60,61:90,91:120,121:150)
ps <- rbind(cbind(rep(p0[1],30),rep(p0[2],30)),
            f1(31:60),
            cbind(rep(p1[1],30),rep(p1[2],30)),
            f2(91:120),
            cbind(rep(p2[1],30),rep(p2[2],30)))



## Resample to two minute intervals and compute zenith angles
tms <- as.POSIXct("2013-04-01","GMT")+(24*60*60)*days
tms.out <- seq(min(tms),max(tms),by=2*60)
zen <- zenith.simulate(tms,ps[,1],ps[,2],tms.out)

d <- data.frame(gmt = tms.out, light = calibfun(zen$zenith))
bad <- is.na(d$light)
seg <- c(0, cumsum(abs(diff(bad))))
##light <- split(d$light[!bad], seg[!bad])
##time <- split(d$gmt[!bad], seg[!bad])

d1 <- data.frame(gmt = d$gmt[!bad],
                 light = d$light[!bad],
                 segment = unclass(factor(seg[!bad])))



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

is.sea <- land.mask(xlim = c(125, 145),
                    ylim= c(-55, -25),
                    n = 4,land = FALSE)

log.prior <- function(p)  {
  f <- is.sea(p)
  ifelse(f|is.na(f),0,-1000)
}



x0 <- cbind(ps[seq(1, nrow(ps), length = max(d1$segment)), ], 0)


fixedx <- seq_len(nrow(x0)) %in% c(1, nrow(x0))

z0 <-  (x0[-nrow(x0),1:2] + x0[-1,1:2])/2
model <- curve.model(d1$gmt, d1$light, d1$segment, calibfun,
                     alpha = c(7, 10), beta = c(150, 80),
                     logp.x = log.prior,
                     logp.z = log.prior,
                     x0 = x0,z0 = z0,fixedx = fixedx)

x.proposal <- mvnorm(S=diag(c(0.005,0.005, 0.005)),n=nrow(x0))
z.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nrow(x0)-1)

fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=100,thin=20)

## testthat
## fixedx[i/n] is really a logical vector of length n, no NAs (and friends)
## dt is all valid, no funky
## all arguments are right length, sign and mode
## calib returns valid values that are increasing over -10:10


all(fixedx | is.sea(chain.last(fit$x)))
all(is.sea(chain.last(fit$z)))

x.proposal <- mvnorm(chain.cov(fit$x),s=0.3, tol = 1e-02)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)

x.proposal <- mvnorm(chain.cov(fit$x),s=0.3, tol = 1e-02)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)




x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)


opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
k <- sample(nrow(x0),20)
matplot(t(fit$x[k,1,]),type="l",lty=1,col="dodgerblue",ylab="Lon")
matplot(t(fit$x[k,2,]),type="l",lty=1,col="firebrick",ylab="Lat")
par(opar)

x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)




## Tune proposals based on previous run

x.proposal <- mvnorm(chain.cov(fit$x),s=0.3, tol = 1e-05)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)

fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=3000,thin=20)

