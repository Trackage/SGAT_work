library(SGAT)

if ( setInternet2(NA)) {
    ##source("C:\\Users\\michae_sum\\Documents\\GitHub\\SGAT\\R\\SGAT.R")
    load("D:\\data\\TAGDATA\\mcmcGeo\\b362_99\\readyRun.Rdata")
    load("D:\\data\\TAGDATA\\mcmcGeo\\b362_99\\old\\metropolisMarch06\\x.Rdata")
} else {
    ##source("C:\\Users\\mdsumner\\Documents\\GitHub\\SGAT\\R\\SGAT.R")
    load("E:\\backup\\Elements\\Azimuth\\DATA\\mcmcGeo\\b362_99\\readyRun.Rdata")
load("E:\\backup\\Elements\\Azimuth\\DATA\\mcmcGeo\\b362_99\\old\\metropolisMarch06\\x.Rdata")
}


## don't ask
body(calib) <- parse(text = ".approxfun(x, y, v, method, yleft, yright, f)")
mkOldNew <- function(calibfun) {
    calibfun
    function(x) calibfun(90 - x)
}

newcalib <- mkOldNew(calib)

##Ellieb362_99 <- d1[,c(2, 3, 6)]
##calibdata <- data.frame(zenith = seq(100, 80, by = 0.2), light = newcalib(seq(100, 80, by = 0.2)))


d1 <- d[!is.na(d$segment), ]
d1 <- d1[d1$depth < 20, ]
d1 <- d1[!is.na(d1$light), ]
d1$segment <- unclass(factor(d1$segment))



x0 <- x

fixedx <- seq_len(nrow(x0)) %in% c(1, nrow(x0))

z0 <-  (x0[-nrow(x0),1:2] + x0[-1,1:2])/2
model <- curve.model(d1$gmt, d1$light, d1$segment, newcalib,
                     alpha = c(7, 10), beta = c(0.5, 0.08),
                     x0 = x0,z0 = z0,
                     fixedx = fixedx)

## find starting points
grid <- list(x = seq(lon.min, lon.max, length = 40),
             y = seq(lat.min, lat.max, length = 30),
             z = array(0, c(40, 30, nrow(x0))))
for (i in seq_along(grid$x)) {
        for (j in seq_along(grid$y)) {
            grid$z[i,j,] <- model$logpx(cbind(rep(grid$x[i], nrow(x0)), grid$y[j], 0))
        }
    }

x0 <- cbind(as.matrix(expand.grid(grid$x, grid$y)[apply(grid$z, 3, which.max), ]), 0)

x.proposal <- mvnorm(S=diag(c(0.005,0.005, 0.05)),n=nrow(x0))
z.proposal <- mvnorm(S=diag(c(0.005,0.005)),n=nrow(x0)-1)
##fit <- stella.metropolis(model,x.proposal,iters=100,thin=20, x0 = x0)

fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=100,thin=20, x0 = x0, z0 = (x0[-nrow(x0),1:2]+ x0[-1,1:2])/2)

## testthat
## fixedx[i/n] is really a logical vector of length n, no NAs (and friends)
## dt is all valid, no funky
## all arguments are right length, sign and mode
## calib returns valid values that are increasing over -10:10


all(fixedx | is.sea(chain.last(fit$x)))
all(is.sea(chain.last(fit$z)))

x.proposal <- mvnorm(chain.cov(fit$x),s=0.3, tol = 1e-06)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)


##fit <- stella.metropolis(model,x.proposal,
##                          x0=chain.last(fit$x),
##                          iters=300,thin=20)





fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)

x.proposal <- mvnorm(chain.cov(fit$x),s=0.3, tol = 1e-05)
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
