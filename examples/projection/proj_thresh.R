## from the SGAT Movement vignette

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
p0 <- c(130,-50)
p1 <- c(140,-80)
p2 <- c(130,-60)
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

library(SGAT)
## Resample to two minute intervals and compute zenith angles
tms <- as.POSIXct("2013-04-01","GMT")+(24*60*60)*days
tms.out <- seq(min(tms),max(tms),by=2*60)
d <- zenith.simulate(tms,ps[,1],ps[,2],tms.out)

## Compute times and locations at which twilight is observed.
twl <- twilight.simulate(d)
twl <- twilight.perturb(twl,rlnorm(nrow(twl),0.3,0.6))

## Initial x locations
x0 <- threshold.path(twl$twilight,twl$rise)$x



model <- threshold.model(twl$twilight,twl$rise,
                         twilight.model="LogNormal",
                         alpha=c(0.3,0.6),beta=c(8,3.5),
                         x0=x0,z0=NULL)




r <- model$residuals(x0)
x0 <- cbind(x0[,1]-ifelse(twl$rise,1,-1)*pmin(0,r-0.1)/1440*360,x0[,2])
library(geosphere)
z0 <- midPoint(x0[-nrow(x0),],x0[-1,])
model <- threshold.model(twl$twilight,twl$rise,
                         twilight.model="LogNormal",
                         alpha=c(0.3,0.6),beta=c(8,3.5),
                         x0=x0,z0=z0)



chx <- NULL
chz <- NULL

## Define initial proposals
x.proposal <- mvnorm(S=diag(c(0.002,0.002)),n=nrow(twl))
z.proposal <- mvnorm(S=diag(c(0.002,0.002)),n=nrow(twl)-1)

## Short test run
fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=300,thin=20)

ith <- 59

chx <- rbind(chx, t(fit$x[ith, , ]))
chz <- rbind(chz, t(fit$z[ith, , ]))

## Tune proposals based on previous run
x.proposal <- mvnorm(chain.cov(fit$x,discard=100),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z,discard=100),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)

chx <- rbind(chx, t(fit$x[ith, , ]))
chz <- rbind(chz, t(fit$z[ith, , ]))

## Tune proposals based on previous run
x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)


chx <- rbind(chx, t(fit$x[ith, , ]))
chz <- rbind(chz, t(fit$z[ith, , ]))

x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=3000,thin=20)

chx <- rbind(chx, t(fit$x[ith, , ]))
chz <- rbind(chz, t(fit$z[ith, , ]))




ith <- 59
opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
matplot(scale(t(fit$x[ith,,]),scale=F),type="l",lty=1,col=c(2,4),
        xlab="",ylab=expression(x[ith]))
matplot(scale(t(fit$z[ith,,]),scale=F),type="l",lty=1,col=c(2,4),
        xlab="",ylab=expression(z[ith]))
par(opar)
