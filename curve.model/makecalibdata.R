
library(SGAT)
if ( setInternet2(NA)) {
    ##source("C:\\Users\\michae_sum\\Documents\\GitHub\\SGAT\\R\\Curve.R")
    load("D:\\data\\TAGDATA\\mcmcGeo\\b362_99\\readyRun.Rdata")
    load("D:\\data\\TAGDATA\\mcmcGeo\\b362_99\\old\\metropolisMarch06\\x.Rdata")
} else {
##    source("C:\\Users\\mdsumner\\Documents\\GitHub\\SGAT\\R\\Curve.R")
    load("E:\\backup\\Elements\\Azimuth\\DATA\\mcmcGeo\\b362_99\\readyRun.Rdata")
load("E:\\backup\\Elements\\Azimuth\\DATA\\mcmcGeo\\b362_99\\old\\metropolisMarch06\\x.Rdata")
}


## don't ask
body(calib) <- parse(text = ".approxfun(x, y, v, method, yleft, yright, f)")
##mkOldNew <- function(calibfun) {
##    calibfun
##    function(x) calibfun(90 - x)
##}

##newcalib <- mkOldNew(calib)

elev <- seq(-12, 12, by = 0.2)
ellie1calib <- data.frame(zenith = 90 - elev, light = calib(elev))

## actual calibration now is
##calib <- approxfun(elliecalib$zenith, elliecalib$light)
##with(elliecalib, plot(zenith, light))
##with(elliecalib, lines(zenith, calib(zenith)))

## data
d1 <- d[!is.na(d$segment), ]
d1 <- d1[!is.na(d1$light), ]
d1$segment <- unclass(factor(d1$segment))
ellie1 <- d1[d1$depth < 15, c("gmt", "light", "segment")]
names(ellie1)[1] <- "time"
rownames(ellie1) <- seq_len(nrow(ellie1))

save(ellie1calib, file = "ellie1calib.rda")
save(ellie1, file = "ellie1.rda")

