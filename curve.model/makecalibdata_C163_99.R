
library(SGAT)


load("E:\\backup\\Elements\\Azimuth\\DATA\\mcmcGeo\\c699_02\\readyRun.Rdata")




## don't ask
body(calib) <- parse(text = ".approxfun(x, y, v, method, yleft, yright, f)")
##mkOldNew <- function(calibfun) {
##    calibfun
##    function(x) calibfun(90 - x)
##}

##newcalib <- mkOldNew(calib)

elev <- seq(-12, 12, by = 0.2)
ellie2calib <- data.frame(zenith = 90 - elev, light = calib(elev))

## actual calibration now is
##calib <- approxfun(elliecalib$zenith, elliecalib$light)
##with(elliecalib, plot(zenith, light))
##with(elliecalib, lines(zenith, calib(zenith)))

## data
##d1 <- d[!is.na(d$segment), ]
##d1 <- d1[!is.na(d1$light), ]
##d1$segment <- unclass(factor(d1$segment))
ellie2 <- d1[d1$depth < 15, c("gmt", "light", "depth", "temp", "segment")]
names(ellie2)[1] <- "time"
rownames(ellie2) <- seq_len(nrow(ellie2))

save(ellie2calib, file = "ellie2calib.rda")
save(ellie2, file = "ellie2.rda")

