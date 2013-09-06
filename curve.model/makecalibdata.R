
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

elev <- seq(-12, 12, by = 0.2)
calib.data <- data.frame(elevation = elev, light = calib(elev), zenith = 90 - elev)
write.csv(calib.data, file = "C:\\Users\\michae_sum\\Documents\\GitHub\\SGAT_work\\curve.model\\calib.data.csv", row.names = FALSE)

