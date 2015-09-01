library(BAStag)

d <- readMTlux("A042_21Mar15_221532.lux")
#twl <- preprocessLight(d, threshold = 96, offset = 16)
#save(twl, file = "twl.RData")
load("twl.RData")
cr <- findCrepuscular(d, twl, min.threshold = 10, max.threshold = 20000)


d$Light <- log(d$Light)

ec <- extractCrepuscular(d, cr)
ec <- subset(ec, Light > 0)

## calibration
calibLoc <- c(-58.83584, -52.21008)
calibDate <- as.POSIXct(c("2014-03-21 00:00:00", "2014-03-25 00:00:00"), tz = "GMT")
calibData <- subset(ec, Date <= calibDate[2])
calibData$Zenith <- zenith(solar(calibData$Date), calibLoc[1], calibLoc[2])
library(mgcv)
fit <- gam(Light ~ s(Zenith, k = 20) + Segment, family = gaussian(link = "identity"), data = calibData)
calibPredict <- with(calibData, data.frame(Zenith = seq(min(Zenith), max(Zenith), length = 500), 
                                      Segment = Segment[1]))
calibPredict$Light <- predict(fit, newdata, type = "response")
calibration <- approxfun(calibPredict$Zenith, calibPredict$Light, rule = 2)



library(raster)
library(maptools)
data(wrld_simpl)


llex <- extent(rep(calibLoc, each = 2) + c(-1, 1, -1, 1) * 10)

#g <- setValues(raster(llex, nrow = 80, ncol = 160), 1)
g <- rasterize(wrld_simpl, raster(llex, nrow = 80, ncol = 160))
g[!is.na(g)] <- 0
g[is.na(g)] <- 1


#x0 <- thresholdLocation(twl$Twilight, twl$Rise)$x
x0 <- thresholdPath(twl$Twilight, twl$Rise)$x
#x0 <- cbind(x0, 0)

model <- essieCurveModel(ec$Date, ec$Light, ec$Segment, calibration,
                         alpha = c(2.2,1.0),
                         beta = c(2.8, 0.12), 
                         x0 = x0)
                     #alpha = c(7, 10), beta = c(8, 3.5), x0 = x0)

###################################################
fit <- essie(model, grid = g)
