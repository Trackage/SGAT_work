load("B362_99.rda")
d$gmt <- seq(as.POSIXct("1999-10-27 07:07:30", tz = "GMT"), length = nrow(d), by = 30)

