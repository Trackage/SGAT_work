




 trim0<- function(x, ...) {
   x[!x > 0] <- NA
   trim(x)
 }

## validity
##  see spatstat/R/verifyclass.S for inspiration
##  any or all NULL images?
##  all pimg objects are consistent
##  all offsets are sensible, etc.


print.Pimage <- function(x, ...) {
  ## this needs to know the x/y/time range, and possibly the sizes of all images, whether any are NULL or funny
    ext <- extent(as.raster(x))
    trange <- format(range(attr(x, "times")))
    Z <- .Z(x)
    cat("Class    :", class(x), c("(Primary/X)", "(Intermediate/Z)")[Z + 1], "\nLength    :", length(x),  "\nTemporal Extent :", trange, "\n")
    ##cat("Time Steps   :")
    ##str(attr(x, "times"))
    print(ext)
    invisible(NULL)
}
summary.Pimage <- function(object, ...) {
  summary(getValues(as.raster(object)))
}

as.POSIXct.Pimage <- function(x, tz = "", ...) {
  .times(x)
}
.times <- function(x) {
  UseMethod(".times")
}
.times.Pimage <- function(x) {
    out <- attr(x, "times")
    if (.Z(x)) out <- out[-length(out)]
    out
}

## must be generic for replacement method to work
".times<-" <- function(x, value) {
  UseMethod(".times<-")
}
".times<-.Pimage" <- function(x, value) {
  attr(x, "times") <- value
  x
}
.Z <- function(x) {
  UseMethod(".Z")
}
.Z.Pimage <- function(x) {
  attr(x, "Z")
}
".Z<-" <- function(x, value) {
  UseMethod(".Z<-")
}
".Z<-.Pimage" <- function(x, value) {
  attr(x, "Z") <- value
  x
}

is.Pimage <- function(x) {inherits(x, "Pimage")}

names.Pimage <- function(x) {
    format(as.POSIXct(x))
}

  ###str.Pimage

plot.Pimage <- function(x, ...) {
  plot(x[seq_along(x)], ...)
}
###"[<-.Pimage"
##"$.Pimage"
##"$<-.Pimage"
#head.Pimage
#tail.Pimage
#range.Pimage
#scale.Pimage ## divide out iterations?
##with.Pimage
##subset.Pimage
#names.Pimage ## formatted date-times
## "names<-.Pimage"
#rev.Pimage  ## should only be in order
## ?? ifelse.Pimage
#median.Pimage
#quantile.Pimage
## transform.Pimage


## coercion
as.Pimage <- function(x, ...) {
    UseMethod("as.Pimage")
}

as.Pimage.Pimage <- function(x) {
    x
}

as.raster.Pimage <- function(x) {
    ## TODO, patch in the data, like tripEstimation::combine()
    ##  include message on object about its origin
    x <- x[[1]]
    raster(nrows = x$xbound[3L], ncols = x$ybound[3L], xmn = x$xbound[1L], xmx = x$xbound[2L], ymn = x$ybound[1L], ymx = x$ybound[2L])
}

bin <- function(x, ...) UseMethod("bin")

bin.Pimage <- function (pimgs, chain, weights = NULL)
{
  if (is.null(weights) && .Z(pimgs)) {
    weights <- c(diff(unclass(.times(pimgs))/3600))
  }
  if (is.null(weights)) weights <- rep(1, length(pimgs))


  if (!(length(weights) == length(pimgs))) stop("length of weights do not match length of p-img list")
  if (nrow(chain) != length(pimgs))
    stop("dimensions of chain do not match length of p-img list")
  #dm <- dim(z)


  for (k in seq_along(weights)) {
    pimgs[[k]] <- bin.pimg(pimgs[[k]], t(chain[k, 1:2, ]), weight = weights[k])
  }
  ## should also have a flag for whether this is initialized/scaled, so iter number is independent
  attr(pimgs, "itersbin") <- attr(pimgs, "itersbin") + dim(chain)[3]
  pimgs
}




