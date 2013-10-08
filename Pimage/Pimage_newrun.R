load("D:\\projects\\SGAT\\Movement_fit.Rdata")
load("E:\\DATA\\SGAT_\\Movement_fit.Rdata")

fit$model$time <- fit$model$twilight

library(SGAT)


## easy spec for a local projection
x <- model.bin(fit, proj = "laea")

## fuller spec, not very sensible
x <- model.bin(fit, bin = "intermediate", proj = "+proj=utm +south +zone=50")

## default grids, different binning
xprimary <- model.bin(fit)
xintermediate <- model.bin(fit, bin = "intermediate")

## these should be equal
cellStats(xprimary[], sum) == prod(dim(fit$x)[-2])
sum(values(xintermediate[])) == difftime(max(fit$model$time), min(fit$model$time), units = "hours")

## accumulate
fit$x[80,1,] <- fit$x[80,1,] + 2
fit$z[80,1,] <- fit$z[80,1,] - 2
xprimary2 <- model.bin(fit, pimg = xprimary)
xintermediate2 <- model.bin(fit, bin = "intermediate", pimg = xintermediate)

## count is doubled
cellStats(xprimary2[], sum)
## time is the same
sum(values(xintermediate2[]))


## ensure that we cannot mix the bin type and Pimage type
xx <- model.bin(fit, pimg = xprimary, bin = "intermediate")
##Error in t(chain[k, 1:2, ]) :
##  error in evaluating the argument 'x' in selecting a method for function 't': Error in chain[k, 1:2, ] : subscript out of bounds
xx <- model.bin(fit, pimg = xintermediate, bin = "primary")




## default whole-world grid
x <- model.bin(fit, grid = raster())

## old-school
img <- list(x = 101:180, y = -55:-21, z = matrix(0, 80, 35))
x <- model.bin(fit, bin = "intermediate", grid = img)

## grid with a projection
prj <-  "+proj=ortho +lon_0=140 +lat_0=-80"
img.p <-  projectExtent(raster(img), prj)
x <- model.bin(fit, grid = img.p)


## try accumulating one
x1 <- model.bin(fit, bin = "intermediate")
x2 <- x1
par(mfrow = c(4, 2))
plot(x1[], col = trip::oc.colors(256))

fittest <- fit
for (i in seq(-3, 3, by = 1)) {

    fittest$z[80,1,] <- fit$z[80,1,] + i
    x2 <- model.bin(fittest, bin = "intermediate", pimg = x2)
    plot(x2[], col = trip::oc.colors(256))
    ##print(sum(x2[][]))
    print(sum(x2[[80]]$image))
}



## try multiple
load("D:\\projects\\SGAT\\Shearwater_fit.Rdata")
fit$model$time <- fit$model$twilight
y <- model.bin(fit, grid = x[1])

multibin <- c(x, y)




img <- list(x = 101:180, y = -55:-21, z = matrix(0, 80, 35))
x <- model.bin(fit, bin = "intermediate", grid = img)

 x[[3]]
$xbound
[1] 101 180  80

$ybound
[1] -55 -21  35

$offset
[1] 29 24

$image
         [,1]       [,2]       [,3]       [,4]       [,5]
[1,]   0.0000    0.00000   11.56764   11.56764    0.00000
[2,] 104.1088 4233.75599 9404.49077 7542.10084 1700.44298
[3,]   0.0000   11.56764   34.70292   69.40584   11.56764

attr(,"class")
[1] "pimg" "list"
>
>
>
>
>
> diff(fit$model$time)[3]
Time difference of 11.56764 hours
> sum(x[[3]]$image)
[1] 23135.28
> sum(x[[3]]$image)/2000
[1] 11.56764
