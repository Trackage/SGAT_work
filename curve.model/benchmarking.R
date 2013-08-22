require(rbenchmark)

nn <- sample(10000, replace = TRUE)
nn <- nn[nn > 100]
l <- lapply(nn, rnorm)


benchmark(ul =  unlist(lapply(l, sum)),
          ulnr = unlist(lapply(l, sum), recursive = FALSE),
          dc = do.call(c, lapply(l, sum)),
          replications = 5000)


