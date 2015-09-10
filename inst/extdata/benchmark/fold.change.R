library(chemhelper)
nvars <- 1000
nsamples <- 50
sample_groups <- 5
data <- replicate(nvars, runif(n=nsamples))

f <- rep_len(1:sample_groups, nsamples)
f <- LETTERS[f]


system.time({
change1 <- fold.change(MAT=data,f=f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y),method=1)
})

system.time({
  change2 <- fold.change(MAT=data,f=f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y),method=2)
})

system.time({
  change3 <- fold.change(MAT=data,f=f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y),method=3)
})


system.time({
  change4 <- fold.change(MAT=data,f=f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y),method=4)
})



identical(change1,change2)
identical(change1,change3)
identical(change1,change4)




