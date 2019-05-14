subcloneDensityEst <- function(GS.data, density.file = "density.txt", density.smooth = 0.1,
                                       density.bw = "nrd0", burn.in = 500, y.max = 5){
  
  ## Input parameters
  # GS.data is the list output from function -- subcloneDirichletGibbs()
  # density.smooth is the parameter used in R's density() function
  
  V.h <- GS.data$V.h
  Phi.h <- GS.data$Phi.h
  no.iter <- nrow(V.h)
  wts <- matrix(NA, nrow = no.iter, ncol = ncol(V.h))
  wts[,1] <- V.h[,1]
  wts[,2] <- V.h[,2] * (1-V.h[,1])
  if(ncol(V.h) > 2) {
    for (i in 3:dim(wts)[2]) 
      wts[,i] <- apply(1-V.h[,1:(i-1)], MARGIN=1, FUN=prod) * V.h[,i]
  }
  xx <- density(c(Phi.h[burn.in-1,]), weights=c(wts[burn.in,]) / sum(c(wts[burn.in,])), bw = density.bw, adjust=density.smooth, from=0, to=1)$x
  
  post.density <- matrix(NA, ncol = no.iter-burn.in+1, nrow = 512) # nrow is equal to the output of function density()
  for(i in burn.in:no.iter)
    post.density[,i - burn.in + 1] <- density(c(Phi.h[i-1,]), weights=c(wts[i,]) / sum(c(wts[i,])), bw = density.bw, adjust=density.smooth, from=0, to=1)$y
  
  # Take the median of density
  yy <- apply(post.density, MARGIN = 1, FUN = quantile, probs = 0.5)
  write.table(cbind(xx,yy),density.file,sep="\t",col.names=c("mutation.burden","median.density"),row.names=F,quote=F)
  
}

