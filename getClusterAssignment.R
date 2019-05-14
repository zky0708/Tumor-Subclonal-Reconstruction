getClusterAssignment <- function(GS.data, density.file = "density.txt", window.size = 20, burn.in = 500, density.threshold = 0.01) {
  
  mutReads <- GS.data$a
  totalReads <- GS.data$d
  no.muts = length(mutReads)
  
  z.i = GS.data$z.i
  V.h = GS.data$V.h
  Phi.h = GS.data$Phi.h
  iter = nrow(z.i)
  
  density = read.table(density.file, row.names = NULL, header = TRUE, sep = "\t")
  
  # Use sliding window to find local optima
  localOptima = NULL
  peak.indices = NULL
  for(i in (1+window.size):(nrow(density)-window.size)) {
    if(density$median.density[i] == max(density$median.density[(i-window.size):(i+window.size)])) {
      localOptima = c(localOptima, density[i,1])
      peak.indices = c(peak.indices, i)
    }
  }
  
  # Deal with cases where no local optima is detected
  if(is.null(localOptima)) {
    peak.indices <- which.max(density$median.density[(1+window.size):(nrow(density)-window.size)])
    localOptima <- density[peak.indices, 1]
  }

  # Remove peaks with very small heights
  densities <- density$median.density[peak.indices]
  id <- which((densities / sum(densities)) > density.threshold)
  localOptima <- localOptima[id]
  peak.indices <- peak.indices[id]
  
  # Assign mutations to clusters
  no.optima = length(localOptima)
  if(no.optima>1){
    boundary = array(NA, no.optima-1)
    for(i in 1:(no.optima-1)) {
      min.density = min(density$median.density[(peak.indices[i]+1):(peak.indices[i+1]-1)])
      min.indices = intersect(which(density$median.density == min.density), (peak.indices[i]+1):(peak.indices[i+1]-1))
      # positions of each minimum between pairs of optima
      boundary[i] = (density[max(min.indices),1] + density[min(min.indices),1])/2
    }
    
    sampledIters = (burn.in+1):iter
    # do not use the initial state
    sampledIters = sampledIters[sampledIters != 1]
    if(length(sampledIters) > 1000) 
      sampledIters = floor(burn.in + (1:1000)*(iter - burn.in) / 1000)
    
    z.i = data.matrix(z.i)
    mutation.preferences = array(0, c(no.muts, no.optima))
    for(s in sampledIters) {
      temp.preferences = array(0, c(no.muts, no.optima))
      for(c in unique(z.i[s,])){
        # use the relative position between Phi.h and boundary to determine the cluster it belongs to
        bestOptimum = sum(Phi.h[s,c] > boundary)+1
        temp.preferences[z.i[s,]==c, bestOptimum] = temp.preferences[z.i[s,]==c, bestOptimum] + 1
      }
      iter.preferences = t(sapply(1:no.muts, function(p,i) {as.integer(p[i,]==max(p[i,]))/sum(p[i,]==max(p[i,]))}, p = temp.preferences))
      mutation.preferences = mutation.preferences + iter.preferences
    }
    mutation.preferences = mutation.preferences/length(sampledIters)
    most.likely.cluster = sapply(1:no.muts, function(i) which.max(mutation.preferences[i,]))
  } else {
    most.likely.cluster = rep(1, no.muts)
  }
  
  return(list(localOptima = localOptima, most.likely.cluster = most.likely.cluster))
  
}
