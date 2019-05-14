subcloneDirichletGibbs <- function(C = 10+1, a, d, copyNumberCoef = rep(2, length(a)), iter = 1000) {

  ## Input parameters
  # C: maximum number of subclones 
  # a: a vector of variant read counts
  # d: a vector of total read counts
  # iter: number of iteration of Gibbs sampler
  
  N <- length(a)
  
  # Hyperparameters for alpha
  A <- B <- 0.01
  
  # Set up data formats for recording iterations
  alpha <- rep(NA, iter)
  pi.h <- matrix(NA, nrow = iter, ncol = C)     # {Phi_k} for k = 1,...,K, i.e., cellular proportion of each cluster
  V.h <- matrix(NA, nrow = iter, ncol = C)      # stick-breaking
  z.i <- matrix(NA, nrow = iter, ncol = N)      # data assignments
  Pr.z <- matrix(NA, nrow = N, ncol = C)        # Pr(z.n = k), used for sample {zn}
  mutBurden <- array(NA, dim = c(iter, C, N))   # expected fraction of reads containing mutation, used as p in Binmomial(n,p)
  
  # Randomize starting positions
  alpha[1] <- 1
  V.h[1,] <- c(rep(0.5, C-1), 1)  # Truncated Dirichlet Process(TDP): last element is set to 1 so that pi.h[,k] = 0 for all k>C
  V.h[1:iter, C] <- rep(1, iter)
  z.i[1,] <- sample.int(C, size = N, replace = TRUE)
  
  lower <- min(copyNumberCoef*a/d)
  upper <- max(copyNumberCoef*a/d)
  difference <- upper - lower
  pi.h[1,] <- runif(C, lower<-(lower-difference/10), upper<-(upper+difference/10))
  for(c in 1:C){
    tmp <- pi.h[1,c]/copyNumberCoef
    tmp[which(tmp<0.000001)] <- 0.000001
    tmp[which(tmp>0.999999)] <- 0.999999
    mutBurden[1,c,] <- tmp
  }
  
  # Update cluster assignment for each mutation
  for(i in 2:iter) {
    
    # independently sample zn from Pr(zn=k)
    for(n in 1:N) {
      Pr.z[n,1] <- log(V.h[i-1,1]) + dbinom(a[n],d[n],prob = mutBurden[i-1,1,n], log = TRUE)
      Pr.z[n,2:C] <- sapply(2:C, function(j) {log(V.h[i-1,j])+sum(log(1-V.h[i-1, 1:(j-1)])) + dbinom(a[n],d[n],prob = mutBurden[i-1,j,n], log = TRUE)})
      Pr.z[n,] <- exp(Pr.z[n,] - logSumExp(Pr.z[n,],na.rm = TRUE))  
      Pr.z[n,is.na(Pr.z[n,])] = 0
    }
    z.i[i,] <- sapply(1:N, function(j) {sum(rmultinom(1,1,Pr.z[j,]) * (1:C))})
    
    # independently sample V.h
    V.h[i,1:(C-1)] <- sapply(1:(C-1), function(j) {rbeta(1, 1+sum(z.i[i,] == j), alpha[i-1]+sum(z.i[i,] > j))})
    V.h[i,which(V.h[i,1:(C-1)] == 1)] <- 0.999    # prevent one stick from taking all the weight
    
    # independently sample pi.h from the posterior distribution in the same family as the base distribution
    pi.h[i,] <- runif(C, lower, upper)
    mutBurden[i,,] <- mutBurden[i-1,,]
    for(c in unique(z.i[i,])) {
      # the mean of gamma distribution is shape/rate.
      pi.h[i,c] <- rgamma(1, shape = sum(a[z.i[i,] == c]), rate = sum((d/copyNumberCoef)[z.i[i,] == c]))
      tmp <- pi.h[i,c]/copyNumberCoef
      tmp[which(tmp<0.000001)] <- 0.000001
      tmp[which(tmp>0.999999)] <- 0.999999
      mutBurden[i,c,] <- tmp
    }
    
    # update alpha
    alpha[i] <- rgamma(1, shape = A+C-1, rate = B-sum(log(1-V.h[i,1:(C-1)])))
    
  }
  return(list(z.i=z.i, V.h = V.h, Phi.h = pi.h, alpha = alpha, a = a, d = d))
}
