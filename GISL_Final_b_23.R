############################# FUNCTIONS ################################
subclone.dirichlet.gibbs <- function(C = 10+1, a, d, copyNumberCoef = rep(2, length(a)), iter = 1000) {
  # C: maximum number of subclones 
  # a: a vector of variant read counts
  # d: a vector of total read counts
  # iter: number of iteration of Gibbs sampler
  
  N <- length(a)
  # print(paste("num.muts=", N, sep=""))
  
  # Hyperparameters for alpha
  A <- B <- 0.01
  
  # set up data formats for recording iterations
  alpha <- rep(NA, iter)
  pi.h <- matrix(NA, nrow = iter, ncol = C)     # {Phi_k} for k = 1,...,K, i.e., cellular proportion of each cluster
  V.h <- matrix(NA, nrow = iter, ncol = C)      # stick-breaking
  z.i <- matrix(NA, nrow = iter, ncol = N)      # data assignments
  Pr.z <- matrix(NA, nrow = N, ncol = C)        # Pr(z.n = k), used for sample {zn}
  mutBurden <- array(NA, dim = c(iter, C, N))   # expected fraction of reads containing mutation, used as p in Binmomial(n,p)
  
  # randomize starting positions
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
  
  # update cluster assignment for each mutation
  for(i in 2:iter) {
    # if(i %% 100 == 0) print(i)
    
    # independently sample zn from Pr(zn=k)
    for(n in 1:N) {
      # use log-space
      # Pr.z[n,1] <- log(V.h[i-1,1]) + a[n]*log(mutBurden[i-1,1,n]) + (d[n]-a[n])*log(1-mutBurden[i-1,1,n])
      Pr.z[n,1] <- log(V.h[i-1,1]) + dbinom(a[n],d[n],prob = mutBurden[i-1,1,n], log = TRUE)
      # Pr.z[n,1] <- log(V.h[i-1,1]) + dnbinom(a[n],d[n],prob = mutBurden[i-1,1,n], log = TRUE)
      # Pr.z[n,2:C] <- sapply(2:C, function(j) {log(V.h[i-1,j])+sum(log(1-V.h[i-1, 1:(j-1)])) + a[n]*log(mutBurden[i-1,j,n]) + (d[n]-a[n])*log(1-mutBurden[i-1,j,n])})
      Pr.z[n,2:C] <- sapply(2:C, function(j) {log(V.h[i-1,j])+sum(log(1-V.h[i-1, 1:(j-1)])) + dbinom(a[n],d[n],prob = mutBurden[i-1,j,n], log = TRUE)})
      # Pr.z[n,2:C] <- sapply(2:C, function(j) {log(V.h[i-1,j])+sum(log(1-V.h[i-1, 1:(j-1)])) + dnbinom(a[n],d[n],prob = mutBurden[i-1,j,n], log = TRUE)})
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
      # mutBurden[i,c,] <- pi.h[i,c]/copyNumberCoef
      tmp <- pi.h[i,c]/copyNumberCoef
      tmp[which(tmp<0.000001)] <- 0.000001
      tmp[which(tmp>0.999999)] <- 0.999999
      mutBurden[i,c,] <- tmp
      # mutBurden[i,c, mutBurden[i,c,]>0.499999] <- 0.499999
    }
    
    # update alpha
    alpha[i] <- rgamma(1, shape = A+C-1, rate = B-sum(log(1-V.h[i,1:(C-1)])))
    
  }
  return(list(z.i=z.i, V.h = V.h, Phi.h = pi.h, alpha = alpha, a = a, d = d))
}

Gibbs.subclone.density.est <- function(GS.data, density.file = "density.txt", density.smooth = 0.1,
                                       density.bw = "nrd0", burn.in = 500, y.max = 5){
  # GS.data is the list output from function -- subclone.dirichlet.gibbs()
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
  
  # take the median of density
  yy <- apply(post.density, MARGIN = 1, FUN = quantile, probs = 0.5)
  write.table(cbind(xx,yy),density.file,sep="\t",col.names=c("mutation.burden","median.density"),row.names=F,quote=F)
}


getClusterAssignment <- function(GS.data, density.file = "density.txt", window.size = 20, burn.in = 500, density.threshold = 0.01) {
  mutReads <- GS.data$a
  totalReads <- GS.data$d
  no.muts = length(mutReads)
  
  z.i = GS.data$z.i
  V.h = GS.data$V.h
  Phi.h = GS.data$Phi.h
  iter = nrow(z.i)
  
  density = read.table(density.file, row.names = NULL, header = TRUE, sep = "\t")
  
  # use sliding window to find local optima
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
  
  # assign mutations to clusters
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


## MAIN CODE ##

args <- commandArgs(TRUE)
suppressMessages(require("matrixStats"))
suppressMessages(require("data.table"))
options(warn = -1)

########## Subchallenge 1 ##########

## Read MuTect somatic mutation files
ssmdat = read.table(args[1],sep='\t',comment.char='#', stringsAsFactors = FALSE)
no.muts <- nrow(ssmdat)
colnames(ssmdat) = c("CHROM","POS","ID","REF", "ALT", "QUAL","FILTER", "INFO","FORMAT","normal", "tumor")
rownames(ssmdat) = paste(ssmdat[,"CHROM"], ssmdat[,"POS"], sep = "_")
# Make sure the characters of CHROM are of all upper case
ssmdat$CHROM <- toupper(ssmdat$CHROM)
# matched normal sample
normal_stat = do.call(rbind, strsplit(as.vector(ssmdat[,"normal"]), split = ":", fixed = TRUE))
colnames(normal_stat) = strsplit(as.vector(unique(ssmdat[, "FORMAT"])), ":")[[1]]
rownames(normal_stat) <- rownames(ssmdat)
# tumor sample
tumor_stat = do.call(rbind, strsplit(as.vector(ssmdat[,"tumor"]), split = ":", fixed = TRUE))
colnames(tumor_stat) = strsplit(as.vector(unique(ssmdat[, "FORMAT"])), ":")[[1]]
rownames(tumor_stat) <- rownames(ssmdat)
# Get read counts of mutation reads and total reads
mutReads <- as.numeric(unlist(lapply(strsplit(as.vector(tumor_stat[,"AD"]),','), '[[', 2)))
totalReads <- as.integer(as.vector(tumor_stat[,'DP']))

## Read Battenberg CNV files
cnvdat <- read.table(args[2], header = TRUE, stringsAsFactors = FALSE)
rownames(cnvdat) = paste(cnvdat[,1], cnvdat[,2], cnvdat[,3], sep = "_")

# Use Battenberg cellularity for SC2 and SC3 
cellularity0 <- read.table(args[3], header = TRUE)$cellularity

# Map mutations to CNV intervals according to their positions
# cnv_ssm = vector("list", nrow(cnvdat)-length(which(cnvdat$chr %in% c("X","Y"))))
cnv_ssm = vector("list", nrow(cnvdat))
names(cnv_ssm) = rownames(cnvdat)[1:length(cnv_ssm)]
for(i in 1:length(cnv_ssm)){
  chr = cnvdat[i,"chr"]
  start = cnvdat[i,"startpos"]
  end = cnvdat[i,"endpos"]
  id = rownames(ssmdat)[which(ssmdat$CHROM == chr)]
  if(chr == "X")
    id <- c(id, rownames(ssmdat)[which(ssmdat$CHROM == "Y")])
  for(j in id)
    if(ssmdat[j,"POS"] >= start && ssmdat[j,"POS"] <= end)
      cnv_ssm[[i]] = c(cnv_ssm[[i]], j)
}

# Get the normal and tumor copy number information for each SSM
# Initialize with normal state
copyNumber.ssm = matrix(NA, nrow = nrow(ssmdat), ncol = 5)
rownames(copyNumber.ssm) = rownames(ssmdat)
colnames(copyNumber.ssm) = c("CNV","normal_cn","tumor_cn", "subclonal", "coef1")
copyNumber.ssm[,"normal_cn"] <- 2 
copyNumber.ssm[,"tumor_cn"] <- 2
copyNumber.ssm[,"subclonal"] <- 0
# Check gender
id.xy <- NULL
if(length(which(ssmdat$CHROM == "Y")) > 0) {
  # Male, so the normal copy number of X nd Y chromosome is 1
  id.xy <- which(ssmdat$CHROM %in% c("X","Y"))
  copyNumber.ssm[id.xy, "normal_cn"] <- 1
  copyNumber.ssm[id.xy,"tumor_cn"] <- 1
}

# Get copy number abbreviation in tumor cells
for(i in 1:length(cnv_ssm)){
  ssms <- cnv_ssm[[i]]
  if(length(ssms) == 0)
    next
  copyNumber.ssm[ssms, "CNV"] <- i
  # copyNumber.ssm[ssms, "tumor_cn"] <- round(cnvdat[i,"ntot"])
  copyNumber.ssm[ssms, "tumor_cn"] <- cnvdat[i,"ntot"]
  if(!is.na(cnvdat[i,"frac2_A"])){
    copyNumber.ssm[ssms, "subclonal"] <- 1
    # copyNumber.ssm[ssms, "tumor_cn"] <- round((cnvdat[i,"nMaj1_A"] + cnvdat[i,"nMin1_A"])*cnvdat[i,"frac1_A"] + (cnvdat[i,"nMaj2_A"] + cnvdat[i,"nMin2_A"])*cnvdat[i,"frac2_A"])
  } else {
    copyNumber.ssm[ssms, "tumor_cn"] <- cnvdat[i,"nMaj1_A"] + cnvdat[i,"nMin1_A"]
  }
}
copyNumber.ssm[intersect(rownames(ssmdat)[id.xy],cnv_ssm[[length(cnv_ssm)]]),"tumor_cn"] <- copyNumber.ssm[intersect(rownames(ssmdat)[id.xy],cnv_ssm[[length(cnv_ssm)]]),"tumor_cn"]/2
copyNumber.ssm[,"coef1"] <- cellularity0*copyNumber.ssm[,"tumor_cn"] + (1 - cellularity0)*copyNumber.ssm[,"normal_cn"]
tmp = which(as.numeric(tumor_stat[,"FA"])*copyNumber.ssm[,"coef1"] >= 1.5)
copyNumber.ssm[tmp, "coef1"] <- copyNumber.ssm[tmp, "coef1"]/2

## Detect possible false positves (FPs) using three conditions
# Condition1: 32 <= BQ <= 39, where BQ is average base quality for reads supporting alleles
id.BQ <- union(which(as.numeric(tumor_stat[,"BQ"]) < 32), which(as.numeric(tumor_stat[,"BQ"]) > 39))
# Condition2: VAF of matched normal samples == 0
id.posNormal <- which(as.numeric(normal_stat[,"FA"])> 0)
# Condition3: not present in dbSNP, i.e. not having reference SNP ID
id.DB <- which(ssmdat[,3] != ".")
id.FP.est <- union(union(id.BQ, id.posNormal),id.DB)

if(length(id.FP.est) > 0) {
  h <- list()
  interval <- 0.05
  xlim <- ceiling(max(as.numeric(tumor_stat[,"FA"])*copyNumber.ssm[, "coef1"]))
  h[["TP"]] <- hist(as.numeric(tumor_stat[-id.FP.est,"FA"])*copyNumber.ssm[-id.FP.est, "coef1"], breaks = seq(0,xlim,interval), plot = FALSE)
  h[["FP"]] <- hist(as.numeric(tumor_stat[id.FP.est,"FA"])*copyNumber.ssm[id.FP.est, "coef1"], breaks = seq(0,xlim,interval), plot = FALSE)
  
  ratio <- (h[["FP"]]$counts) / (h[["TP"]]$counts)
  ratio[which(is.na(ratio))] = 0
  ratio[which(ratio == Inf)] = 10
  for(i in 2:length(ratio)){
    if(ratio[i] < 1)
      break
  }
  if(ratio[i-1]>=1){
    max.FP <- interval*(i-1)
  } else
    max.FP <- 0
  id1.to.remove <- which(as.numeric(tumor_stat[,"FA"])*copyNumber.ssm[,"coef1"] <= max.FP)
  
  # Estimate the size of false postive SSMs
  num.FP <- round((sum(as.numeric(tumor_stat[id.FP.est,"FA"])*copyNumber.ssm[id.FP.est,"coef1"] <= max.FP) + length(id.FP.est)) / 2)
} else {
  id1.to.remove <- NULL
  num.FP <- 0
}

## Select mutations satisfying total copy number == 1 in tumor cells 
## so that only one possible case for copy number of variant allele, i.e. C_var = C_tot = 1
id1 <- which(round(copyNumber.ssm[,"tumor_cn"], digits = 1) == 1)
if(length(id1) > 30) {

  set.seed(123456) 
  GS.data.binomial.0 <- subclone.dirichlet.gibbs(C=10+1, iter = 5000,
                                                 a = mutReads[id1], d = totalReads[id1], 
                                                 copyNumberCoef = copyNumber.ssm[id1,"coef1"])
  GS.data.binomial <- GS.data.binomial.0
  Gibbs.subclone.density.est(GS.data.binomial)
  cluster.assignment <- getClusterAssignment(GS.data.binomial, density.threshold = 0.01)
  occupied.clusters = sort(unique(cluster.assignment$most.likely.cluster))
  optima.0 = cluster.assignment$localOptima[occupied.clusters]
  num.clusters <- length(optima.0)

} else {
  num.clusters <- 8
  optima.0 <- NULL
}

if(length(which(optima.0 >= 0.1)) == 1 && sum(copyNumber.ssm[,"subclonal"]) == 0) {
  # If satisfying these two conditions, we claim there is only one subpopulation in the sample, 
  # and hence we can answer all the subchallenges.
  # where we ignore subpopulation with cellular proportion smaller than 10% due to possible false positive effect 
  
  ## Difference 1 from 'final1'
  optima = max(optima.0)
  cellularity = optima
  no.clusters = 1
  assignments <- rep(1, nrow(ssmdat))
  CCM <- matrix(1, nrow = no.muts, ncol = no.muts)
  nodes <- matrix(0,nrow = 1, ncol = 1)
  
  write.table(assignments,"GISL_subchallenge2A.txt",row.names=F,col.names=F,quote=F,sep="\t")
  write.table(nodes,"GISL_subchallenge3A.txt",row.names=T,col.names=F,quote=F,sep="\t")
  
  # Deal with possible memory problem in case of large matrix
  if(no.muts > 5000) {
    times <- (no.muts^2) %/% (5000^2)
    interval <- no.muts %/% times
    for(i in 1:times) {
      fwrite(as.data.frame(matrix(1, nrow = interval, ncol = no.muts)), "GISL_subchallenge2B.txt", append = TRUE, col.names = F, quote = F, sep = "\t")
      fwrite(as.data.frame(matrix(0, nrow = interval, ncol = no.muts)), "GISL_subchallenge3B.txt", append = TRUE, col.names = F, quote = F, sep = "\t")
      gc(verbose = F)
    }
    fwrite(as.data.frame(matrix(1, nrow = no.muts %% times, ncol = no.muts)), "GISL_subchallenge2B.txt", append = TRUE, col.names = F, quote = F, sep = "\t")
    fwrite(as.data.frame(matrix(0, nrow = no.muts %% times, ncol = no.muts)), "GISL_subchallenge3B.txt", append = TRUE, col.names = F, quote = F, sep = "\t")
  } else {
    fwrite(as.data.frame(matrix(1,nrow = no.muts,ncol = no.muts)), "GISL_subchallenge2B.txt", col.names = F, quote = F, sep = "\t")
    gc(verbose = FALSE)
    
    fwrite(as.data.frame(matrix(0,nrow = no.muts,ncol = no.muts)), "GISL_subchallenge3B.txt", col.names = F, quote = F, sep = "\t")
    gc(verbose = FALSE)
  }
  
} else {
  
  set.seed(123456) 
  ## Remove the possible FPs and use the estimated limit (flexibility allowed by adding pseudocount == 2)
  GS.data.binomial.1 <- subclone.dirichlet.gibbs(C = (num.clusters+2)+1, iter = 10000, 
                                                 a = mutReads[setdiff(1:no.muts,id1.to.remove)], d = totalReads[setdiff(1:nrow(ssmdat),id1.to.remove)], 
                                                 copyNumberCoef = copyNumber.ssm[setdiff(1:nrow(ssmdat),id1.to.remove),"coef1"])
  GS.data.binomial <- GS.data.binomial.1
  Gibbs.subclone.density.est(GS.data.binomial)
  cluster.assignment <- getClusterAssignment(GS.data.binomial, density.threshold = 0.01)
  occupied.clusters = sort(unique(cluster.assignment$most.likely.cluster), decreasing = TRUE)
  optima.1 = sort(cluster.assignment$localOptima[occupied.clusters], decreasing = TRUE)
  ## Difference 1 from 'final1'
  cellularity <- max(max(optima.1), max(optima.0))
  no.clusters <- length(optima.1)
  optima <- c(cellularity, optima.1[2:no.clusters])
  
  
  ########## Subchallenge 3A: Constructing phylogenetic trees ##########
  
  nodes <- matrix(nrow = length(optima), ncol = 5)
  # Assuming each node has at most two children
  colnames(nodes) <- c("proportion", "proportion_left", "parent","lchild","rchild")
  nodes[,"proportion"] <- optima
  nodes[,"proportion_left"] <- optima
  nodes[1, "parent"] <- 0
  
  unassigned <- which(is.na(nodes[,"parent"]))
  while(length(unassigned)) {
    add.to <- setdiff(which(is.na(nodes[,"rchild"])),unassigned)
    for(i in add.to) {
      if(nodes[unassigned[1],"proportion"] <= nodes[i,"proportion_left"]){
        nodes[unassigned[1], "parent"] <- i
        nodes[i,"proportion_left"] <- nodes[i,"proportion_left"]-nodes[unassigned[1],"proportion"]
        if(is.na(nodes[i,"lchild"]))
          nodes[i,"lchild"] <- unassigned[1]
        else
          nodes[i,"rchild"] <- unassigned[1]
        break
      }
    }
    unassigned <- which(is.na(nodes[,"parent"]))
  }
  
  write.table(nodes[,"parent"],"GISL_subchallenge3A.txt",row.names=T,col.names=F,quote=F,sep="\t")
  
  
  ########## Subchallenge 2 & 3B: Determine mutation assignments ##########

  ## Difference 2 from 'final1'
  # Determine mutation assignments using directly the results of Dirichlet process
  assignments <- rep(no.clusters, no.muts)
  assignments[setdiff(1:no.muts,id1.to.remove)] <- match(cluster.assignment$most.likely.cluster,occupied.clusters)
  
  write.table(assignments,"GISL_subchallenge2A.txt",row.names=F,col.names=F,quote=F,sep="\t")
  
  ## Co-clustering matrix for SC_2B
  ## Difference 3 from 'final1'
  # Compute probability of co-clustering matrix by binomial distribution
  new_coef1 <- cellularity*copyNumber.ssm[,"tumor_cn"] + (1 - cellularity)*copyNumber.ssm[,"normal_cn"]
  probs <- matrix(0, nrow = no.muts, ncol = no.clusters)
  for(i in 1:no.muts) {
    p <- optima/new_coef1[i]
    tmp <- sapply(p, function(x) dbinom(x = mutReads[i], size = totalReads[i], prob = x))
    probs[i,] <- tmp/sum(tmp)
  }
  
  # Co-clustering matrix for SC_2B
  CCM <- matrix(0, nrow = no.muts, ncol = no.muts)
  diag(CCM) <- 1
  for(i in 1:(no.muts-1)){
    CCM[i, (i+1):no.muts] <- sapply((i+1):no.muts, function(x) min(1,sum(probs[i,]*probs[x,])))
    CCM[(i+1):no.muts, i] <- CCM[i, (i+1):no.muts]
  }
  CCM <- round(CCM, digits = 4)
  
  fwrite(as.data.frame(CCM), "GISL_subchallenge2B.txt", col.names = F, quote = F, sep = "\t")
  
  # Find ancestors for each cluster
  ancestors <- vector("list", no.clusters)
  for(i in 2:no.clusters) {
    parent <- nodes[i,"parent"]
    # parent node's ancestors are also its ancestors
    ancestors[[i]] <- c(ancestors[[parent]], parent)
  }
  descendants <- sapply(1:no.clusters, function(i) which(sapply(ancestors, function(x) (i %in% x))))
  
  ADM <- matrix(0, nrow = no.muts, ncol = no.muts)
  anc.clusters <- unique(unlist(ancestors))
  dec.clusters <- unique(unlist(descendants))
  for(i in 1:(no.muts-1)) {
    for(j in (i+1):no.muts) {
      p1 = 0; p2 = 0
      for(c in anc.clusters) {
        p1 <- p1 + probs[i, c]*sum(probs[j, descendants[[c]]])
      }
      for(d in dec.clusters) {
        p2 <- p2 + probs[i, d]*sum(probs[j, ancestors[[d]]])
      }
      ADM[i,j] <- min(1,round(p1, digits = 4))
      ADM[j,i] <- max(min(round(p2, digits = 4),1-ADM[i,j]-CCM[i,j]),0)
    }
  }
  
  
  rm(CCM); gc(verbose = FALSE)
  
  fwrite(as.data.frame(ADM), "GISL_subchallenge3B.txt", col.names = F, quote = F, sep = "\t")
  rm(ADM); gc(verbose = FALSE)
  
}

