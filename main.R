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

# Use Battenberg data to estimate cellularity0 because we cannot directly use the given one for SC1 
# Note: it is not the final output

id2 <- which(cnvdat[,"pval"] < 0.05)
if(length(id2) > 0) {
  id2 <- id2[which.max(cnvdat[id2,"BAF"])]
  cellularity0 <- (2*cnvdat[id2,"BAF"] - 1) / (cnvdat[id2, "nMaj1_A"]*cnvdat[id2, "frac1_A"] + cnvdat[id2, "nMaj2_A"]*cnvdat[id2, "frac2_A"] - 1 - 
                                                 cnvdat[id2, "BAF"]*((cnvdat[id2, "nMaj1_A"]+cnvdat[id2,"nMin1_A"])*cnvdat[id2, "frac1_A"] + (cnvdat[id2, "nMaj2_A"] + cnvdat[id2, "nMin2_A"])*cnvdat[id2, "frac2_A"]-2))
} else {
  id1 <- which(cnvdat[,"pval"] == 1)
  id1 <- id1[which.max(cnvdat[id1,"BAF"])]
  cellularity0 <- (2*cnvdat[id1,"BAF"] - 1) / (cnvdat[id1, "nMaj1_A"] - 1 - cnvdat[id1, "BAF"]*(cnvdat[id1,"nMaj1_A"]+cnvdat[id1,"nMin1_A"]-2))
}

# Map mutations to CNV intervals according to their positions
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

  write.table(cellularity,"GISL_subchallenge1A.txt",row.names=F,col.names=F,quote=F,sep="\t")
  write.table(no.clusters,"GISL_subchallenge1B.txt",row.names=F,col.names=F,quote=F,sep="\t")
  write.table(cbind(1:(no.clusters+1),c(no.muts-num.FP, num.FP),c(round(optima, digits = 6),0)),"GISL_subchallenge1C.txt",row.names=F,col.names=F,quote=F,sep="\t")

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
  
  # Output results for SC1
  write.table(cellularity,"GISL_subchallenge1A.txt",row.names=F,col.names=F,quote=F,sep="\t")
  write.table(no.clusters,"GISL_subchallenge1B.txt",row.names=F,col.names=F,quote=F,sep="\t")
  
  ## Difference 2 from 'final1'
  # Determine mutation assignments using directly the results of Dirichlet process
  assignments <- rep(no.clusters, no.muts)
  assignments[setdiff(1:no.muts,id1.to.remove)] <- match(cluster.assignment$most.likely.cluster,occupied.clusters)
  
  tmp = c(table(assignments), num.FP)
  tmp[no.clusters] = max(0, tmp[no.clusters] - num.FP)
  write.table(cbind(1:(no.clusters+1),tmp,c(round(optima, digits = 6),0)),"GISL_subchallenge1C.txt",row.names=F,col.names=F,quote=F,sep="\t")
  
}

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

