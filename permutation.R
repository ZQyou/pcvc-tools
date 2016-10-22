" The script 'permutation' performs the permutation test over the pairs of genomes
  failed passing a given threshold of sig, as a negative control. "

source("lib.R")
source("class.R")
source("pcprofiles.R")

permutation <- function(pc_csv_file, thresh=1.0, num_random=1){

vr <- ViralRef()
vr <- pcprofiles(pc_csv_file,thresh,"",doNeg=TRUE)

### Permutation Test

#cat("Number of contig:",vr@num_contig,"\n")
contig_pairs <- combn(vr@num_contig,2)
perm_vec     <- which(vr@sig > 0 & vr@sig < thresh)
num_perm     <- length(perm_vec)
perm_pairs   <- t(contig_pairs[,perm_vec])
#cat(sprintf("Number of pairs with Sig. less than %.1f: %d\n",thresh,num_perm))

cat("\nStarting randomization test .....\n")
# Randomization 
for(r in 1:num_random){
  cat(sprintf("  # of Randomization: %d\n",r))
  ptm           <- proc.time()
  vr@mat        <- apply(vr@mat,2,sample)
  common_pc     <- t(vr@mat) %*% vr@mat
  num_contig_pc <- colSums(vr@mat)
  #print_mem_usage(object.size(common_pc))
  print_run_time(ptm)

  #cat("----- Calcuating Sig -----\n")
  perm_sig  <- matrix(0,num_perm)
  #perm_sig  <- matrix(0,num_perm,4)
  num_pc    <- vr@num_pc
  logT      <- log10(choose(vr@num_contig,2))

  #ptm <- proc.time()
  for(k in 1:num_perm){
    i <- perm_pairs[k,1]
    j <- perm_pairs[k,2] 
    a <- num_contig_pc[i]
    b <- num_contig_pc[j]
    c <- common_pc[i,j]
    pval <- dhyper(c,a,(num_pc-b),b)
    sig  <- -log10(pval)-logT

    perm_sig[k]  <- sig
    #perm_sig[k,]  <- c(a,b,c,sig)
  }
  #print_run_time(ptm)

  cat("  Found",length(which(perm_sig > thresh)),"pair(s) with a Sig. above the threshold\n")
  #length(which(perm_sig[,3] > 0))

  #perm_csv <- cbind(perm_pairs,perm_sig)
  #colnames(perm_csv) <- c("contig1_id","contig2_id","a","b","c","sig")
  #perm_csv_file <- paste(paste("perm",r,sep="_"),"csv",sep=".")
  #cat(sprintf("Writing the sig data for the permutated matrix into %s\n",perm_csv_file))
  #write.csv(perm_csv, file=perm_csv_file, row.names=FALSE)
}

}
