" The script 'pcprofiles' computes the similarity network on the contigs with 
  the significant value above than a given threshold "

source("lib.R")
source("class.R")

#options(echo=TRUE) # if you want see commands in output file
#args <- commandArgs(trailingOnly=TRUE)
#print(args)

pcprofiles <- function(pc_csv_file, thresh=1.0, sig_csv_file="sig.csv", doNeg=FALSE){

#if(length(args) < 1){
#  cat(sprintf("Usage: Rscript current_script  'pc_csv' '['threshold = 1.0'] ['sig_csv_file']\n"))
#  stop("At least one argument must be supplied.\n", call.=FALSE)
#}else if(length(args) == 1){ 
#  thresh  <- 1.0
#  sig_csv_file <- "sig.csv"
#}else if(length(args) == 2){ 
#  if(!is.numeric(args[1])){
#    stop("The second argument must be a number.\n", call.=FALSE)
#  }
#  sig_csv_file <- "sig.csv"
#}else{
#  thresh  <- as.double(args[2])
#  sig_csv_file <- args[3]
#}
#  pc_csv_file  <- args[1]

vr <- ViralRef()

### Read the genome data
cat("\nReading the viral reference data .....\n")
pc_csv <- read.csv(pc_csv_file, header=TRUE)
#cat(sprintf("Size of genome data (from %s): %d\n",pc_csv_file,dim(pc_csv)[1]))
contig <- pc_csv$contig_id
pc     <- pc_csv$pc_id
rm(pc_csv)

# Sorting
pc           <- pc[order(contig)]
contig       <- contig[order(contig)]
uniq_contig  <- unique(contig)	# Identify unique contig
uniq_pc      <- unique(pc)
uniq_pc      <- uniq_pc[order(uniq_pc, decreasing=TRUE)]
uniq_pc      <- uniq_pc[uniq_pc != ""]	# remove singletons

vr@num_contig     <- length(uniq_contig)
vr@num_pc         <- length(uniq_pc)
vr@num_singleton  <- length(which(pc == ""))
logT              <- log10(choose(vr@num_contig,2))

cat(sprintf("  Total number of contigs: %d\n",vr@num_contig))
cat(sprintf("  Total number of protein clusters: %d\n",vr@num_pc))
cat(sprintf("  Total number of singletons: %d\n",vr@num_singleton))

### Build the matrix
cat("\nBuilding PC matrix .....\n")
ptm <- proc.time()
vr_mat <- matrix(0,vr@num_pc,vr@num_contig)
num_contig_singletron <- matrix(0,vr@num_contig)
for(i in 1:vr@num_contig){
  tmp <- pc[which(contig %in% uniq_contig[i])]
  num_contig_singletron[i] <- length(which(tmp == ""))
  vr_mat[match(unique(tmp),uniq_pc),i] <- 1
}
vr@mat <- vr_mat
rm(vr_mat)
vr@pc     <- uniq_pc
vr@contig <- uniq_contig
rm(pc,contig,uniq_contig,uniq_pc)

num_contig_pc  <- colSums(vr@mat)
num_contig_pcs <- num_contig_pc + num_contig_singletron
common_pc      <- t(vr@mat) %*% vr@mat 
print_mem_usage(object.size(common_pc))

print_run_time(ptm)

contig_pairs     <- combn(vr@num_contig,2)
num_contig_pairs <- dim(contig_pairs)[2]
cat(sprintf("\nCalculating the Sig. values for %d pairs of genomes\n",num_contig_pairs))
n <- vr@num_pc + vr@num_singleton
vr_sig  <- matrix(0,num_contig_pairs,5)
ptm <- proc.time()
for (k in 1:num_contig_pairs){
  i <- contig_pairs[1,k]
  j <- contig_pairs[2,k]
  a <- num_contig_pcs[i]
  b <- num_contig_pcs[j]
  c <- common_pc[i,j]
  pval <- dhyper(c,a,(n-b),b)
  sig <- -log10(pval)-logT
   
  vr_sig[k,] <- c(a,b,c,pval,sig)
}
print_mem_usage(object.size(vr_sig))
print_run_time(ptm)

if(doNeg){

  ### Kick out the negative group
  vr@sig     <- vr_sig[,5]
  sig_vec    <- which(vr_sig[,5] < thresh & vr_sig[,5] > 0)
  num_sig    <- length(sig_vec)
  sig_pairs  <- t(contig_pairs[,sig_vec])
  vr_sig     <- vr_sig[sig_vec,]
  cat(sprintf("\nTotal number of pairs with a Sig. value less than %.1f: %d\n",thresh,num_sig))
  sig_csv <- cbind(sig_pairs,vr_sig[,c(1:3,5)])
  colnames(sig_csv) <- c("contig1_id","contig2_id","a","b","c","sig")
  write.csv(sig_csv,file="perm0.csv",row.names=FALSE)

  cat("Return ViralRef object for the permtuation test\n")
  return(vr)

}else{

  ### Output similar pairs 
  ptm <- proc.time()
  sig_vec    <- which(vr_sig[,5] > thresh)
  num_sig    <- length(sig_vec)
  sig_pairs  <- t(contig_pairs[,sig_vec])
  sig_contig <- cbind(paste(vr@contig[sig_pairs[,1]]),paste(vr@contig[sig_pairs[,2]]))
  vr_sig     <- vr_sig[sig_vec,]
  # set the Sig. maximum to 300.0
  vr_sig[which(vr_sig[,5] > 300),5] <- 300
  #cat(sprintf("Total number of pairs with a Sig. value larger than %.1f: %d\n",thresh,num_sig))
  cat(sprintf("\nKicking out %d pairs with a Sig. value larger than %.1f\n",num_sig,thresh))
  cat(sprintf("  Output is %s\n",sig_csv_file))
  sig_csv <- cbind(sig_pairs,vr_sig)
  sig_csv <- cbind(sig_contig,vr_sig)
  colnames(sig_csv) <- c("contig1_id","contig2_id","a","b","c","pval","sig")
  write.csv(sig_csv,file=sig_csv_file,row.names=FALSE)
  print_run_time(ptm)

}

}

