
source("lib.R")
#source("class.R")

#VCscript <- function(pc_csv_file){

#pc_csv_file <- "vc_pc_name.csv"
pc_csv_file <- "remark_vc.csv"

### Read the genome data
cat("\nReading the viral reference data .....\n")
pc_csv <- read.csv(pc_csv_file, header=TRUE)
#cat(sprintf("Size of genome data (from %s): %d\n",pc_csv_file,dim(pc_csv)[1]))
contig <- pc_csv$contig_id
vc     <- pc_csv$VC
pc     <- pc_csv$pc_id
rm(pc_csv)

# Sorting
contig       <- contig[order(contig)]
vc           <- vc[order(contig)]
pc           <- pc[order(contig)]

uniq_contig  <- unique(contig)  # Identify unique contig
uniq_vc      <- unique(vc)
uniq_pc      <- unique(pc)
uniq_pc      <- uniq_pc[order(uniq_pc)]
uniq_pc      <- uniq_pc[uniq_pc != ""]  # remove singletons

num_contig     <- length(uniq_contig)
num_vc         <- length(uniq_vc)
num_pc         <- length(uniq_pc)
num_singleton  <- length(which(pc == ""))

cat(sprintf("  Total number of contigs: %d\n",num_contig))
cat(sprintf("  Total number of protein clusters: %d\n",num_pc))
cat(sprintf("  Total number of virus clusters: %d\n",num_vc))
cat(sprintf("  Total number of singletons: %d\n",num_singleton))

### Build the matrix
cat("\nBuilding PC matrix .....\n")
ptm <- proc.time()
pc_contig <- matrix(0,num_pc,num_contig)
vc_contig <- matrix(0,num_vc,num_contig)
num_contig_singletron <- matrix(0,num_contig)
for(i in 1:num_contig){
  tmp <- which(contig %in% uniq_contig[i])
  num_contig_singletron[i] <- length(which(pc[tmp] == ""))
  pc_contig[match(unique(pc[tmp]),uniq_pc),i] <- 1
  vc_contig[match(unique(vc[tmp]),uniq_vc),i] <- 1
}
num_contig_pc  <- colSums(pc_contig)
num_contig_pcs <- num_contig_pc + num_contig_singletron
common_pc      <- t(pc_contig) %*% pc_contig
print_run_time(ptm)

contig_pairs     <- combn(num_contig,2)
num_contig_pairs <- dim(contig_pairs)[2]
prcnt_pc  <- matrix(0,num_contig_pairs,2)
ptm <- proc.time()
for (k in 1:num_contig_pairs){
  i <- contig_pairs[1,k]
  j <- contig_pairs[2,k]
  prcnt_pc[k,1] <- common_pc[i,j]/num_contig_pcs[i]
  prcnt_pc[k,2] <- common_pc[i,j]/num_contig_pcs[j]
}
prcnt_pc <- prcnt_pc * 100.0
print_run_time(ptm)


### intra & inter VC percentages
contig_vc_id <- which(vc_contig==1,arr.ind=T)[,1]
#try(
  if(length(contig_vc_id) > num_contig){
    stop("STOP!!!! Number of VC belongings is larger than number of contigs")
  }
#)

#prcnt_vc <- matrix(0,num_contig,2)
prcnt_vc <- matrix(0,num_contig,num_vc)
contig_vc <- matrix(0,num_contig)
common_pc <- common_pc - diag(num_contig_pc)
num_shared_pc <- colSums(common_pc)
ptm <- proc.time()
for (i in 1:num_contig){
  contig_vc[i] <- paste(uniq_vc[contig_vc_id[i]])
  if(num_shared_pc[i] == 0) next
  #prcnt_vc[i,1] <- sum(vc_contig[contig_vc_id[i],]*common_pc[i,])/num_shared_pc[i]
  tmp <- vc_contig[contig_vc_id[i],]*common_pc[i,]
  vec_nz <- which(tmp!=0)
  vc_id <- contig_vc_id[vec_nz]
  for(j in 1:length(vc_id)){
    prcnt_vc[i,vc_id[j]] <- prcnt_vc[i,vc_id[j]] + tmp[vec_nz[j]]/num_shared_pc[i]
  }
}
neg_vc_contig <- matrix(1,(dim(vc_contig)[1]),(dim(vc_contig)[2])) - vc_contig
for (i in 1:num_contig){
  if(num_shared_pc[i] == 0) next
  #prcnt_vc[i,2] <- sum(neg_vc_contig[contig_vc_id[i],]*common_pc[i,])/num_shared_pc[i]
  tmp <- neg_vc_contig[contig_vc_id[i],]*common_pc[i,]
  vec_nz <- which(tmp!=0)
  vc_id <- contig_vc_id[vec_nz]
  for(j in 1:length(vc_id)){
    prcnt_vc[i,vc_id[j]] <- prcnt_vc[i,vc_id[j]] + tmp[vec_nz[j]]/num_shared_pc[i]
  }
}
prcnt_vc <- prcnt_vc * 100.0
print_run_time(ptm)
#rm(common_pc)

###
### Script 2 and 3
###
pc_vc <- matrix(0,num_pc,num_vc)
for(i in 1:num_vc){
  tmp <- pc[which(vc %in% uniq_vc[i])]
  pc_vc[match(unique(tmp),uniq_pc),i] <- 1
}

### Find hallmark genes 
vc_pairs <- combn(num_vc,2)
num_vc_pairs <- dim(vc_pairs)[2]
hallmark_pc <- matrix(0,num_vc_pairs,num_pc)
for (k in 1:num_vc_pairs){
  i <- vc_pairs[1,k]
  j <- vc_pairs[2,k]
  hallmark_pc[k,] <- pc_vc[,i]*pc_vc[,j]
}

### Find signature genes
signature_pc <- t(pc_vc)
for (k in 1:num_vc_pairs){
  i <- vc_pairs[1,k]
  j <- vc_pairs[2,k]
  signature_pc[i,] <- signature_pc[i,] - hallmark_pc[k,]
  signature_pc[j,] <- signature_pc[j,] - hallmark_pc[k,]
} 

###
### Make some output
###
prcnt_contig <- cbind(paste(contig_vc[contig_pairs[1,]]),paste(contig_vc[contig_pairs[2,]]),
		      paste(uniq_contig[contig_pairs[1,]]),paste(uniq_contig[contig_pairs[2,]]))
prcnt_pc[which(prcnt_pc==0)] <- ""
prcnt_csv <- cbind(prcnt_contig,prcnt_pc)
colnames(prcnt_csv) <- c("VC1","VC2","contig1","contig2","% of contig2","% of contig1")
write.table(prcnt_csv,file="percentage.csv",quote=F,sep=',',row.names=F)

prcnt_contig <- cbind(paste(uniq_contig))
prcnt_vc[which(prcnt_vc==0)] <- ""
prcnt_csv <- cbind(prcnt_contig,contig_vc,prcnt_vc)
colnames(prcnt_csv) <- c("contig","VC",rbind(paste(uniq_vc)))
write.table(prcnt_csv,file="percentage_vc.csv",quote=F,sep=',',row.names=F)

hallmark_count <- rowSums(hallmark_pc)
hallmark_vc <- cbind(paste(uniq_vc[vc_pairs[1,]]),paste(uniq_vc[vc_pairs[2,]]))
hallmark_pc_id <- apply(hallmark_pc,1,function(x) paste(uniq_pc[which(x==1)],collapse=','))
#hallmark_pc_id <- matrix(0,num_vc_pairs)
#for (k in 1:num_vc_pairs){
#  hallmark_pc_id[k] <- paste(uniq_pc[which(hallmark_pc[k,]==1)],collapse=',')
#}
hallmark_tsv <- cbind(hallmark_vc,hallmark_count,hallmark_pc_id)
colnames(hallmark_tsv) <- c("vc1","vc2","count","pcs")
write.table(hallmark_tsv,file="hallmark.tsv",quote=F,sep='\t',row.names=F)

signature_count <- rowSums(signature_pc>0)
signature_vc <- cbind(paste(uniq_vc))
signature_pc_id <- apply(signature_pc,1,function(x) paste(uniq_pc[which(x==1)],collapse=','))
#signature_pc_id <- matrix(0,num_vc)
#for (k in 1:num_vc){
#  signature_pc_id[k] <- paste(uniq_pc[which(signature_pc[k,]==1)],collapse=',')
#}
signature_tsv <- cbind(signature_vc,signature_count,signature_pc_id)
colnames(signature_tsv) <- c("vc","count","pcs")
write.table(signature_tsv,file="signature.tsv",quote=F,sep='\t',row.names=F)

#}
