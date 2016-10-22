print_run_time <- function(ptm){
  rtm <- proc.time() - ptm
  cat(sprintf("  Elapsed time: %.1f secs; User time: %.1f secs\n",rtm[3],rtm[1]))
}

print_mem_usage <- function(obj_size){
  cat("  Memory used:",format(obj_size, units="Kb"),"\n")
}

getSig <- function(a,b,c){
  n    <- 23021
  logT <- log10(choose(2008,2))
  pval <- dhyper(c,a,(n-b),b)
  return (-log10(pval)-logT)
}
