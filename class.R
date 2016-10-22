
ViralRef <- setClass(
  "ViralRef",
  slots = c(
	    mat="matrix",sig="vector",
	    num_pc="integer",num_contig="integer",num_singleton="integer",
	    contig="factor",pc="factor"	## should contain only characters
	   )
  )

