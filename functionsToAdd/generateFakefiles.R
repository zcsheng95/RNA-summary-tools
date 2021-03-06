set.seed(1234)
ngene <- 2000
reads <- 55
depth <- 100000
samples <- 5
rootdir <- "/home/zs68@dhe.duke.edu/Documents/dcibioinformatics/dcibioinformatics/inst/extdata"


### Create Sample folders
# system(paste0('cd ',rootdir, '\ rm *'))
for(i in 1:samples){
  folder <- file.path(rootdir,paste0("sample",i),"star")
  if(!dir.exists(folder))
  dir.create(folder)
}
#### Create fake files ReadsPerGene.out.tab
rnames <- c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous", paste0("gene",1:ngene,".1"))

for(i in 1:samples){
  filenames <- "ReadsPerGene.out.tab"
  outfile <- file.path(rootdir,paste0("sample",i),"star",filenames)
  ## first 4 rows
  p <- runif(n = 4, min =0,1)
  ## rest 2000 genes
  x <- rmultinom(3, depth, c(0.4*(p/sum(p)),rep.int(0.6/ngene, ngene)))
  row.names(x) <- rnames
  write.table(x, file = outfile,quote = F, col.names = F, sep = '\t')
}


#### Create fake files sample.metrics
pirnames <- c("PF_BASES	PF_ALIGNED_BASES	RIBOSOMAL_BASES	CODING_BASES	UTR_BASES	INTRONIC_BASES	INTERGENIC_BASES	IGNORED_READS	CORRECT_STRAND_READS	INCORRECT_STRAND_READS	NUM_R1_TRANSCRIPT_STRAND_READS	NUM_R2_TRANSCRIPT_STRAND_READS	NUM_UNEXPLAINED_READS	PCT_R1_TRANSCRIPT_STRAND_READS	PCT_R2_TRANSCRIPT_STRAND_READS	PCT_RIBOSOMAL_BASES	PCT_CODING_BASES	PCT_UTR_BASES	PCT_INTRONIC_BASES	PCT_INTERGENIC_BASES	PCT_MRNA_BASES	PCT_USABLE_BASES	PCT_CORRECT_STRAND_READS	MEDIAN_CV_COVERAGE	MEDIAN_5PRIME_BIAS	MEDIAN_3PRIME_BIAS	MEDIAN_5PRIME_TO_3PRIME_BIAS	SAMPLE	LIBRARY	READ_GROUP")
pirnames <- strsplit(pirnames,split = "\t")[[1]]
for (i in 1:samples) {
  filenames <- paste0("sample",i,".metrics")
  outfile <- file.path(rootdir, paste0("sample",i),filenames)
  file.create(outfile)
  #### Writing Lines to files
  sink(outfile)
  for(i in 1:4){cat("## Line",i,"\n")}
  cat("\n")
  cat("## METRICS CLASS        picard.analysis.RnaSeqMetrics\n")
  total <- depth * reads
  probs <- vector(mode = "numeric", length = 28)
  x <- c(0.05, 0.65, 0.20, 0.07, 0.015) + runif(5,-0.015,0.015) 
  probs[1:5] <- x/sum(x)
  aligned <- ceiling(runif(n = 1, min = 0.97 * total, max = 0.99 * total)) 
  bases <- c(total, aligned,floor(probs * aligned))
  out <- matrix(bases, nrow = 1, ncol = 30)
  cat(pirnames, sep = '\t')
  cat('\n')
  cat(out, sep = '\t')
  cat('\n\n')
  cat("## HISTOGRAM    java.lang.Integer\n## Some Data")
  sink()
}

