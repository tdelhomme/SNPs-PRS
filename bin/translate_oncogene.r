args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}

argsL <- as.list(do.call("cbind", parseArgs(args))[c(F,T)])
names(argsL) <- as.list(do.call("cbind", parseArgs(args))[c(T,F)])
args <- argsL;rm(argsL)

if(! is.null(args$help)) {
  cat("
      Mandatory arguments:
      --SNPs                     - input txt file containing SNPs in rows, samples in columns
      --mut_table                - MC3 mutation table
      --samples_table            - Input txt file ID and sample_id columns
      Optional arguments:
      --nb_cpu                   - number of cpu to use
      --no_cancer_type           - T/F indicating that we have merged all cancer types (default=FALSE, not merged)
      --help \n\n")
  q(save="no")
}

if(is.null(args$SNPs)) {stop("Option --SNPs should be provided")} else{SNPs=args$SNPs}
if(is.null(args$mut_table)) {stop("Option --mut_table should be provided")} else{mut_table=args$mut_table}
if(is.null(args$samples_table)) {stop("Option --samples_table should be provided")} else{samples_table=args$samples_table}
if(is.null(args$nb_cpu)) {nb_cpu=1} else{nb_cpu=as.numeric(args$nb_cpu)}
if(is.null(args$no_cancer_type)) {no_cancer_type=FALSE} else{no_cancer_type=TRUE}
print(paste(date(), "  INFO: number of CPUs: ", nb_cpu))

library(foreach)
library(parallel)

### loading the data ###
sm_table = read.table(samples_table, h = T, stringsAsFactors = F)
input_snps = read.table(SNPs, h = T, stringsAsFactors = F)
muts = read.table(mut_table, h = T, stringsAsFactors = F)

## keep SNPs for samples with the target tissue type we are analyzing ###
idsm = unlist(lapply(colnames(input_snps[, 4:ncol(input_snps)]), function(s) unlist(strsplit(s, "_"))[1]))
if(! no_cancer_type){
  input_snps = input_snps[, c(1:3, 3 + which(idsm %in% 
                                               sm_table[which(sm_table$cancer_type == unique(muts$cancer_type)), "ID"]))]
} # we don't need to do something when no_cancer_type because we keep only samples in our data (l.58)
top_genes = unique(muts$Hugo_Symbol)

print(paste(date(), "  INFO: starting to translate the table"))

# on table per each gene, to do not loose genotypes of not-mutated samples
cl <- parallel::makeCluster(nb_cpu)
doParallel::registerDoParallel(cl)
dat <- foreach(g = top_genes, .combine = cbind)  %dopar% {
  snps_sm = unlist(lapply(colnames(input_snps)[4:ncol(input_snps)], function(s){
    unlist(strsplit(s, "_"))[1]}))
  tcga_sm = unlist(lapply(snps_sm, function(s) sm_table[which(sm_table$ID == s), "sample_id"]))
  mut_samples = muts[which(muts$Hugo_Symbol == g), "SM"]
  res = data.frame(tcga_sm %in% mut_samples) # vector of true/false for being mutated
  colnames(res) = g
  res
}
parallel::stopCluster(cl)

r = t(as.matrix(input_snps[,4:ncol(input_snps)]))
colnames(r) = input_snps$rs
return_dat = cbind( r, dat)

print(paste(date(), "  INFO: writing the tables"))
write.table(return_dat,
            file = paste(tools::file_path_sans_ext(mut_table), "_translate.txt",sep=""),
            col.names = NA, row.names = T, quote = F, sep = "\t")

print(paste(date(), "  INFO: completed"))
