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
      --translate_table          - input txt file containing for each sample the SNPs data + gene data in columns
      --tag                      - string corresponding to cancer type
      --help \n\n")
  q(save="no")
}

if(is.null(args$translate_table)) {stop("Option --translate_table should be provided")} else{translate_table=args$translate_table}
if(is.null(args$tag)) {stop("Option --tag should be provided")} else{cancertype = as.character(args$tag)}

library(epitools)

###
### computation of the PRS per gene per sample ###
###

print(paste(date(), "  INFO: load the translated SNPs and convert them to numeric values"))
trans_snps = read.table(translate_table, h = T, stringsAsFactors = F)
snps = grep("SNP", colnames(trans_snps))
genes = colnames(trans_snps)[which(!grepl("SNP", colnames(trans_snps)))]
num_snps = apply( as.matrix(trans_snps[,snps]), 2, function(r)as.numeric(as.factor(r)) )
num_snps[which(num_snps == 3)] = 2 # consider AA and (AT,TT)
rownames(num_snps) = rownames(trans_snps)
train = sample(1:nrow(num_snps), nrow(num_snps)*0.65) # use 65% of samples for the beta estimation
notrain = setdiff( (1:nrow(num_snps)), train)

library(foreach)
cl <- parallel::makeCluster(10) # this should be an input parameter
doParallel::registerDoParallel(cl)

res = foreach( gene = genes, .combine = 'cbind' ) %dopar% {
  library(epitools)
  print(paste(date(), "  INFO: starting gene : ", gene, sep=""))
  # compute the betas
  betas = apply(num_snps[train,snps], 2, function(c){
    minor_geno = as.numeric(names(table(c))[which.min(table(c))]) # define the alt geno (=exposure) 
    # we need ordered vector because of the dat(table(...)) part, i.e. first line should be unexposed
    c2 = rep("1REF", length(c))
    c2[which(c == minor_geno)] = "2ALT"
    dat = table(data.frame(c2, trans_snps[train, gene]))
    if(nrow(dat) < 2){ OR=1 } else { # if we only see one genotype, then we can have any effect (log(1)=0)
      ORs = oddsratio.wald(dat)$measure
      OR = ORs[,"estimate"]["2ALT"]
    }
    return(log(OR))
  })
  # compute the PRS, i.e. a vector that contains the PRS value for each sample
  betas[which(is.infinite(betas))] = NA
  betas
}

PRS = apply(num_snps[notrain,snps], 1, function(r){
    prs_genes = r * res
    apply(prs_genes, 2, function(c) sum(as.numeric(c), na.rm=T))
})
rownames(PRS) = genes

assign(paste(cancertype,"__betas",sep=""), res)
assign(paste(cancertype,"__PRS",sep=""), PRS)
assign(paste(cancertype,"__testlabel",sep=""), trans_snps[notrain,genes])

save(list=paste(cancertype,"__betas", sep=""), file=paste(cancertype,"_betas.Rdata",sep=""))
save(list=paste(cancertype,"__PRS", sep=""), file=paste(cancertype,"_PRS.Rdata",sep=""))
save(list=paste(cancertype,"__testlabel", sep=""), file=paste(cancertype,"_labels.Rdata",sep=""))
