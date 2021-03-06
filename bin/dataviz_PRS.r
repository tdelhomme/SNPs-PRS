library(ROCR)
library(ggpubr)
library(pROC)

labels_files = list.files(path = ".", pattern = "labels")
cancertypes = unlist(lapply(labels_files, function(x) unlist(strsplit(x,"_translate_labels"))[1]))

for(ct in cancertypes){
  pdf(paste("ROC_",ct,".pdf",sep=""),12,6)
  for(oncotype in c("MUTATED", "NONMUTATED")){
    if(oncotype == "MUTATED") { ff = "1" ; tt = "2" }
    if(oncotype == "NONMUTATED") { tt = "1" ; ff = "2" }
    load(paste(ct, "_translate_labels.Rdata",sep=""))
    load(paste(ct, "_translate_PRS.Rdata",sep=""))
    genes =  rownames(get(ls(pattern = paste(ct,"_translate__PRS",sep=""))))
    dat = data.frame(gene = genes, cancertype = ct, AUC = NA, PREC = NA, REC=NA)
    ### plot precision/recall and ROC ###
    print(paste(date(), "  INFO: plot ROC per gene for cancer type: ", ct, sep=""))
    colors <- rainbow(length(genes))
    par(mfrow=c(1,2))
    plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),ylab="Precision",xlab="Recall",bty='n',main=oncotype)
    for(gene in genes){
      score = get(paste(ct,"_translate__PRS",sep=""))[gene,]
      test.labelsold = get(paste(ct,"_translate__testlabel",sep=""))[,gene]
      # we can't use FALSE/TRUE as labels when negative scores because it considers s<t as false
      # and in our case we can have true <t (because we labelled with t/f)
      test.labels=rep(tt, length(test.labelsold))
      test.labels[which(test.labelsold == FALSE)] = ff
      pred <- prediction(score, test.labels)
      perf <- performance(pred, "prec", "rec")
      roc.x <- unlist(perf@x.values)
      roc.y <- unlist(perf@y.values)
      lines(roc.y ~ roc.x, col = colors[which(genes == gene)], lwd = 2)
      dat[match(gene, genes),"PREC"] = ifelse(length(which(roc.x > 0.1))==0,
                                              0, max(roc.y[which(roc.x > 0.1)], na.rm=T))
      dat[match(gene, genes),"REC"] = ifelse(length(which(roc.y >= 0.75))==0,
                                             0, max(roc.x[which(roc.y >= 0.75)], na.rm=T))
    }
    plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),ylab="Sensitivity",xlab="Specificity",bty='n',main=oncotype)
    for(gene in genes){
      score = get(paste(ct,"_translate__PRS",sep=""))[gene,]
      test.labelsold = get(paste(ct,"_translate__testlabel",sep=""))[,gene]
      # we can't use FALSE/TRUE as labels when negative scores because it considers s<t as false
      # and in our case we can have true <t (because we labelled with t/f)
      test.labels=rep(tt, length(test.labelsold))
      test.labels[which(test.labelsold == FALSE)] = ff
      pred <- prediction(score, test.labels)
      perf <- performance(pred, "sens", "spec")
      roc.x <- unlist(perf@x.values)
      roc.y <- unlist(perf@y.values)
      lines(roc.y ~ roc.x, col = colors[which(genes == gene)], lwd = 2)
      dat[match(gene, genes),"AUC"] = round(auc(test.labels, score,
                                                levels = c(ff,tt), direction = "<"),3)
    }
    legend("bottomleft", legend = genes, col = colors, pch = 19, horiz = F, cex=0.8)
    if(match(ct, cancertypes) == 1) {
      assign(paste("all_dat_", oncotype, sep=""),dat) } else {
        assign(paste("all_dat_", oncotype, sep=""), rbind(get(paste("all_dat_", oncotype, sep="")), dat)) }
  }
  dev.off()
}

# use the correct cancer type ID
for(oncotype in c("MUTATED", "NONMUTATED")){
  all_dat = get(paste("all_dat_", oncotype, sep=""))
  all_dat[grep("Brain", all_dat$cancertype), "tissue"] = "LGG"
  all_dat[grep("Breast", all_dat$cancertype), "tissue"] = "BRCA"
  all_dat[grep("Kidney", all_dat$cancertype), "tissue"] = "KIRCP"
  all_dat[grep("Lung_Aden", all_dat$cancertype), "tissue"] = "LUAD"
  all_dat[grep("Skin", all_dat$cancertype), "tissue"] = "SKCM"
  all_dat[grep("Lung_Squam", all_dat$cancertype), "tissue"] = "LUSC"
  
  ### plot boxplots ###
  p1 <- ggboxplot(all_dat, x = "tissue",
                  y = "AUC",
                  color = "tissue", palette = "jco",
                  add = "jitter", outlier.shape = NA,           # Add jittered points (remove outliers otherwise duplicated)
                  #add.params = list(size = 0.5, jitter = 0.2),  # Point size and the amount of jittering
                  label = "gene",                               # column containing point labels
                  label.select = list(criteria = "`y` >= 0.75"),  # Select some labels to display
                  font.label = list(size = 9, face = "italic"), # label font
                  repel = TRUE,                                 # Avoid label text overplotting
                  legend = "NONE", ylim = c(0,1)
  )
  
  p2 <- ggboxplot(all_dat, x = "tissue",
                  y = "PREC", ylab = "Max Precision(Recall>0.1)",
                  color = "tissue", palette = "jco",
                  add = "jitter",  outlier.shape = NA,          # Add jittered points (remove outliers otherwise duplicated)
                  #add.params = list(size = 0.1, jitter = 0.2),  # Point size and the amount of jittering
                  label = "gene",                               # column containing point labels
                  label.select = list(criteria = "`y` >= 0.75"),  # Select some labels to display
                  font.label = list(size = 9, face = "italic"), # label font
                  repel = TRUE,                                 # Avoid label text overplotting
                  legend = "NONE", ylim = c(0,1)
  ) + geom_hline(yintercept = 0.75, linetype = 2)
  
  p3 <- ggboxplot(all_dat, x = "tissue",
                  y = "REC", ylab = "Max Recall(Precision>0.75)",
                  color = "tissue", palette = "jco",
                  add = "jitter",  outlier.shape = NA,          # Add jittered points (remove outliers otherwise duplicated)
                  #add.params = list(size = 0.1, jitter = 0.2),  # Point size and the amount of jittering
                  label = "gene",                               # column containing point labels
                  label.select = list(criteria = "`y` > 0.1"),  # Select some labels to display
                  font.label = list(size = 9, face = "italic"), # label font
                  repel = TRUE,                                 # Avoid label text overplotting
                  legend = "NONE", ylim = c(0,1)
  ) + geom_hline(yintercept = 0.1, linetype = 2)
  
  p4 <- ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3, nrow = 1)
  ggsave(paste("boxplot_auc_prec_",oncotype,".pdf",sep=""), plot = p4, width = 15, height = 4)
}
