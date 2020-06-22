library(ROCR)
library(ggpubr)
library(pROC)

labels_files = list.files(path = ".", pattern = "PRS")
cancertypes = unlist(lapply(labels_files, function(x) unlist(strsplit(x,"_translate_labels"))[1]))

for(ct in cancertypes){
  load(paste(ct, "_translate_labels.Rdata",sep=""))
  load(paste(ct, "_translate_PRS.Rdata",sep=""))
  genes = unlist(lapply(ls(pattern = paste(ct,"_translate__testlabel",sep="")), 
                        function(x) unlist(strsplit(x,"testlabel"))[2]))
  dat = data.frame(gene = genes, cancertype = ct, AUC = NA, PREC = NA, REC=NA)
  ### plot precision/recall and ROC ###
  print(paste(date(), "  INFO: plot ROC per gene for cancer type: ", ct, sep=""))
  colors <- rainbow(length(genes))
  pdf(paste("ROC_",ct,".pdf",sep=""),12,6)
  par(mfrow=c(1,2))
  plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),ylab="Precision",xlab="Recall",bty='n')
  for(gene in genes){
    score = get(paste(ct,"_translate__PRS",gene,sep=""))
    test.labels = get(paste(ct,"_translate__testlabel",gene,sep=""))
    # we can't use FALSE/TRUE as labels when negative scores because it considers s<t as false
    # and in our case we can have true <t (because we labelled with t/f)
    test.labels=rep("GENE", length(test.labelsold))
    test.labels[which(test.labelsold == FALSE)] = "NOGENE"
    pred <- prediction(score, test.labels)
    perf <- performance(pred, "prec", "rec")
    roc.x <- unlist(perf@x.values)
    roc.y <- unlist(perf@y.values)
    lines(roc.y ~ roc.x, col = colors[which(genes == gene)], lwd = 2)
    dat[match(gene, genes),"PREC"] = ifelse(length(which(roc.x > 0.1))==0,
                                            0, max(roc.y[which(roc.x > 0.1)], na.rm=T))
    dat[match(gene, genes),"REC"] = ifelse(length(which(roc.y > 0.7))==0,
                                           0, max(roc.x[which(roc.y > 0.7)], na.rm=T))
  }
  plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),ylab="Sensitivity",xlab="Specificity",bty='n')
  for(gene in genes){
    score = get(paste(ct,"_translate__PRS",gene,sep=""))
    test.labels = get(paste(ct,"_translate__testlabel",gene,sep=""))
    pred <- prediction(score, test.labels)
    perf <- performance(pred, "sens", "spec")
    roc.x <- unlist(perf@x.values)
    roc.y <- unlist(perf@y.values)
    lines(roc.y ~ roc.x, col = colors[which(genes == gene)], lwd = 2)
    dat[match(gene, genes),"AUC"] = round(auc(test.labels, score,
                                              levels = c(FALSE,TRUE), direction = "<"),3)
  }
  legend("bottomleft", legend = genes, col = colors, pch = 19, horiz = F, cex=0.5)
  dev.off()
  if(match(ct, cancertypes) == 1) {all_dat = dat} else {all_dat = rbind(all_dat, dat)}
}

# use the correct cancer type ID
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
                label.select = list(criteria = "`y` > 0.7"),  # Select some labels to display
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
                label.select = list(criteria = "`y` > 0.7"),  # Select some labels to display
                font.label = list(size = 9, face = "italic"), # label font
                repel = TRUE,                                 # Avoid label text overplotting
                legend = "NONE", ylim = c(0,1)
) + geom_hline(yintercept = 0.7, linetype = 2)

p3 <- ggboxplot(all_dat, x = "tissue",
                y = "REC", ylab = "Max Recall(Precision>0.7)",
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
ggsave("boxplot_auc_prec.pdf", plot = p4, width = 15, height = 4)