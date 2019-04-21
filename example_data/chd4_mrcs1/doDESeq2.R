
library(DESeq2)

setwd(".")

qval = 0.05
fc_filter =  1 # 2 fold

count_table <- read.delim("raw_table.tsv", row.names=1) # 'row.names must equal first item in [0]
drops <- c('name')
count_table = count_table[, !(names(count_table) %in% drops)]

conds = c('Luc', 'Luc', 'Luc', 'Mcrs1', 'Mcrs1', 'Chd4', 'Chd4')
sam_names = c('shLuc.rep1',	'shLuc.rep2',	'shLuc.rep3',	'shMcrs1.rep1',	'shMcrs1.rep2',	'shChd4.rep1',	'shChd4.rep2')

design = data.frame(shrna=conds)
rownames(design) = sam_names

contrasts1 = c('Luc.vs.Chd4', 'Luc.vs.Mcrs1')
contrasts2 = c('Luc', 'Luc')
contrasts3 = c('Chd4', 'Mcrs1')
contrasts4 = c(0.05, 0.05) # Not used 

contrasts = data.frame(c1=contrasts2, c2=contrasts3, q=contrasts4, stringsAsFactors=FALSE)
rownames(contrasts) = contrasts1

cdsFull <- DESeqDataSetFromMatrix(countData=count_table, 
                                  colData=design, 
                                  design=~shrna)

colData(cdsFull)$shrna = relevel(colData(cdsFull)$shrna, "Luc")
sizeFactors(cdsFull) = rep(1, dim(design)[1]) # Already normalised GC counts, load them as 1
dds2 = estimateDispersions(cdsFull)
dds = nbinomWaldTest(dds2)
#dds = DESeq(cdsFull, fitType='local')

res <- results(dds)#, contrast=c(1,0,0,1,0,0,-1,0,0))
res_names = resultsNames(dds)

pdf(file="DESeq2_results.pdf", width=8, height=7)
par(mfrow=c(1,1))
rest = c("type", "up", "dn")

plotDispEsts(dds)

for (i in 1:dim(contrasts)[1]) {
  print(contrasts[i,])
  
  c1 = contrasts[i,]$c1
  c2 = contrasts[i,]$c2
  res = results(dds, contrast=c("shrna", c1, c2))
  
  #write.table(res, file=paste(row.names(contrasts[i,]), ".all.tsv", sep=""), sep="\t", col.names=NA)

  res = res[complete.cases(res),] # Filter NAs
  pres = res[res$padj<qval,]
  pres = pres[abs(pres$log2FoldChange) > fc_filter, ] # select fold-change here
  dn = pres[pres$log2FoldChange > 0 ,] # Yes, this is weird, but in DESeq2 the up and dn are not
  up = pres[pres$log2FoldChange < -0,] # in the orientation that you might expect.
  
  print(paste("up:", dim(up)[1], "dn:", dim(dn)[1]))
  sub_name = paste(row.names(contrasts[i,]), dim(pres)[1], "DE genes", sep=" ")
  
  # Roll my own plotMA()
  plot(log(res[,1]), res[,2], ylim=c(-6, 6))
  abline(h=c(-fc_filter, fc_filter), col = "blue")
  abline(h=0, col="grey")
  points(log(up[, 1]), up[ ,2], col="red", pch=16, cex=0.9)
  points(log(dn[, 1]), dn[ ,2], col="green", pch=16, cex=0.9)
  title(row.names(contrasts[i,]), sub=sub_name)
  
  #write.table(pres, file=paste("", row.names(contrasts[i,]), ".all.tsv", sep=""), sep="\t", col.names=NA)
  write.table(up, file=paste("", row.names(contrasts[i,]), ".up.tsv", sep=""), sep="\t", col.names=NA)
  write.table(dn, file=paste("", row.names(contrasts[i,]), ".dn.tsv", sep=""), sep="\t", col.names=NA)
}

dev.off()

