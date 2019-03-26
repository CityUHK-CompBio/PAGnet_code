##########################Based on RNA-seq, identified DEGs using DESeq############################
rm(list=ls())
setwd("~/PA/RNA-seq/counts/")

counts <- read.delim("HLR_counts.txt",header=T,sep="\t")
colnames(counts) <- c("gene","liang_wt_R1","liang_wt_R2",
                      "liang_pqsR(mvfR)_R1","liang_pasR(mvfR)_R2",
                      "liang_rhlR_R1","liang_rhlR_R2",
                      "liang_PA4595_R1","liang_PA4595_R2",
                      "liang_vqsM_R1", "liang_vqsM_R2",
                      "Haas_wt_R1", "Haas_wt_R2",
                      "Haas_rpoN_R1", "Haas_rpoN_R2",
                      "Haas_qscR_R1", "Haas_qscR_R2")
rownames(counts) <- counts[,1]
counts <- counts[,-1]
###all 0###

library(DESeq)

mat <- counts[,c("liang_wt_R1","liang_wt_R2",
                 "liang_vqsM_R1", "liang_vqsM_R2")]
row0 <- c()
for(i in 1 : dim(mat)[1]){
  if(sum(mat[i,] == "0") == 2){
    row0 <- c(row0,i)
  }
}
mat <- mat[-row0,]
mat <- mat + 1
type <- factor(c( "control", "control",
                  "knock","knock"))
database <- round(as.matrix(mat))
cds <- newCountDataSet(mat,type)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds,"control","knock")

write.csv(res,file="vqsM_DEGs_all.csv")
write.table(res,file="vqsM_DEGs_all.txt",sep="\t",row.names=T,col.names=T,quote=F)

degs <- res[which((res$pval < 0.05)&(abs(res$log2FoldChange)>1) ),]
write.table(degs,"vqsM_DEGs.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.csv(degs,"vqsM_DEGs.csv")

###########################################################################