##########################PAVIRnet#########################################
########################Hypergeometric test###########################

rm(list=ls())
setwd("~/PA/HypergeometricTest/")

TF <- as.matrix(read.delim("TF_list.txt",sep="\t",header=F))[,2]

universe <- c()
net_list <- list.files("regulons/")
for(i in 1 : length(net_list)){
  net <- as.matrix(read.delim(paste("regulons/",net_list[i],sep=""),sep="\t",header=F))
  target.genes <- unique(net[,2])
  universe<- c(universe,target.genes)
}
universe <- union(TF,unique(universe))

universe.number <- length(unique(universe))

hit_list <- as.matrix(read.delim("hits/T6SS_genes_id.txt",sep="\t",header=F))
total.hits <- intersect(unique(hit_list[,2]),universe)
total.hits.number <- length(total.hits)
mra_results <- c()
for(i in 1 : length(net_list)){
  net <- as.matrix(read.delim(paste("regulons/",net_list[i],sep=""),sep="\t",header=F))
  TF <- unlist(strsplit(net_list[i],"_"))[1]
  target.genes <- unique(net[,2])
  target.genes.number <- length(target.genes)
  observed.Hits <- length(intersect(target.genes,hit_list))
  pval <- phyper(observed.Hits - 1, m = total.hits.number, n = universe.number - total.hits.number, k = target.genes.number, lower.tail=F )
  mra_results <- rbind(mra_results,data.frame(TF,universe.number,target.genes.number,total.hits.number,observed.Hits,pval))
}
pvals.adj <- p.adjust(mra_results[,6], method="BH")
mra_results <- cbind(mra_results,pvals.adj)

mra_results <- mra_results[order(mra_results[,6]),]
write.csv(mra_results, file="T6SS_mra_results_all.csv")

##################################################################################
