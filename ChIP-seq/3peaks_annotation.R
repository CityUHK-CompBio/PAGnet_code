#######################################PAVIRnet##########################################
#################################peak annotation########################################
##########################1. annotated peaks to PA genome################################
rm(list=ls())
library(ChIPpeakAnno)
library(GenomicFeatures)
library(GenomicRanges)
setwd("~/workspace/DrDeng/PAO1_v4/")

peak_list <- list.files("2annotate_promoters/narrowPeak/")

sss <- read.delim("1redefine_promoter/PAO1_gene_region.txt",sep="\t",header=F)
sss<-sss[,c(1,2,3,4,5)]

colnames(sss)<-c("chr","start","end","strand","id")
ge <- as(sss,"GRanges")
for(s in peak_list){
  pe <- read.delim(paste("2annotate_promoters/narrowPeak/",s,sep=""),sep="\t",header=F)
  sn <- rep("chromosome",length(pe[,1]))
  peaks1 <- GRanges(seqnames=as.matrix(sn),
                    ranges=IRanges(start=as.numeric(as.matrix(pe[,2])),
                                   end= as.numeric(as.matrix(pe[,3])),
                                   names=as.matrix(as.matrix(pe[,4]))))
  
  annotatedPeak <- annotatePeakInBatch(peaks1, AnnotationData=ge)
  df <- data.frame(as.data.frame(annotatedPeak), as.data.frame(values(ge[as.numeric(values(annotatedPeak)$feature),])))
  write.table(df, paste("4direct_anno/anno/",s,"_annotation.txt",sep=""),quote=FALSE, row.names=FALSE, sep="\t")
}

##############################################################################################################
#####################################2 selected overlapStart regions#########################################
rm(list=ls())
setwd("~/workspace/DrDeng/PAO1_v4/4direct_anno/")
flist <- list.files("anno")
for(i in 1 : length(flist)){
  files <- read.delim(paste("anno/",flist[i],sep=""),header=T,sep="\t")
  fname <- unlist(strsplit(flist[i],"_"))[1]
  for(j in 1 : dim(files)[1]){
    if(files[j,"insideFeature"] == "overlapStart"){
      write.table(files[j,],paste("overlapStart/",flist[i],"_overlapStart.txt",sep=""),row.names=F,col.names=F,quote=F,append=T,sep="\t")
      write.table(files[j,"id"],paste("ChIP_targets/",fname,"_ChIP_targets.txt",sep=""),row.names=F,col.names=F,quote=F,append=T,sep="\n")
    }
  }
}

############################################################################################################
######################################3 Considerated operon#################################################
rm(list=ls())
setwd("~/workspace/DrDeng/PAO1_v4/")
library(GenomicRanges)
library(ChIPpeakAnno)
prom_file <- read.delim("1redefine_promoter/PAO1_gene_with_promoter.txt",sep="\t",header=F)
head(prom_file)
prom_region <- prom_file[,c(1,6,7,5)]
colnames(prom_region) <- c("chr","start","end","name")
prom_region <- as(prom_region,"GRanges")

for(list in list.files("2annotate_promoters/narrowPeak")){
  state_file <- read.delim(paste("2annotate_promoters/narrowPeak/",list,sep=""),header=F,sep="\t")
  state_cand <- state_file[,c(2,3)]
  state_chrom <- rep("chromosome",dim(state_cand)[1])
  state_cand <- cbind(state_chrom,state_cand)
  colnames(state_cand) <- c("chr","start","end")
  state_region <- as(state_cand,"GRanges")
  ov <- findOverlaps(prom_region,state_region,select="first")
  inter <- which(ov != "NA")
  pre_promoter <- as.data.frame(prom_region[inter,])
  #pre_enhancer <- state_cand[-inter,]
  #write.table(pre_enhancer,paste("./enhancer/enhancer_",list,sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  write.table(pre_promoter,paste("2annotate_promoters/promoter_",list,sep=""),row.names=F,col.names=F,quote=F,sep="\t")
}

#################################################################################################################
#######################################4 combine gene information and peak information###########################

rm(list=ls())
setwd("~/workspace/DrDeng/PAO1_v4/")

prom_file <- read.delim("1redefine_promoter/PAO1_gene_with_promoter.txt",sep="\t",header=F)
head(prom_file)
prom_region <- prom_file[,c(1,6,7,5)]
colnames(prom_region) <- c("chr","start","end","name")

for(list in list.files("2annotate_promoters/narrowPeak/")){
  state_file <- read.delim(paste("2annotate_promoters/narrowPeak/",list,sep=""),header=F,sep="\t")
  state_cand <- state_file[,c(1,2,3)]
  colnames(state_cand) <- c("chr","start","end")
  for(i in 1 : dim(state_cand)[1]){
    peak_start <- as.numeric(as.vector(state_cand[i,2]))
    peak_end <- as.numeric(as.vector(state_cand[i,3]))
    peak_region <- c(peak_start:peak_end)
    for(j in 1 : dim(prom_region)[1]){
      prom_start <- as.numeric(as.vector(prom_region[j,2]))
      prom_end <- as.numeric(as.vector(prom_region[j,3]))
      prom_name <- prom_region[j,4]
      promoter_region <- c(prom_start : prom_end)
      if(length(intersect(promoter_region,peak_region) != 0)){
        want <- data.frame(state_file[i,],prom_file[j,])
        write.table(want,paste("2annotate_promoters/promoter_with_info/promoter_with_info_",list,sep=""),row.names=F,col.names=F,quote=F,sep="\t",append=T)
      }
    }
  }
}
#################################################5 Find peak and gene pairs###################################################################
rm(list=ls())
setwd("~/workspace/DrDeng/PAO1_v4/")
for(list in list.files("2annotate_promoters/promoter_with_info/")){
  file <- read.delim(paste("2annotate_promoters/promoter_with_info/",list,sep=""),header=F,sep="\t")
  genes <- unique(as.matrix(file[,15]))
  for(i in 1 : length(genes)){
    genei <- genes[i]
    log2FE <- 0
    peak <- c()
    num <- 0
    for(j in 1 : dim(file)[1]){
      if(as.matrix(file[j,15]) == genei){
        log2FE <- log2FE + as.numeric(as.vector(file[j,7]))
        peak <- paste(peak,as.matrix(file[j,4]),sep=" ")
        num <- num + 1
      }
    }
    med <- log2FE / num
    want <- data.frame(genei,log2FE,med,num,peak)
    write.table(want,paste("3gene_peak_pairs/pairs/log2FE_",list,sep=""),sep="\t",row.names=F,col.names=F,quote=F,append=T)
  }
}


####################################################################################################################################
###############################################Obtained final peak annotated genes (upstream and operon)#############################
rm(list=ls())
setwd("~/workspace/DrDeng/PAO1_v4/")
flist <- list.files("3gene_peak_pairs/pairs")
for(i in 1 : length(flist)){
  fname <- unlist(strsplit(flist[i],"_"))[5]
  files <- as.matrix(read.delim(paste("3gene_peak_pairs/pairs/",flist[i],sep=""),header=F,sep="\t"))
  targets <- unique(files[,1])
  ifile <- as.matrix(read.delim(paste("4direct_anno/ChIP_targets/",fname,"_ChIP_targets.txt",sep=""),sep="\n",header=F))
  target <- unique(union(targets,ifile))
  write.table(target,paste("5merge_targets/",fname,"_merge_targets.txt",sep=""),sep="\n",row.names=F,col.names=F,quote=F)
}

################################################################################################################################
