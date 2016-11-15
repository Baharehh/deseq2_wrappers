library("DESeq2")
library("RColorBrewer")
library("gplots")
library("dplyr")
library("BiocParallel")
library("ggplot2")
library(data.table)

# Metadata

padjCutoff <- 0.01 # If this is not 0.1, we have to tweak DESeq's filtering.
register(MulticoreParam(20))
parallelFlag <- TRUE

countData=data.frame(read.table("counts.txt",header=T,sep='\t'))
colData=data.frame(read.table("counts.metadata",header=T,sep='\t'))

ddsFullCountTable<-DESeqDataSetFromMatrix(
	countData=countData,
	colData=colData,
	design= ~Condition+Density+CellCyclePoint+CellCyclePoint*Condition+Density*Condition+CellCyclePoint*Density+Condition*CellCyclePoint*Density) 
dds<-DESeq(ddsFullCountTable,parallel=parallelFlag)

browser() 	

# Generate differential openness matrix
# Each column is a pair of (treated, controls) conditions 
# as specified in conditionsToCompare
numCols <- nrow(conditionsToCompare)
numRows <- nrow(countData)
diffMat <- matrix(, ncol = numCols, nrow = numRows)
confidenceMatt<-matrix(,ncol=numCols, nrow=numRows)
foldChangeMat<-matrix(,ncol=numCols,nrow=numRows) 


for (i in 1:numCols){
  res <- results(dds, 
                 contrast = c("condition", 
                              conditionsToCompare[i, 'treated'], 
                              conditionsToCompare[i, 'controls']), 
                 parallel = parallelFlag)
  low_deltas=which(abs(res$log2FoldChange)<0.15)
  res$log2FoldChange[low_deltas]=0
  diffMat[, i] <- (res$padj <= padjCutoff) * sign(res$log2FoldChange)
  confidenceMatt[,i]<-res$padj
  foldChangeMat[,i]<-res$log2FoldChange
}

#summary(res)
#plotMA(res)

# Replace NAs and filter out all rows of diffMat that are all zero
# Remember to take abs so that +1s and -1s don't cancel each other out
diffMat[is.na(diffMat)] <- 0
confidenceMatt[is.na(confidenceMatt)]<-1
foldChangeMat[is.na(foldChangeMat)]<-0 
#idxToKeep <- rowSums(abs(diffMat)) > 0
#diffMat <- diffMat[idxToKeep, ]
colnames(diffMat) <- conditionsToCompare[, 'both']
colnames(confidenceMatt)<-conditionsToCompare[,'both'] 
colnames(foldChangeMat)<-conditionsToCompare[,'both'] 

peakNames <- read.table(
  sprintf("ppr.merged.bed", peakCutoff),header=TRUE,
  colClasses = c("character"))
alldata=data.frame(peakNames,as.data.frame(diffMat))
allconfidence=data.frame(peakNames,as.data.frame(confidenceMatt))
allfoldchange=data.frame(peakNames,as.data.frame(foldChangeMat))

# Write diffMat to disk
write.table(alldata,
            file = sprintf(
              "diffMat_multifactor.csv",
              padjCutoff, 
              peakCutoff),sep='\t')
#Write confidenceMatt to disk
write.table(allconfidence,
            file = sprintf(
              "confidenceMat_multifactor.csv",
              padjCutoff, 
              peakCutoff),sep='\t')
write.table(allfoldchange,
            file = sprintf(
              "foldChange_multifactor.csv",
              padjCutoff, 
              peakCutoff),sep='\t')
