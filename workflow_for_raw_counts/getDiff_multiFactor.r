library("DESeq2")
library("RColorBrewer")
library("gplots")
library("dplyr")
library("BiocParallel")
library("ggplot2")
library(data.table)

# Metadata

padjCutoff <- 0.01/20 
foldChangeCutoff <- 0.10
register(MulticoreParam(32))
parallelFlag <- TRUE

countData=data.frame(read.table("counts.txt",header=T,sep='\t'))
colData=data.frame(read.table("counts.metadata",header=T,sep='\t'))
colData$Density=relevel(colData$Density,ref="LD")


#effect of high vs low density
#we use density & condition as factors and include the interaction term 
ddsFullCountTable<-DESeqDataSetFromMatrix(
	countData=countData,
	colData=colData,
	design= ~Density+Condition+Density:Condition)
dds<-DESeq(ddsFullCountTable,parallel=parallelFlag,betaPrior=FALSE)


res_density=results(dds, name="Density_HD_vs_LD")
png("maplot.density.png")
plotMA(res_density,ylim=c(-1,1))
dev.off()


#Treatment effect in low density samples
res_treatment_ld= results(dds, contrast=c("Condition","treated","control"))
png("maplot.treatment.ld.png")
plotMA(res_treatment_ld,ylim=c(-1,1))
dev.off() 

#Treatment effect for high density samples
res_treatment_hd= results(dds, list( c("Condition_treated_vs_control","DensityHD.Conditiontreated") ))
png("maplot.treatment.hd.png")
plotMA(res_treatment_hd,ylim=c(-1,1))
dev.off()

#THIS IS AN EXAMPLE OF A SINGLE FACTOR ANALYSIS -- SAMPLE IS THE SINGLE INDEPENDENT VARIABLE. 
ddsIndividualFullCountTable<-DESeqDataSetFromMatrix(
	countData=countData,
	colData=colData,
	design= ~Sample)
ddsIndividual<-DESeq(ddsIndividualFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
#create matrices to store the results from the PER-SAMPLE univariate comparison
#get the individual comparisons
conditionsToCompare <- data.frame(
  matrix(
    c('earlyG1.HD.treated', 'earlyG1.HD.controls',
      'earlyG1.LD.treated', 'earlyG1.LD.controls',
      'lateG1.HD.treated','lateG1.HD.controls',
      'lateG1.LD.treated','lateG1.LD.controls',
      'SG2M.HD.treated','SG2M.HD.controls',
      'SG2M.LD.treated','SG2M.LD.controls',      
      'lateG1.HD.treated','earlyG1.HD.treated',
      'lateG1.LD.treated','earlyG1.LD.treated',
      'SG2M.HD.treated','lateG1.HD.treated',
      'SG2M.LD.treated','lateG1.LD.treated',
      'SG2M.HD.treated','earlyG1.HD.treated',
      'SG2M.LD.treated','earlyG1.LD.treated'), 
    ncol = 2,
    byrow = TRUE),
  stringsAsFactors = FALSE)

numCols <- 15
numRows <- nrow(countData)
diffMat <- matrix(, ncol = numCols, nrow = numRows)
confidenceMatt<-matrix(,ncol=numCols, nrow=numRows)
diffMat<-matrix(,ncol=numCols,nrow=numRows) 
diffMat[,1]=(res_density$padj<=padjCutoff)*(abs(res_density$log2FoldChange) >= foldChangeCutoff)*(res_density$log2FoldChange) 
diffMat[,2]=(res_treatment_ld$padj <=padjCutoff)*(abs(res_treatment_ld$log2FoldChange)>=foldChangeCutoff)*(res_treatment_ld$log2FoldChange) 
diffMat[,3]=(res_treatment_hd$padj <=padjCutoff)*(abs(res_treatment_hd$log2FoldChange)>=foldChangeCutoff)*(res_treatment_hd$log2FoldChange) 

colNameEntries=c("HD_vs_LD","LD_treated_vs_control","HD_treated_vs_control")
for (i in 1:12){
  res <- results(ddsIndividual, 
                 contrast = c("Sample", 
                              conditionsToCompare[i, 1], 
                              conditionsToCompare[i, 2]), 
                 parallel = parallelFlag)
  colNameEntries=append(colNameEntries,paste(conditionsToCompare[i,1],"_vs_",conditionsToCompare[i,2],sep=""))
  #plot the MAplot
  #png(paste(conditionsToCompare[i,1],"_vs_",conditionsToCompare[i,2],".png",sep=""))
  #plotMA(res,ylim=c(-1,1))
  #dev.off() 
  diffMat[, 3+i] <- (res$padj <= padjCutoff)*(abs(res$log2FoldChange) >=foldChangeCutoff)*(res$log2FoldChange)
}
# Replace NAs 
# Remember to take abs so that +1s and -1s don't cancel each other out
diffMat[is.na(diffMat)] <- 0
colnames(diffMat) <- colNameEntries 

peakNames <- read.table("ppr.merged.bed",header=TRUE,colClasses = c("character"))
alldata=data.frame(peakNames,as.data.frame(diffMat))
# Write diffMat to disk
write.table(alldata,file = "diffMat_multifactor.STRINGENT.csv",sep='\t')


#IF YOU WANT TO QC THE REPLICATES, YOU CAN PLOT RLOG-TRANSFORMED REPLICATE PAIRS AND CHECK CORRELATION
log-transform and plot pairs of replicates
rld<-rlog(ddsIndividual)
for(i in seq(1,24,2){
      png(paste(as.character(i),"_",as.character(i+1),".png",sep=""))
      plot( assay(rld)[, i:i+1], col="#00000020", pch=20, cex=0.3, xlab=colData[i,"Sample"],ylab=colData[i+1,"Sample"])
      dev.off()
     }	

