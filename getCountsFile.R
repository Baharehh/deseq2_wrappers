##INPUTS 
data<-read.table('/srv/scratch/annashch/stemchell/het/output/remapped_merged_peaks/all.idr.bed.merged',header=F) 
pr_files=c(
"/srv/scratch/annashch/stemcells/het/output/M5_remapped/output/signal/macs2/rep2/M5_rep4_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/M5_remapped/output/signal/macs2/rep1/M5_rep2_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/Hk_remapped/output/signal/macs2/rep2/Hk_rep2_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/Hk_remapped/output/signal/macs2/rep1/Hk_rep1_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/CC_remapped/output/signal/macs2/rep2/CC_rep2_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/CC_remapped/output/signal/macs2/rep1/CC_rep1_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/3hr_remapped/output/signal/macs2/rep2/3hr_rep3_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/3hr_remapped/output/signal/macs2/rep1/3hr_rep2_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/16hr_remapped2/output/signal/macs2/rep1/16hr_rep2_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/16hr_remapped2/output/signal/macs2/rep1/16hr_rep1_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/16hr_remapped2/output/signal/macs2/rep3/16hr_rep3_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/48hr_remapped/output/signal/macs2/rep2/48hr_rep3_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/48hr_remapped/output/signal/macs2/rep1/48hr_rep2_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/H1_remapped/output/signal/macs2/rep2/H1_rep3_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/H1_remapped/output/signal/macs2/rep1/H1_rep1_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig",
"/srv/scratch/annashch/stemcells/het/output/H1_remapped2/output/signal/macs2/rep1/H1_rep2_merged_R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig")
names(pr_files)=c("M5_r4","M5_r2","Hk_r2","Hk_r1","CC_r2","CC_r1","3hr_r3","3hr_r2","16hr_r2","16hr_r1","16hr_r3","48hr_r3","48hr_r2","H1_r3","H1_r1","H1_r2") 

outputfname="/srv/scratch/annashch/stemcells/het/pseudorep_counts.regression.csv"
readlength=51 

#############################################################################

library('GenomicRanges') 
library('rtracklayer') 


colnames(data)<-c("chr","start","end")
bed<-with(data,GRanges(chr,IRanges(start+1,end)))
bw=BigWigFileList(pr_files) 
counts <- matrix(NA, nrow = length(bed), ncol = length(bw))
colnames(counts) <- names(bw)
chromnames=levels(seqnames(bed))

for(i in seq_len(length(bw))) 
{
    last=1 	 
    coverage <- import(bw[[i]], as = 'RleList')
    for(j in chromnames)
    {
    range_vals=ranges(bed[seqnames(bed)==j])
    cur_coverage=coverage[[j]] 
    if(is.null(cur_coverage))
    {
    counts[last:(last+length(range_vals)-1),i]=matrix(0,nrow=length(range_vals),ncol=1)
    }
    else
    {
    newvals=sum(Views(cur_coverage, ranges(bed[seqnames(bed)==j])))
    counts[last:(last+length(newvals)-1), i] <-newvals 
    }   
    last=last+length(range_vals) 
    }

}
## Divide by read length and round to integer numbers
counts <- round(counts / readlength, 0)
write.csv(counts,file=outputfname)

