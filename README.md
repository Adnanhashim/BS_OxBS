# Bisulfite (BS) and oxidataive Bisulfite (oxBS) Sequencing

# AIM

* Check how the lack of DNMT3B alter the epigenetic profile of cells
* Differential landscape of CpG and non-CpG methylation in hESCs caused by absence of DNMT3B. 

## Background & introduction

DNA methylation is an epigenetic modification that modulates chromatin structure, and transcriptional control, and regulates the activation and repression of genes involved in stem cell function and lineage specification. DNA methylation is catalyzed by the DNA methyltransferase (DNMTs) protein family: DNMT1 is responsible for maintenance of DNA methylation, while DNMT3A and DNMT3B act as de novo methyltransferases. Methylated cytosines are found mainly at CpG dinucleotides, but can also be present at non-CpG sites (CHG and CHH; where H corresponds to A, T or C). CpG and non-CpG methylation are found throughout the genome, including repetitive sequences, enhancers, promoters, and gene bodies, and both are involved in regulation of gene expression. Non-CpG methylation is mostly restricted to specific cell types, such as pluripotent stem cells, oocytes, neurons and glial cells. Accordingly, methylated non-CpG sites function as key epigenetic marks in embryonic stem cells and regulate cell type specific functions. Here using whole genome bisulfite and RNA sequencing, we investigated whether and how the absence of the de novo DNA methyltransferase (DNMT) 3B affects CpG and non-CpG methylation and hydroxy-methylation in human embryonic cells (ESC). We provide a single base resolution map of CpG and non-CpG methylated sites in hESCs in the presence or absence of DNMT3B, and further show that loss of DNMT3B not only leads to hypomethylation of CpG but also non-CpG sites


### Bisulfite (BS) and oxidataive Bisulfite (oxBS) Sequencing:


## Sequencing & sample informations

5-methylcytosine (5mc) is well known epigenetic mark that modulates chromatin structure, and transcriptional control, and regulates the activation and repression of genes involved in stem cell function and lineage specification. Abnormal DNA methytlation (hypomethylation or hypermethylation) has been observed in many diseases and cancers. 5mc can be oxidized to 5hmc. Bisulfite sequencing (BS) cannot distinguish between 5mc and 5hmc therefore oxidative bisulfite sequenicng (oxBS) is develeoped that adds an additional enzymatically oxidative step to bisulfite sequencing to discriminate between 5mc and 5hmc marks at single-base resolution in genomic DNA. 

* BS-seq can not dicriminate between 5mc and 5hmc  
  BS seq = 5mc + 5hmc

* In oxBS addtional oxidation step is added to BS-seq to distinguish between 5mc and 5hmc so as a result we get true 5mc.
  oxBS seq = 5mc (true 5mc)

* To obatin the true 5hmc level in geneomic DNA we do subtraction of BS-seq and oxBS-seq.
  BS seq - oxBS seq = 5hmc (true 5hmc)

![booth_et_al](https://github.com/Adnanhashim/BS_OxBS/blob/master/booth_et_al_science_crop.png)
                                                               (Booth et al. 2012, Science)





## Data Analysis:

### Quality control: FASTQC

FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

```bash
fastqc *.fastq -o 1_fastQC
```

### Quality control: cegx-bsExpress (spike-in oligonucleotides controls)
cegx_bsExpress (https://github.com/cegx-ds/cegx_bsExpress) is a tool for quality control of BS-seq and oxBS-seq sequeincing libraries. Spike-in oligonucleotides controls (SQ1hmc, Sq3hmC and SQ6hmc) where amounts of cytosine modificaton are known are used in BS-Seq and oxBS-Seq to access the successful BS and oxBS experiments. cegx_bsExpress can used for BS-seq and oxBS-seq data analysis from raw (fastq files) to genome-wide methylation call but I am using cegx-bsExpress just for the analysis of spike-in oligonuclotide controls. 

```bash
mkdir 1_bsExpress
nohup sh -c 'for fq1 in *R1_001.fastq.gz; do
fq2=${fq1/_R1_/_R2_}
out=`basename $fq1 _R1_001.fastq.gz`
#out1=${out}_R1.fq.gz
#out2=${out}_R2.fq.gz
#echo $fq1, $fq2, $out, $out1, $out2

if [[ ! -f "$fq1" ]]; then echo "WRONG"; break; fi
if [[ ! -f "$fq2" ]]; then echo "WRONG"; break; fi

echo -e "Ruuning Pre­alignment CEGX tailed Sequencing Controls: cegx_bsExpress" 
mkdir 1_bsExpress/$out/
oxBS_controls=/lsc/common/Adnantools/bsExpress-0.4.1b/control_reference/oxBS_controls-v1.0.fa 
bsExpress -i $fq1 $fq2 -r $oxBS_controls -p runqc --outdir 1_bsExpress/$out/
done' > 1_bsExpress/nohup_1_bsExpress_spikein_alignment.txt

```
cegx_bsExpress generates follwoing outputs 
* *.coverage.pdf (Histogram of the read coverage C to T at each known cytosine modfication in spike-in oligonucleotides control sequences)  
* *.conversion.pdf (Histogram of the percentage of unconverted C to T at each known cytosine modfication in spike-in oligonucleotides control sequences)
* *.oxqc_summary.txt (conversion information of C to T) 


Example of Wildtype (WT) quantitative assesment of 5mc and 5hmc percentages per base call in spike-in oliginucleotide controls (SQ1hmc, SQ3hmc and SQ6hmc) 
![WT_spikein_controls](https://github.com/Adnanhashim/BS_OxBS/blob/master/spike-in_control.png)

  
* WT BS summary: 

Spike-in | Conversion | Percentage of conversion | Converted reads | Total reads 
-------- | -----------| ------------------------ |  -------------- | -----------
SQ1hmC | 5hmC | 96.61 | 114 | 118
SQ1hmC | 5mC | 91.42 | 4019 | 4396
SQ1hmC |	C | 0.33 | 10 | 3053
SQ3hmC |	5hmC | 97.25 | 1060 | 1090
SQ3hmC |	5mC | 94.21 |	3790 | 4023
SQ3hmC |	C | 0.42 |	17 | 4087
SQ6hmC |	5hmC |	96.31 |	1226 | 1273
SQ6hmC |	5mC |	68.6 |	830 | 1210
SQ6hmC |	C |	0.47 |	6 | 1266
SQC |	5mC |	96.59 |	6831 | 7072
SQC |	C |	0.16 |	92 | 57487
SQmC | 5mC | 92.58 |	1173 | 1267
all |	5hmC |96.74 |	2400 | 2481
all	| 5mC |	92.63 |	16643 | 17968
all	| C | 0.19 | 125 | 65893  


* WT oxBS summary:
 
Spike-in | Conversion | Percentage of conversion | Converted reads | Total reads 
-------- | -----------| ------------------------ |  -------------- | -----------
SQ1hmC | 5hmC |	3.33 | 2 | 60
SQ1hmC |	5mC |	89.89 |	818 |	910
SQ1hmC |	C |	0	| 0 |	674
SQ3hmC |	5hmC |	4.38 |	12 |	274
SQ3hmC |	5mC |	92.34 |	651 |	705
SQ3hmC |	C |	0 |	0 |	658
SQ6hmC |	5hmC | 3.93 |	28 |	712
SQ6hmC |	5mC	| 90.32 |	541 |	599
SQ6hmC |	C	| 0.37 | 2 |	535
SQC |	5mC	| 95.91 |	2605 |	2716
SQC |	C	| 0.13 |	29 |	22303
SQmC |	5mC |	90.51 |	687 |	759
all |	5hmC |	4.02 |	42 |	1046
all |	5mC |	93.2 |	5302 |	5689
all |	C |	0.13 |	31 |	24170


### Adapter trimming

TrimGlore (https://github.com/FelixKrueger/TrimGalore) was used to trim adapter sequences. 

```bash
nohup sh -c 'for fq1 in *R1.fq.gz; do
fq2=${fq1/R1/R2}
out=`basename $fq1 _R1.fq.gz`
mkdir ./2_adapter_trim
echo $fq1, $fq2, $out
 
if [[ ! -f "$fq1" ]]; then echo "WRONG"; break; fi
if [[ ! -f "$fq2" ]]; then echo "WRONG"; break; fi

trim_galore --paired --illumina $fq1 $fq2 -o ./2_adapter_trim ;done' > ./2_adapter_trim/nohup_trim_glore.txt
```
Example of Wildtype (WT) before and after adapter trimming.
![WT_TrimGlore](https://github.com/Adnanhashim/BS_OxBS/blob/master/trimglore.png)

  
### Alignment: bismark 
Bismark (https://github.com/FelixKrueger/Bismark) was used to align bisulfite and oxidative bisulfite treated sequencing reads to the human genome (hg38).

```bash

mkdir 2_bismark_alignment
nohup sh -c 'for fq1 in *_R1_001_trim.fastq.gz; do
fq2=${fq1/_R1_/_R2_}
out=`basename $fq1 _R1_001_trim.fastq.gz`
#echo $fq1, $fq2, $out

if [[ ! -f "$fq1" ]]; then echo "WRONG"; break; fi
if [[ ! -f "$fq2" ]]; then echo "WRONG"; break; fi

echo -e "Alignment: Paired­End alignment of reads to a reference genome"
mkdir 2_bismark_alignment/$out/

bismark --multicore 12 --bowtie2 --genome_folder /lsc/common/Adnantools/bismark_indexed_genome_hg38_bowtie2 -1 $fq1 -2 $fq2 --output_dir 2_bismark_alignment/$out/ --temp_dir /storage/scratch/Adnan/tmp_dir_bismark &&

#sam to sorted BAM
samtools view -bS 2_bismark_alignment/$out/${out}.sam | samtools sort -@ 16 - -o  2_bismark_alignment/$out/${out}_sorted.bam 

#index of BAM 
samtools index 2_bismark_alignment/$out/${out}_sorted.bam 2_bismark_alignment/$out/${out}_sorted.bam.bai
rm 2_bismark_alignment/$out/${out}.sam; done' > nohup_1_bsExpress_bismark_alignment.txt
```

### Merge replicates & remove duplicates:

Replicates were merged for further analysis. 

```bash
samtools merge  -@ 32 -O BAM *merged.bam *.Rep1.bam *.Rep2.bam
```


### Methylation call

```R
setwd("./")

#load libraries

library(methylKit)
library(genomation)
library(GenomicRanges)
library("EnrichedHeatmap")
library(ggplot2)
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(cluster) # diana
library(circlize) # colorRamp2
library(reshape2)
library(rtracklayer)

# Methylation call at CpGs
processBismarkAln(location=list("WT_bs_1-WT-bs-r1_S11_5-WT-bs-r2_S15_sorted.bam","WT_oxbs_2-WT-oxbs-r1_S12_6-WT-oxbs-r2_S16_sorted.bam","B_KO_bs_3-B-bs-r1_S13_7-B-bs-r2_S17_sorted.bam","B_KO_oxbs_4-B-oxbs-r1_S14_8-B-oxbs-r2_S10_sorted.bam"),sample.id= list("WT_bs","WT_oxbs","B_KO_bs","4B_KO_oxbs" ), assembly="hg38",save.folder = "methyl_call/", save.context = c("CpG"), read.context = "CpG",nolap = TRUE, mincov = 10, minqual = 20, phred64 = FALSE,treatment=c(1,1,0,0), save.db = TRUE)

more WT_oxbs_CpG.txt | egrep -iv "chrY|chrUn*|*random" > WT_oxbs_CpG_mod.txt 
more BKO_oxbs_CpG.txt | egrep -iv "chrY|chrUn*|*random" > BKO_oxbs_CpG_mod.txt

# First of need to convert adjusted 5hmc into the methRead format
more Bko_adjusted_true_5hmC.txt | sed 's/\"//g' | grep -v ^chr| cut -f2,3,5-8|sed 's/\t/\./' |sed 's/\+/F/'|sed 's/\-/R/' |awk -v OFS="\t" '{print $1,$1,$2,$3,100*($4/$3),100*($5/$3)}'| sed 's/\./\t/2' > Bko_adjusted_true_5hmC_tmp.txt
more WT_adjusted_true_5hmC.txt | sed 's/\"//g' | grep -v ^chr| cut -f2,3,5-8|sed 's/\t/\./' |sed 's/\+/F/'|sed 's/\-/R/' |awk -v OFS="\t" '{print $1,$1,$2,$3,100*($4/$3),100*($5/$3)}'| sed 's/\./\t/2' > WT_adjusted_true_5hmC_tmp.txt

echo -e "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT" | cat - WT_adjusted_true_5hmC_tmp.txt > WT_adjusted_true_5hmC_methylkit_input_final.txt 
echo -e "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT" | cat - Bko_adjusted_true_5hmC_tmp.txt > Bko_adjusted_true_5hmC_methylkit_input_final.txt

#True 5mc

dir.create("./CpG")
cpg.files=list("../CpG/5mc/WT_oxbs_CpG_mod.txt","../CpG/5mc/BKO_oxbs_CpG_mod.txt")
mobj.cpg =methRead(cpg.files,sample.id=list("WT","BKO"),assembly="hg38",treatment=c(0,1),context="CpG") #CpG methylation call

filt.mobj.cpg = filterByCoverage(mobj.cpg, lo.count=10, hi.perc=99.9)
save(mobj.cpg, filt.mobj.cpg, file="methylRawList.RData")
#load("methylRawList.RData")

sink("CpG/methylKit.cpg.statistics.txt")
     for (i in c(1:length(mobj.cpg))){
          sample.id = mobj.cpg[[i]]@sample.id
          pPass = nrow(filt.mobj.cpg[[i]])/nrow(mobj.cpg[[i]])
          print(paste(sample.id, ":", sprintf("%.2f",pPass*100),"%"))
          getMethylationStats(filt.mobj.cpg[[i]], plot=F, both.strands=F) # % methylation distribution
          png(paste0("CpG/",sample.id,".cpg.methylation.png"),height=3000,width=3000,res=300) 
          getMethylationStats(filt.mobj.cpg[[i]], plot=T, both.strands=F) # Plot
          dev.off()
          getCoverageStats(filt.mobj.cpg[[i]], plot=F, both.strands=F)
          png(paste0("CpG/",sample.id,".cpg.ori-coverage.png"),height=3000,width=3000,res=300)
          getCoverageStats(mobj.cpg[[i]], plot=T, both.strands=F)
          dev.off()
          png(paste0("CpG/",sample.id,".cpg.coverage.png"),height=3000,width=3000,res=300)
          getCoverageStats(filt.mobj.cpg[[i]], plot=T, both.strands=F)
          dev.off()
    }
sink()

#Differential methylation (hypo_hyper) CpG sites/bases
meth.cpg = unite(filt.mobj.cpg) #common sites
write.table(meth.cpg,"CpG/WT_VS_BKO_filtered_unite_CPG.txt",sep="\t")  #WT_VS_BKO
mDiff.cpg = calculateDiffMeth(meth.cpg)
write.table(meth.cpg,"CpG/output_WT_VS_BKO_Diff_all_CPG.txt",sep="\t")  #WT_VS_BKO
mDiff25p.hyper.cpg = getMethylDiff(mDiff.cpg, difference=25, qvalue=0.01, type="hyper")
mDiff25p.hypo.cpg = getMethylDiff(mDiff.cpg, difference=25, qvalue=0.01, type="hypo")

write.table(mDiff25p.hyper.cpg,"CpG/output_WT_VS_BKO_Diff_hyper_CPG.txt",sep="\t")  #WT_VS_BKO
write.table(mDiff25p.hypo.cpg,"CpG/output_WT_VS_BKO_Diff_hypo_CPG.txt",sep="\t")

mDiff25p.hyper.chg = getMethylDiff(mDiff.chg, difference=25, qvalue=0.01, type="hyper")
mDiff25p.hypo.chg = getMethylDiff(mDiff.chg, difference=25, qvalue=0.01, type="hypo")

write.table(mDiff25p.hyper.chg,"CHG/output_WT_VS_BKO_Diff_hyper_CHG.txt",sep="\t")  #WT_VS_BKO
write.table(mDiff25p.hypo.chg,"CHG/output_WT_VS_BKO_Diff_hypo_CHG.txt",sep="\t")

mDiff25p.hyper.chh = getMethylDiff(mDiff.chh, difference=25, qvalue=0.01, type="hyper")
mDiff25p.hypo.chh = getMethylDiff(mDiff.chh, difference=25, qvalue=0.01, type="hypo")

write.table(mDiff25p.hyper.chh,"CHH/output_WT_VS_BKO_Diff_hyper_CHH.txt",sep="\t")  #WT_VS_BKO
write.table(mDiff25p.hypo.chh,"CHH/output_WT_VS_BKO_Diff_hypo_CHH.txt",sep="\t")

mDiff25p.cpg = getMethylDiff(mDiff.cpg, difference=25, qvalue=0.01)
mDiff25p.chg = getMethylDiff(mDiff.chg, difference=25, qvalue=0.01)
mDiff25p.chh = getMethylDiff(mDiff.chh, difference=25, qvalue=0.01)

write.table(mDiff25p.cpg,"CpG/output_WT_VS_BKO_Diff_significant_CPG.txt",sep="\t")  #WT_VS_BKO
write.table(mDiff25p.chg,"CHG/output_WT_VS_BKO_Diff_significant_CHG.txt",sep="\t")
write.table(mDiff25p.chh,"CHH/output_WT_VS_BKO_Diff_significant_CHH.txt",sep="\t")
  

sink("CpG/methylKit.DMperChr.cpg.txt")
diffMethPerChr(mDiff.cpg, meth.cutoff=25, qvalue.cutoff=0.01, plot=FALSE)
sink()
png("CpG/methylKit.DMperChr.cpg.png",height=3000,width=3000,res=300) 
diffMethPerChr(mDiff.cpg, meth.cutoff=25, qvalue.cutoff=0.01, plot=TRUE)
legend("topright", title="CpG", legend=c("hyper","hypo"), fill=c("magenta","aquamarine4"))
dev.off()

gene.obj=readTranscriptFeatures("../hg38_RefSeq.bed")

annotateWithGeneParts(as(mDiff25p.cpg,"GRanges"),gene.obj)
annotateWithGeneParts(as(mDiff25p.hyper.cpg,"GRanges"),gene.obj)
annotateWithGeneParts(as(mDiff25p.hypo.cpg,"GRanges"),gene.obj)

diffAnn.hypo.cpg=annotateWithGeneParts(as(mDiff25p.hypo.cpg,"GRanges"),gene.obj)
diffAnn.hyper.cpg=annotateWithGeneParts(as(mDiff25p.hyper.cpg,"GRanges"),gene.obj)
write.table(getAssociationWithTSS(diffAnn.hypo.cpg),"CpG/diffAnn_hypo_CpG_associtaion_within_TSS.txt",sep="\t")
write.table(getAssociationWithTSS(diffAnn.hyper.cpg),"CpG/diffAnn_hyper_CpG_associtaion_within_TSS.txt",sep="\t")

getTargetAnnotationStats(diffAnn.hypo.cpg,percentage=TRUE,precedence=TRUE)
getTargetAnnotationStats(diffAnn.hyper.cpg,percentage=TRUE,precedence=TRUE)

png("CpG/DE_hypo_CpG_methylation_annotation_genelevel.png", height = 3000, width = 4800, res = 300)
plotTargetAnnotation(diffAnn.hypo.cpg,precedence=TRUE,main="differential hypo methylation annotation")
dev.off()

png("CpG/DE_hyper_methylation_annotation_genelevel.png", height = 3000, width = 4800, res = 300)
plotTargetAnnotation(diffAnn.hyper.cpg,precedence=TRUE,main="differential hyper methylation annotation")
dev.off()

###################  Differential methylation (hypo_hyper) regions

		###################  CPG

tiles.cpg=tileMethylCounts(meth.cpg,win.size=1000,step.size=1000)
write.table(tiles.cpg,"CpG/WT_VS_BKO_filtered_unite_regions_1000.txt",sep="\t")

tileDiff.cpg = calculateDiffMeth(tiles.cpg)
write.table(tileDiff.cpg,"CpG/output_WT_VS_BKO_all_Regions_CPG.txt",sep="\t")
tileDiff25p.cpg = getMethylDiff(tileDiff.cpg, difference=25,qvalue=0.01)
tileDiff25p.cpg.hyper = getMethylDiff(tileDiff.cpg, difference=25,qvalue=0.01,type="hyper")
tileDiff25p.cpg.hypo = getMethylDiff(tileDiff.cpg, difference=25,qvalue=0.01,type="hypo")
write.table(tileDiff25p.cpg,"CpG/output_WT_VS_BKO_significant_Regions_CPG.txt",sep="\t")
write.table(tileDiff25p.cpg.hyper,"CpG/output_WT_VS_BKO_significan_hyper_Regions_CPG.txt",sep="\t")
write.table(tileDiff25p.cpg.hypo,"CpG/output_WT_VS_BKO_significan_hypo_Regions_CPG.txt",sep="\t")

png("CpG/percentage_DE_hyper_and_hypo_methylated_regions_per_chr_CPG.png", height = 3000, width = 4800, res = 300)
diffMethPerChr(tileDiff.cpg, meth.cutoff=25, qvalue.cutoff=0.01, plot=TRUE)
legend("topright", title="CpG", legend=c("hyper","hypo"),fill=c("magenta","aquamarine4"))
dev.off()

##Annotation of DMR
gene.obj=readTranscriptFeatures("../hg38_RefSeq.bed")
cpg.obj=readFeatureFlank("../cpgi.hg38.bed",feature.flank.name=c("CpGi","shores"))

annotateWithGeneParts(as(tileDiff5p.cpg,"GRanges"),gene.obj)
annotateWithGeneParts(as(tileDiff5p.cpg.hyper,"GRanges"),gene.obj)
annotateWithGeneParts(as(tileDiff5p.cpg.hypo,"GRanges"),gene.obj)

diffAnn.hypo.cpg.tile=annotateWithGeneParts(as(tileDiff5p.cpg.hypo,"GRanges"),gene.obj)
diffAnn.hyper.cpg.tile=annotateWithGeneParts(as(tileDiff5p.cpg.hyper,"GRanges"),gene.obj)
write.table(getAssociationWithTSS(diffAnn.hypo.cpg.tile),"CpG/DMR/diffAnn_hypo_CpG_associtaion_within_TSS.txt",sep="\t")
write.table(getAssociationWithTSS(diffAnn.hyper.cpg.tile),"CpG/DMR/diffAnn_hyper_CpG_associtaion_within_TSS.txt",sep="\t")

getTargetAnnotationStats(diffAnn.hypo.cpg.tile,percentage=TRUE,precedence=TRUE)
getTargetAnnotationStats(diffAnn.hyper.cpg.tile,percentage=TRUE,precedence=TRUE)

png("CpG/DMR/DE_hypo_CpG_methylation_annotation_genelevel.png", height = 3000, width = 4800, pointsize = 12,res = 300)
plotTargetAnnotation(diffAnn.hypo.cpg.tile,precedence=TRUE,main="Hyper-DMR Gene annotation")
dev.off()

png("CpG/DMR/DE_hyper_methylation_annotation_genelevel.png", height = 3000, width = 4800,pointsize = 12, res = 300)
plotTargetAnnotation(diffAnn.hyper.cpg.tile,precedence=TRUE,main="Hyper-DMR Gene annotation")
dev.off()

diffCpGann.hypo.cpg.tile=annotateWithFeatureFlank(as(tileDiff5p.cpg.hypo,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
diffCpGann.hyper.cpg.tile=annotateWithFeatureFlank(as(tileDiff5p.cpg.hyper,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")

png("CpG/DMR/DE_hypo_regions_CPG_methylation_annotation_CpGi_shores_otherRegions.png", height = 3000, width = 4800,pointsize = 12, res = 300)
plotTargetAnnotation(diffCpGann.hypo.cpg.tile,col=c("green","gray","white"),main="Hypo-DMR CpGi/shores annotation")
dev.off()

png("CpG/DMR/DE_hyper_regions_CPG_methylation_annotation_CpGi_shores_otherRegions.png", height = 3000, width = 4800,pointsize = 12, res = 300)
plotTargetAnnotation(diffCpGann.hyper.cpg.tile,col=c("green","gray","white"),main="Hyper-DMR CpGi/shores annotation")
dev.off()



################### Annoation_gene_level	
	
	###################  Average plot

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
transcripts.GR = transcripts(txdb)
elementMetadata(transcripts.GR) = data.frame(name = elementMetadata(transcripts.GR)[,2], stringsAsFactors=FALSE)

promoters.GR = promoters(txdb, upstream=1000, downstream=1000)
elementMetadata(promoters.GR) = data.frame(name = elementMetadata(promoters.GR)[,2], stringsAsFactors=FALSE)
tss.GR = promoters(txdb, upstream=0, downstream=1)
elementMetadata(tss.GR) = data.frame(name = elementMetadata (tss.GR)[,2], stringsAsFactors=FALSE)
exons.GR = unlist(exonsBy(txdb, "tx", use.names=TRUE))

elementMetadata(exons.GR) = data.frame(name = names(exons.GR), stringsAsFactors=FALSE)
names(exons.GR) = NULL
introns.GR = unlist(intronsByTranscript(txdb, use.names=TRUE))

elementMetadata(introns.GR) = data.frame(name = names(introns.GR), stringsAsFactors=FALSE)
names(introns.GR) = NULL
gene.obj_1 = GRangesList("exons" = exons.GR, "introns" = introns.GR, "promoters" = promoters.GR, "TSSes" = tss.GR, "transcripts" = transcripts.GR)
#gene.obj_1

		###################  CPG
perc.meth.cpg = percMethylation(meth.cpg)
write.table(perc.meth.cpg,"CpG/percentage_cpg.txt",sep="\t")
 
GR_unite_mc.cpg = granges(as(meth.cpg, "GRanges"))
GR_unite_WT.cpg = granges(as(meth.cpg, "GRanges"))
GR_unite_BKO.cpg = granges(as(meth.cpg, "GRanges"))
elementMetadata(GR_unite_WT.cpg)$meth.cpg = perc.meth.cpg[,1]
elementMetadata(GR_unite_BKO.cpg)$meth.cpg = perc.meth.cpg[,2]

extend = 5000

TSS_GR_unite_WT.cpg = normalizeToMatrix(GR_unite_WT.cpg, gene.obj$TSSes, value_column = "meth.cpg", mean_mode = "absolute", extend = extend, empty_value = NA)
TSS_GR_unite_BKO.cpg = normalizeToMatrix(GR_unite_BKO.cpg, gene.obj$TSSes, value_column = "meth.cpg", mean_mode = "absolute", extend = extend, empty_value = NA)

mean_TSS_GR_unite_WT.cpg=colMeans(TSS_GR_unite_WT.cpg, na.rm = TRUE)
mean_TSS_GR_unite_BKO.cpg=colMeans(TSS_GR_unite_BKO.cpg, na.rm = TRUE)
dfTSS.cpg = data.frame(Sample = c("WT"," KO"))


dfTSS.cpg[,names(mean_TSS_GR_unite_WT.cpg)] <- rbind(mean_TSS_GR_unite_WT.cpg, mean_TSS_GR_unite_BKO.cpg)
dfTSS.cpg = melt(dfTSS.cpg, id = c("Sample"))
write.table(dfTSS.cpg,"CpG/TSS_average_plot_CPG.txt",sep="\t")
png("CpG/Averageplot_TSS_CpG_methylation_CPG.png",height = 3000, width = 4800, res = 300)
ggplot(dfTSS.cpg, aes(variable, value, group = Sample))+ geom_line(aes(color = Sample), size = 0.6) + theme_light() + geom_line(aes(color = Sample), size = 0.6) + theme_light() + scale_x_discrete(breaks = c("u1","d50"), labels=c("-5000","+5000")) + theme(legend.position="right", axis.title.x = element_text(face = "bold", size = 16), axis.title.y = element_text(face = "bold", size = 16),axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 18)) + labs(title = "Average CpG methylation at promoter region", x ="TSS", y = "% methylation") + annotate("rect", xmin = 50, xmax = 51, ymin = 0, ymax = max(dfTSS.cpg$value)+0.25, alpha = .2)
dev.off()

Tx_WT_cpg = normalizeToMatrix(GR_unite_WT.cpg, gene.obj_1$transcripts, value_column = "meth.cpg", mean_mode = "absolute", extend = extend, empty_value = NA, target_ratio = 0.4)
Tx_BKO_cpg = normalizeToMatrix(GR_unite_BKO.cpg, gene.obj_1$transcripts, value_column = "meth.cpg", mean_mode = "absolute", extend = extend, empty_value = NA, target_ratio = 0.4)

mean_Tx_WT_cpg=colMeans(Tx_WT_cpg, na.rm = TRUE)
mean_Tx_BKO_cpg=colMeans(Tx_BKO_cpg, na.rm = TRUE)
dfTx.cpg = data.frame(Sample = c("WT"," BKO"))
dfTx.cpg[,names(mean_Tx_WT_cpg)] <- rbind(mean_Tx_WT_cpg, mean_Tx_BKO_cpg)
dfTx.cpg = melt(dfTx.cpg, id = c("Sample"))
write.table(dfTx.cpg,"CpG/Transcribed_regions_average_plot_CPG.txt",sep="\t")
png("CpG/Averageplot_transcribed_region_CpG_methylation.png",height = 3000, width = 4800, res = 300)
ggplot(dfTx.cpg, aes(variable, value, group = Sample)) +geom_line(aes(color = Sample), size = 0.6) + theme_light() + scale_x_discrete(breaks = c("u1","d50"), labels= c("-5000","+5000")) + theme(legend.position="right",axis.title.x = element_text(face = "bold", size = 16),axis.title.y = element_text(face = "bold", size = 16),axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),panel.grid.minor = element_blank(),plot.title = element_text(face = "bold", size = 18)) + labs(title = "Average CpG methylation at transcribed region ",x = "Transcribed Region", y = "% methylation") + annotate("rect", xmin = 51, xmax = 117, ymin = 0, ymax = max(dfTx.cpg$value)+0.25, alpha = .2)
dev.off()
