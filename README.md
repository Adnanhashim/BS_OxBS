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


