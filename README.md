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

![booth_et_al](https://github.com/Adnanhashim/BS_OxBS/blob/master/Booth_et_al_science.png)
                                                               (Booth et al. 2012, Science)





## Data Analysis:

### Quality control: FASTQC

FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

```bash
fastqc *.fastq -o 1_fastQC
```

### Removal of Adapter sequences: 

### Methylation Anlaysis (5mc):
