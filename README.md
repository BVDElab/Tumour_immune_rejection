# Tumour immune rejection triggered by activation of a2-adrenergic receptors
Analyses supporting the article "Tumour immune rejection triggered by activation of a2-adrenergic receptors"


## single-cell RNA sequencing
2 samples were processed: clonidine and untreated for MC38-OVA mice with tumour. Tumour tissue was extracted and sent for sequencing. Samples were processed with Scanpy, using the external algorithm PhenoGraph for subpopulation detection. Clusters were manually annotated, differential testing was performed with wilcoxon rank sum tests implemented in Diffxpy. Enrichment analysis has been conducted with GSEApy.

To run, execute the Python script Main.py. Extract the two ZIP files from the GEO project and place the corresponding folders in "Data".
You will also need to place the set enrichment definitions in a folder named "Bader_GSEA_GMTs". These you can find here: http://baderlab.org/GeneSets


### Dependencies
diffxpy 0.7.4

gseapy 0.10.5

matplotlib 3.2.1

numpy 1.21.0

pandas 1.1.5

scanpy 1.7.0

seaborn 0.11.1


### Data sources
All datasets have been uploaded under accession GSE....




## RNA sequencing
9 bulk samples of mouse (biological triplicates) were analyses in a reference design, using WT-control as reference.
The follow sample identifiers were present:


 ID | Type | Treatment | SampleGroup
--- | --- | --- | ---
964-1	| WT	| Clonidine	| WT-Clonidine
964-2	| WT	| Clonidine	| WT-Clonidine
964-3	| WT	| Clonidine	| WT-Clonidine
964-4	| WT	| Control	| 	WT-Control
964-5	| WT	| Control	| WT-Control
964-6	| WT	| Control	| WT-Control
964-10	| KO	| Clonidine	| KO-Treated
964-11	| KO	| Clonidine	| KO-Treated
964-12	| KO	| Clonidine	| KO-Treated
964-13	| KO	| Control	| KO-Control
964-14	| KO	| Control	| KO-Control
964-15	| KO	| Control	| KO-Control

Samples were processed on a Windows machine, using Ubuntu 20.02 through WSL2.

```Running trim_galore
find  path_2_fastqfiles  -name "*_1.fastq.gz" | cut -d "_" -f1 | parallel -j 40 trim_galore --illumina --paired --fastqc -o trim_galore/ {}\_1.fastq.gz {}\_2.fastq.gz
```

A pre-built index file for mouse GRCm38.96 was downloaded from https://github.com/pachterlab/kallisto-transcriptome-indices/releases
Subsequently, all paired fastq files were mapped with (default settings):

```Running Kallisto
kallisto quant -i mousereferencefile -o output folder -t 60 fastq1.fastq fastq2.fastq
```

To process the data and obtain the intermediate files available here, set your work directory to the folder Counts_Kallisto in "Data":
```R set work directory
setwd("~Your/Path/To/RNAseq/")
```
Then run the corresponding R script. It will list all files. If all corresponding count files were found, it will print "TRUE" as output.

Raw data has been uploaded on GEO under the accession GSExxxxx. Quantification results by Kallisto can be found in the RNAseq/Data folder.


### Used versions

Apleglm 1.12.0

DESeq 1.30.1

tximport 1.18.0

TXDB Mus Musculus 3.10



## Survival analysis
### Overview
Script to generate the Kaplan-Meier curves provided as Supplementary figure 4D.
Overall survival (OS) obtained from the Xena database was matched to the processed data from the TCGA consortium for the paper "The immune landscape of cancer".

### Versions used

• lifelines 0.25.6

• matplotlib 3.2.1

• numpy 1.21.0

• pandas 1.1.5

• seaborn 0.11.1 



### Data sources
Place the following files in the corresponding "Data" folder:

Cancer immune atlas (Thorsson et al., 2018): https://www.sciencedirect.com/science/article/pii/S1074761318301213

• RNA (final): https://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611

Xena (Goldman et al., 2020): https://doi.org/10.1038/s41587-020-0546-8)

### Code for extended figures

- On Zenodo with DOI 10.5281/zenodo.7673031

• Survival: https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA_survival_data
