library(GenomicFeatures)
library(tximport)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(DESeq2)

dir <-file.path(getwd(), "Data", "Counts_Kallisto")
list.files(dir)

samples<-list.files(dir)
files <- file.path(dir, samples, "abundance.tsv")
names(files) <- samples
all(file.exists(files))

txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene<-na.omit(tx2gene)

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
txi$counts <-round(txi$counts,0) #Need to round to integers to get rid of the multimapping for deseq2
write.table(txi$counts, "tximport_counts.txt")

#Load metadata
samples <- read.table(file.path("SampleDesc.txt"), header=TRUE)

#Needed for sorting
samples$Sample <- factor(samples$Sample)
samples$SampleType <- factor(samples$SampleType)
samples$Pairwise <- factor(samples$Pairwise)
samples$State <- factor(samples$State)

#samples<-samples[order(samples$Sample),]
samples<-samples[match(c(colnames(txi$counts)), samples$Sample),]
rownames(samples) <- samples$Sample

#Changes pairwise
dds <- DESeqDataSetFromMatrix(countData = txi$counts,
                              colData = samples,
                              design= ~ Pairwise)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$Pairwise <- relevel(dds$Pairwise, ref = "WT-Control")
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)

#Pairwise contrasts can be made using (coefficients in resultsNames):
#e.g. WT treated vs untreated
BC <- results(dds, contrast=c("Pairwise", "WT-Clonidine", "WT-Control"),
  independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE, name="BC")
BC <- lfcShrink(dds, coef="Pairwise_WT.Clonidine_vs_WT.Control", type="apeglm")
write.table(BC, "Pairwise_WT-Clonidine_vs_WT-Control.txt", sep="\t")


dds <- DESeqDataSetFromMatrix(countData = txi$counts,
                              colData = samples,
                              design= ~ Pairwise)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$Pairwise <- relevel(dds$Pairwise, ref = "KO-Control")
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)

BD <- results(dds, contrast=c("Pairwise", "KO-Clonidine", "KO-Control"),
              independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE, name="BC")
BD <- lfcShrink(dds, coef="Intercept", type="apeglm")
write.table(BD, "Pairwise_KO-Clonidine_vs_KO-Control.txt", sep="\t")
rm(dds)
rm(dds)
