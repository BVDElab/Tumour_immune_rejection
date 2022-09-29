library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(scuttle)
library(AnnotationHub)
library(tidyverse)
library(robustbase)
library(pheatmap)
library(celldex)
library(SingleR)
library(bluster)

## ----------------------------------
## Prepare data

if (!file.exists("S04_sce_z46.rds")) {
    z46 <- read10xCounts(c("data/scd_220628_filtered/Z4", "data/scd_220628_filtered/Z6"))
    z46$Group <- ifelse(grepl("Z4", z46$Sample), "WT.CTRL", "WT.CLD")

    ## Add chromosomes (sequence names) for MT QC
    ah <- AnnotationHub()
    ## query(ah, "EnsDb.Mmusculus")
    ## ens <- ah[["AH104895"]] ## Ensembl 107 - 74 missing
    ens <- ah[["AH100674"]] ## Ensembl 106 - 66 missing

    gi <- select(ens, keytype = "GENEID", keys = rownames(z46),
                 column = c("GENEID", "SEQNAME", "SYMBOL"))
    rownames(gi) <- gi$GENEID
    rd <- rowData(z46)
    rd <- merge(rd, gi,
                by.x = "ID", by.y = "GENEID")
    rownames(rd) <- rd$ID

    ## Keep only genes that have a location. Those that don't are
    ## probably retired.
    z46 <- z46[rownames(z46) %in% rd[[1]], ]
    rowData(z46) <- rd[rownames(z46), ]
    rowData(z46)$is_mito <- rowData(z46)$SEQNAME == "MT"

    ## Remove features that haven't been detected
    ## table(rowSums(counts(z46)) > 1)
    ## FALSE  TRUE
    ##  9838 22381
    z46 <- z46[rowSums(counts(z46)) > 1, ]

    z46 <- scater::logNormCounts(z46)
    z46 <- scater::runPCA(z46, exprs_values = 2)

    stopifnot(validObject(z46))
    ## Save object
    saveRDS(z46, file = "S04_sce_z46.rds")

} else {
    z46 <- readRDS("S04_sce_z46.rds")
}

table(z46$Group)
tapply(colSums(assay(z46)), z46$Group, summary)

## ----------------------------------
## Per cell quality control

z46 <- addPerCellQCMetrics(z46, percent.top = 25,
                           subsets = list(Mito = rowData(z46)$is_mito))


df <- perCellQCMetrics(z46, percent.top = 25,
                       subsets = list(Mito = rowData(z46)$is_mito))

reasons <- perCellQCFilters(df, sub.fields = "subsets_Mito_percent")

colSums(as.matrix(reasons))
## low_lib_size            low_n_features high_subsets_Mito_percent
##            0                         0                      3729
##      discard
##         3729

attributes(reasons$high_subsets_Mito_percent)

tapply(colData(z46)$subsets_Mito_percent, z46$Group, function(x) table(x > 15))
## $WT.CLD

## FALSE  TRUE
## 14658  1750

## $WT.CTRL

## FALSE  TRUE
##  6401  2018

stats <- cbind(log10(df$sum), log10(df$detected),
               df$subsets_Mito_percent)
outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")

reasons <- cbind(reasons, multi.outlier)
table(reasons[, 4:5])

colData(z46)$discard <- reasons$discard

plotPCA(z46, colour_by = "Group")

plotPCA(z46, colour_by = "discard", shape_by = "Group", point_size = 3)

gridExtra::grid.arrange(
               plotColData(z46, x = "Group", y = "sum", colour_by = "discard") +
               scale_y_log10() + ggtitle("Library size"),
               plotColData(z46, x = "Group", y = "detected", colour_by = "discard") +
               scale_y_log10() + ggtitle("Detected genes"),
               plotColData(z46, x = "Group", y = "subsets_Mito_percent", colour_by = "discard")
               + ggtitle("Percentage of MT genes"))

plotColData(z46, x = "sum", y = "subsets_Mito_percent",
            colour_by = "discard")
plotColData(z46, x = "sum", y = "subsets_Mito_percent",
            colour_by = "Group")


## perFeatureQCMetrics(z46, subsets = list(CLD = grep("CLD", z46$Group),
##                                         CTRL = grep("CTRL", z46$Group)))

## There's a lot cells and features with very few counts.

## Exploring upperleft cluster
z46$upcl <- reducedDims(z46)[[1]][, 1] < -8 & reducedDims(z46)[[1]][, 2] > 8

plotPCA(z46, colour_by = "discard")
plotPCA(z46, colour_by = "upcl")

## ## Filter out poor quality cells
## z46fl <- z46[, !z46$discard]

## NB: one could argue that it would be good (and easier) to proceed
## with the filtered cells from here on. See second normalisation
## section below.

## Filtering from Python scripts:
## - cells
sel1 <- z46$subsets_Mito_percent < 5    ## percent mito count < 5%
sel2 <- colSums(counts(z46) > 0) > 200  ## at least 200 expressed genes per cell
## - genes
sel3 <- rowSums(counts(z46) > 0) > 100  ## gene must be expressed in at least 100 genes
sel4 <- rowSums(counts(z46) > 1) < 2500 ## remove cells with > 2500 expressed genes

table(sel1, sel2)
##         sel2
## sel1    FALSE  TRUE
##   FALSE     9 11507
##   TRUE      3 13308

table(sel3, sel4)
##         sel4
## sel3    FALSE  TRUE
##   FALSE     0  8929
##   TRUE   1398 12054

## Not filtering genes, as this seems to have an important effect of
## the HGV modelling later.

z46fl <- z46[, z46$subsets_Mito_percent < 10]

table(z46fl$Group)

## ----------------------------------
## Normalisation after cell filtering

z46fl <- scater::logNormCounts(z46fl) ## overwrite logcounts

set.seed(100)
z46fl$qclfl  <- quickCluster(z46fl)

z46fl <- scater::runPCA(z46fl, exprs_values = 2)

plotPCA(z46fl, colour_by = "Group")
plotPCA(z46fl, colour_by = "qclfl")

z46fl <- computeSumFactors(z46fl, cluster = z46$qclfl, min.mean=0.1)
z46fl <- logNormCounts(z46fl) ## overwrites the logcounts assay

plotPCA(z46fl, colour_by = "Group")
plotPCA(z46fl, colour_by = "qclfl")

## ----------------------------------
## Feature selection

hvg <- modelGeneVar(z46fl)

fit_hvg <- metadata(hvg)
plot(fit_hvg$mean, fit_hvg$var,
     xlab = "Mean of log-expression",
    ylab = "Variance of log-expression")
curve(fit_hvg$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

hvg[order(hvg$bio, decreasing = TRUE), ]

plot(fit_hvg$mean, fit_hvg$var,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression",
     col = ifelse(hvg$p.value < 0.1, "red", "black"))
curve(fit_hvg$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

hvg_pois <- modelGeneVarByPoisson(z46fl)

plot(hvg_pois$mean, hvg_pois$total,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")
curve(metadata(hvg_pois)$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

topHvg <- getTopHVGs(hvg, n=1000)

rowSubset(z46fl) <- topHvg

## ----------------------------------
## PCA

z46fl <- runPCA(z46fl, subset_row = topHvg) ## Overwrites PCA
## z46fl <- runPCA(z46fl)
plotPCA(z46fl, colour_by = "Group")

## other function

z46fl <- fixedPCA(z46fl, rank = 100, subset.row = topHvg)

per_var <- attr(reducedDim(z46fl), "percentVar")
plot(per_var, log="y", xlab="PC", ylab="Variance explained (%)")
plot(cumsum(per_var), ylab = "Cumulative variance")

## Litte variance explained in first dimensions, and then only slowly
## increasing - 100 first PCs explain only over 40% variance.

plotReducedDim(z46fl, dimred = "PCA", colour_by="Group")

plotReducedDim(z46fl, dimred = "PCA",
               ncomponents = 4,
               colour_by = "Group")

set.seed(123)
z46fl <- runTSNE(z46fl, dimred = "PCA")
z46fl <- runUMAP(z46fl, dimred = "PCA")


gridExtra::grid.arrange(
               plotReducedDim(z46fl, dimred = "PCA", colour_by = "Group"),
               plotReducedDim(z46fl, dimred = "TSNE", colour_by = "Group"),
               plotReducedDim(z46fl, dimred = "UMAP", colour_by = "Group"))

## save to save long runs
saveRDS(z46fl, file = "S04_sce_z46fl.rds")

## ----------------------------------
## Clustering

## Parameters here used to get around 25 clusters, to compared them
## with the clusters/markers from the first experiment.
nn.clusters <-
    clusterCells(z46fl, use.dimred = "PCA",
                 BLUSPARAM = bluster::SNNGraphParam(k = 3,
                                                    type = "jaccard",
                                                    cluster.fun = "louvain"))
nlevels(nn.clusters)
table(nn.clusters)

plotReducedDim(z46fl, "TSNE", colour_by = I(nn.clusters))

set.seed(100)
kgraph.clusters <- clusterCells(z46fl,
                                use.dimred = "PCA",
                                BLUSPARAM = TwoStepParam(
                                    first = KmeansParam(centers = 2000,
                                                        iter.max = 100L),
                                    second = NNGraphParam(k=10)
                                ))

nlevels(kgraph.clusters)
table(kgraph.clusters)

ctrl <- z46fl$Group == "WT.CTRL"
cld <- z46fl$Group == "WT.CLD"

gridExtra::grid.arrange(
               plotTSNE(z46fl[, ctrl],
                        colour_by = I(kgraph.clusters[ctrl]),
                        text_by=I(kgraph.clusters[ctrl])) +
               ggtitle("CTRL (kgraph)") + theme (legend.position="none"),
               plotTSNE(z46fl[, cld],
                        colour_by = I(kgraph.clusters[cld]),
                        text_by=I(kgraph.clusters[cld])) +
               ggtitle("CLD (kgraph)") + theme (legend.position="none"),
               plotTSNE(z46fl[, ctrl],
                        colour_by = I(nn.clusters[ctrl]),
                        text_by=I(nn.clusters[ctrl])) +
               ggtitle("CTRL (nn)") + theme (legend.position="none"),
               plotTSNE(z46fl[, cld],
                        colour_by = I(nn.clusters[cld]),
                        text_by=I(nn.clusters[cld])) +
               ggtitle("CLD (nn)") + theme (legend.position="none"))

# Using a large pseudo-count for a smoother color transition
# between 0 and 1 cell in each 'tab'.
tab <- table(paste("nn", nn.clusters),
             paste("kgraph", kgraph.clusters))

ivw <- pheatmap(log10(tab+10),
                main = "NN vs Kgraph",
                color = viridis::viridis(100),
                silent = TRUE)

ivw

## use kgraph, as they have better resolution in the smaller
## peripheral clusters. NN splits the large central cloud in many
## non-distinct clusters.

colLabels(z46fl) <-
    colData(z46fl)$kgraph.clusters <- kgraph.clusters

## --------------------------------------------------------------
## Cell type annotation. Using the first experiment as reference.

sce <- readRDS("sce.rds")
rownames(z46fl) <- rowData(z46fl)$ID
rownames(sce) <- rowData(sce)$gene_ids

sce$labels2 <- as.character(sce$Pheno_name)
i <- which(sce$labels2 == "Macrophage")
sce$labels2[i] <- paste0("Macrophage_",
                         sce$pheno_leiden[i])

if (!file.exists("S04-outputs/pred_full.rds")) {
    pred_full <- SingleR(test=z46fl, ref=sce, labels=sce$labels2)
    saveRDS(pred_full, file = "S04-outputs/pred_full.rds")
} else
    pred_full <- readRDS("S04-outputs/pred_full.rds")
}

z46fl$preds_full <- pred_full$labels

gridExtra::grid.arrange(
               plotTSNE(z46fl[, ctrl],
                        colour_by = "preds_full",
                        text_by = "preds_full") +
               ggtitle("CTRL"),
               plotTSNE(z46fl[, cld],
                        colour_by = "preds_full",
                        text_by = "preds_full") +
               ggtitle("CND"),
               ncol = 1)

tab <- table(colData(z46fl)[, c("preds_full", "Group")])
round(sweep(tab, 2, as.numeric(table(z46fl$Group)), "/"), 3)

## proportions of immune cells

sel <- z46fl$preds_full != "Tumor cells"
tab2 <- table(colData(z46fl)[sel, c("preds_full", "Group")])
tab2 <- round(sweep(tab2, 2, as.numeric(table(z46fl$Group[sel])), "/"), 3)
tab2

## proportions of macrophages

## -------------------------------------------------------------------
## Correspondance between original macrophage subgroups and those in
## here.
##
##     NEW       |    ORIG
## --------------|--------------
## Macrophage_10 | Macrophage_5
## Macrophage_2  | Macrophage_2
## Macrophage_3  | Macrophage_1
## Macrophage_8  | Macrophage_3
## Macrophage_9  | Macrophage_4

sel_macro <- grepl("Macrophage", z46fl$preds_full)
tab_macro <- table(colData(z46fl)[sel_macro, c("preds_full", "Group")])
round(sweep(tab_macro, 2, as.numeric(table(z46fl$Group[sel_macro])), "/"), 3)

write.csv(tab2, file = "S04-outputs/cell_proportions.csv")

## Repeat figures without tumour cells
gridExtra::grid.arrange(
               plotTSNE(z46fl[, sel & ctrl],
                        colour_by = "preds_full",
                        text_by = "preds_full") +
               ggtitle("CTRL"),
               plotTSNE(z46fl[, sel & cld],
                        colour_by = "preds_full",
                        text_by = "preds_full") +
               ggtitle("CND"),
               ncol = 1)
