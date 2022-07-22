library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
require(data.table)
library(RColorBrewer)
library(ggrepel)
library(UpSetR)
library(MASS)

# Data analysis based on `S03.py`

genes_interest <- c("Cd74", "Cxcl2", "H2-Aa", "H2-Ab1", 
                    "H2-D1", "H2-DMa", "H2-DMb1", "H2-Eb1",
                    "H2-K1", "H2-T23")
# merged data
data_path <- paste0("results/merged_filtered/merged_adata.h5ad")

adata_X_merged <- fread(file.path(data_path, "X.csv"), header = FALSE)
adata_var_merged <- fread(file.path(data_path, "var.csv"), header = TRUE)
adata_obs_merged <- fread(file.path(data_path, "obs.csv"), header = TRUE)
range(adata_X_merged[1:100,1:100])

colnames(adata_X_merged) <- adata_var_merged$V1
rownames(adata_X_merged) <- adata_obs_merged$V1
table(adata_obs_merged$treatment)

adata_obs_merged$Pheno_name <- as.factor(adata_obs_merged$Pheno_name)

library(SingleCellExperiment)
sce = SingleCellExperiment(assays = list(logcounts = t(adata_X_merged)))
colData(sce) = DataFrame(adata_obs_merged)
rowData(sce) = DataFrame(adata_var_merged)
colnames(sce) = colData(sce)$V1

m = sce[, sce$Pheno_name == "Macrophage"]

#######################
## Calculate results
#######################

tests <- list("treated_untreated" = c("untreated", "treated"),
              "KO_CLD_Control"=c("KO_Control", "KO_CLD"),
              "KO_GBZ_Control"=c("KO_Control", "KO_GBZ"), 
              "WT_CLD_Control"=c("WT_Control", "WT_CLD"),
              "WT_GBZ_Control"=c("WT_Control", "WT_GBZ"))

## res = vector("list", length = length(tests))
## names(res) = names(tests)
## for (i in 1:length(tests)) {
##     gr1 = which(m$treatment == tests[[i]][1])
##     gr2 = which(m$treatment == tests[[i]][2])
    
##     lfc = rowMeans(assay(m)[, gr2]) - rowMeans(assay(m)[, gr1])
##     pvals <- mapply(assay(m)[, gr2],
##                     assay(m)[, gr1],
##                     FUN = function(x,y) wilcox.test(x,y, paired = FALSE)$p.value)  
##     qvals <- p.adjust(pvals, method = "BH")
##     df <- data.frame(lfc, pvals, qvals) %>% 
##         rownames_to_column("gene")
##     df$genes_of_interest <- df$gene %in% genes_interest    
##     decideTests <- rep(0, nrow(df))
##     decideTests[df$qvals < 0.05 & df$lfc < 0.5] = -1
##     decideTests[df$qvals < 0.05 & df$lfc > 0.5] = 1
##     df$decideTests <- decideTests
##     res[[i]] = df  
## }


########################
## Load diffxpy results
########################

tsv_files <- list.files("results/merged_filtered/Macrophage/", )
data_path <- "results/merged_filtered/Macrophage/"

list_dt <- list()
test_title <- test <- c()
for (i in 1:length(tsv_files)){
  test[i] <- sub("_Diff.*","",tsv_files[i])
  test[i] <- sub("Macrophage_","",test[i])
  strsplit_res <- strsplit(test[i],"_")[[1]]
  test_title[i] <- paste0(strsplit_res[1],"_",strsplit_res[3],
                          " vs ",strsplit_res[1],"_",strsplit_res[2],
                          collapse = "")
  df = read.table(file = file.path(data_path,tsv_files[i]), 
                  sep = '\t', header = TRUE)
  df$genes_of_interest <- df$gene %in% genes_interest    
  decideTests <- rep(0, nrow(df))
  decideTests[df$qval < 0.05 & df$log2fc < -0.5] = -1
  decideTests[df$qval < 0.05 & df$log2fc > 0.5] = 1
  df$decideTests <- decideTests  
  list_dt[[i]] <- df
}
names(list_dt) <- test
test_title[test_title=="treated_NA vs treated_untreated"] = "treated vs untreated"
resdx = list_dt

for (i in 1:length(resdx)) {
    write.csv(resdx[[i]][, -1], file = paste0("stats-", names(resdx)[i], ".csv"),
              row.names = FALSE)
}

##################
## volcano plot
##################

## for (i in 1:length(tests)){  
##   df <- dt_list[[names(tests)[i]]]
##   plots[[names(tests)[i]]] <- 
##     ggplot(df, aes(x = log2fc, y = -log10(qvalue)))  +
##     geom_hline(yintercept = -log10(0.05), color="cornflowerblue",
##                linetype = "dashed")+
##     geom_vline(xintercept = c(-1,1), color="chartreuse3",
##                linetype = "dashed")+
##     geom_point(aes(color=genes_of_interest), shape=16, alpha=0.7) +
##     geom_point(data=df[df$gene %in%genes_interest,], 
##                fill="red", shape=21) +
##     scale_color_manual(values = c("grey", "red")) + 
##     geom_text_repel(data=df[df$gene %in%genes_interest,], 
##                     aes(label = gene), max.overlaps = 100)+
##     ggtitle(names(tests)[i]) + theme_minimal() + xlim(-max(abs(df$log2fc)),
##                                                       max(abs(df$log2fc))) +
##     theme(legend.position = "none")
## }


## plots <- list()
## for (i in 1:length(tsv_files)){
##   dt <- list_dt[[i]]
##   sqval <- sort(dt$qval)
##   id <- order(dt$qval)[1]
##   dt$qval[id] <- sqval[2]
##   dt$genes_of_interest <- dt$gene %in% genes_interest
##   plots[[i]] <- ggplot(dt, aes(x = log2fc, y = -log10(qval)))  +
##     geom_hline(yintercept = -log10(0.05), color="cornflowerblue",
##                linetype = "dashed")+
##     geom_vline(xintercept = c(-1,1), color="chartreuse3",
##                linetype = "dashed")+
##     geom_point(aes(color=genes_of_interest), shape=16, alpha=0.7) +
##     geom_point(data=dt[dt$gene %in%genes_interest,], 
##                fill="red", shape=21) +
##     scale_color_manual(values = c("grey", "red")) + 
##     geom_text_repel(data=dt[dt$gene %in%genes_interest,], 
##                     aes(label = gene), max.overlaps = 100)+
##     ggtitle(test_title[i]) + theme_minimal() + xlim(-max(abs(dt$log2fc)),
##                                                     max(abs(dt$log2fc))) +
##     theme(legend.position = "none")
## }

## plots

## for (i in 1:length(tsv_files)){
##   svg(paste0("figures/merged_filtered/volcanoPlot",test_title[i],
##              ".svg"), height = 5, width = 5,  family = "sans")
##   print(plots[[i]])
##   dev.off()
## }


##################
## UpSet plot 
##################

list_genes_up <- lapply(list_dt, function(x) x$gene[x$decideTests == 1])
list_genes_down <- lapply(list_dt, function(x) x$gene[x$decideTests == -1])
names(list_genes_down) <- names(list_genes_up) <- names(resdx)

up <- fromList(list_genes_up)
down <- fromList(list_genes_down)

pdf("upset-up.pdf")
UpSetR::upset(up, nsets = 5, keep.order = TRUE, sets = names(list_genes_up))
grid.text("logFC > 0.5 & q-value < 0.05 ",
          x = 0.75, 
          y = 0.95, gp=gpar(fontsize=13))
dev.off()

pdf("upset-down.pdf")
UpSetR::upset(down, nsets = 5, keep.order = TRUE, sets = names(list_genes_down))
grid.text("logFC < -0.5 & q-value < 0.05 ",
          x = 0.75, 
          y = 0.95, gp=gpar(fontsize=13))
dev.off()

#######################
##  lfc correlations ##
#######################

lfcmat = Reduce(cbind, lapply(resdx, "[[", "log2fc"))
colnames(lfcmat) = names(resdx)
rownames(lfcmat) = resdx[[1]]$gene

lfcmat[abs(lfcmat) < 0.5] = NA

write.csv(cor(lfcmat, use = "pairwise.complete.obs"), file = "cor-lfc.csv")

## |                  | KO_Control_CLD| KO_Control_GBZ| treated_untreated| WT_Control_CLD| WT_Control_GBZ|
## |:-----------------|--------------:|--------------:|-----------------:|--------------:|--------------:|
## |KO_Control_CLD    |          1.000|          0.981|             0.260|         -0.421|          0.902|
## |KO_Control_GBZ    |          0.981|          1.000|            -0.100|          0.596|          0.991|
## |treated_untreated |          0.260|         -0.100|             1.000|         -0.803|         -0.230|
## |WT_Control_CLD    |         -0.421|          0.596|            -0.803|          1.000|         -0.177|
## |WT_Control_GBZ    |          0.902|          0.991|            -0.230|         -0.177|          1.000|

colnames(m) = sub("-", "_", colnames(m))
m$V1 = gsub("-", "_", m$V1)



pdf("violins.pdf", width = 16, height = 8)
data.frame(assay(m[genes_interest, ])) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    pivot_longer(values_to = "exprs",
                 names_to = "V1",
                 -gene) %>%
    full_join(data.frame(colData(m)[, c("V1", "treatment")])) %>%
    ggplot(aes(x = treatment, y = exprs)) +
    geom_violin() +
    facet_grid(~ gene) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()    


pdf("violins-no-gbz.pdf", width = 16, height = 8)
data.frame(assay(m[genes_interest, ])) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    pivot_longer(values_to = "exprs",
                 names_to = "V1",
                 -gene) %>%
    full_join(data.frame(colData(m)[, c("V1", "treatment")])) %>%
    filter(treatment %in% c("WT_CLD", "WT_Control", "KO_CLD", "KO_Control")) %>%
    ggplot(aes(x = treatment, y = exprs)) +
    geom_violin() +
    facet_grid(~ gene) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()    

