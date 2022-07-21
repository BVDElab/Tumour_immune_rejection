knitr::opts_chunk$set(echo = TRUE)
library(pander)
library(tidyr)
library(tidyverse)
library(DT)
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
rm(list=ls())

# Data analysis based on `S03.py`

########################
###       DEA        ###
########################
genes_interest <- c("Cd74", "Cxcl2", "H2-Aa", "H2-Ab1", 
                    "H2-D1", "H2-DMa", "H2-DMb1", "H2-Eb1",
                    "H2-K1", "H2-T23")
# merged data
data_path <- paste0("results/merged_filtered/merged_adata.h5ad")

adata_X_merged <- fread(file.path(data_path, "X.csv"), header = F)
adata_var_merged <- fread(file.path(data_path, "var.csv"), header = TRUE)
adata_obs_merged <- fread(file.path(data_path, "obs.csv"), header = TRUE)
range(adata_X_merged[1:100,1:100])

colnames(adata_X_merged) <- adata_var_merged$V1
rownames(adata_X_merged) <- adata_obs_merged$V1
table(adata_obs_merged$treatment)

adata_obs_merged$Pheno_name <- as.factor(adata_obs_merged$Pheno_name)

## volcano plot
##################

### Manually =================


identical(rownames(adata_X_merged),adata_obs_merged$V1)
adata_X_Macrophage <- adata_X_merged[which(adata_obs_merged$Pheno_name=="Macrophage"),]
adata_obs_Macrophage <- adata_obs_merged[which(adata_obs_merged$Pheno_name=="Macrophage"),]

dff <- adata_obs_Macrophage %>% dplyr::select(treatment)
df <- data.frame(adata_X_Macrophage,dff, check.names = FALSE)
res <- df %>% group_by(treatment) %>%
  summarize_all(~mean(.))
res %>% dplyr::select(all_of(c("treatment",genes_interest)))

res2 <- df %>% group_by(treatment) %>%
  summarize_all(~median(.))
res2 %>% dplyr::select(all_of(c("treatment",genes_interest)))

adata_X_Macrophage <- adata_X_merged[which(adata_obs_merged$Pheno_name=="Macrophage"),]
adata_obs_Macrophage <- adata_obs_merged[which(adata_obs_merged$Pheno_name=="Macrophage"),]

apply(adata_X_Macrophage[adata_obs_Macrophage$treatment=="WT_CLD","Cxcl2"],2,mean, na.rm=TRUE)
apply(adata_X_Macrophage[adata_obs_Macrophage$treatment=="WT_Control","Cxcl2"],2,mean, na.rm=TRUE)

apply(adata_X_Macrophage[adata_obs_Macrophage$treatment=="treated","Cxcl2"],2,mean, na.rm=TRUE)
apply(adata_X_Macrophage[adata_obs_Macrophage$treatment=="untreated","Cxcl2"],2,mean, na.rm=TRUE)

plots <- dt_list <- list()

tests <- list("treated_untreated" = c("untreated", "treated"),
              "KO_CLD_Control"=c("KO_Control", "KO_CLD"), "KO_GBZ_Control"=c("KO_Control","KO_GBZ"), 
              "WT_CLD_Control"=c("WT_Control","WT_CLD"), "WT_GBZ_Control"=c("WT_Control","WT_GBZ"))

for (i in 1:length(tests)){
  id0 <- which(adata_obs_Macrophage$treatment==tests[[i]][1])
  id1 <- which(adata_obs_Macrophage$treatment==tests[[i]][2])
  log2fc <- apply(adata_X_Macrophage[id1,],2,mean, na.rm=TRUE) - 
    apply(adata_X_Macrophage[id0,],2,mean, na.rm=TRUE)
  
  pvals <- mapply(adata_X_Macrophage[id0,],adata_X_Macrophage[id1,], 
                  FUN = function(x,y) wilcox.test(x,y, paired=FALSE)$p.value)
  
  qvalue <- p.adjust(pvals, method = "BH")
  df <- data.frame(log2fc, pvals, qvalue) %>% 
    rownames_to_column("gene")
  df$genes_of_interest <- df$gene %in% genes_interest
  
  decideTests <- rep(0, nrow(df))
  decideTests[df$qvalue<0.05 & abs(df$log2fc)<0.5] = -1
  decideTests[df$qvalue<0.05 & abs(df$log2fc)>0.5] = 1
  df$decideTests <- decideTests
  
  dt_list[[names(tests)[i]]]  <- df
  
}


sqval <- sort(dt_list[["KO_CLD_Control"]]$qvalue)
id <- order(dt_list[["KO_CLD_Control"]]$qvalue)[1]
dt_list[["KO_CLD_Control"]]$qvalue[id] <- sqval[2]


for (i in 1:length(tests)){
  
  df <- dt_list[[names(tests)[i]]]
  plots[[names(tests)[i]]] <- 
    ggplot(df, aes(x = log2fc, y = -log10(qvalue)))  +
    geom_hline(yintercept = -log10(0.05), color="cornflowerblue",
               linetype = "dashed")+
    geom_vline(xintercept = c(-1,1), color="chartreuse3",
               linetype = "dashed")+
    geom_point(aes(color=genes_of_interest), shape=16, alpha=0.7) +
    geom_point(data=df[df$gene %in%genes_interest,], 
               fill="red", shape=21) +
    scale_color_manual(values = c("grey", "red")) + 
    geom_text_repel(data=df[df$gene %in%genes_interest,], 
                    aes(label = gene), max.overlaps = 100)+
    ggtitle(names(tests)[i]) + theme_minimal() + xlim(-max(abs(df$log2fc)),
                                                      max(abs(df$log2fc))) +
    theme(legend.position = "none")
  
}

plots



### Based on diffxpy =================

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
  list_dt[[i]] <- read.table(file = file.path(data_path,tsv_files[i]), 
                             sep = '\t', header = TRUE)

}
names(list_dt) <- test

test_title[test_title=="treated_NA vs treated_untreated"] = "treated vs untreated"

## volcano plot 
plots <- list()
for (i in 1:length(tsv_files)){
  dt <- list_dt[[i]]
  sqval <- sort(dt$qval)
  id <- order(dt$qval)[1]
  dt$qval[id] <- sqval[2]
  
  dt$genes_of_interest <- dt$gene %in% genes_interest
  decideTests <- rep(0, nrow(dt))
  decideTests[dt$qval<0.05 & dt$log2fc<0] = -1
  decideTests[dt$qval<0.05 & dt$log2fc>0] = 1
  dt$decideTests <- decideTests
  list_dt[[i]]$decideTests <- decideTests
  
  
  plots[[i]] <- ggplot(dt, aes(x = log2fc, y = -log10(qval)))  +
    geom_hline(yintercept = -log10(0.05), color="cornflowerblue",
               linetype = "dashed")+
    geom_vline(xintercept = c(-1,1), color="chartreuse3",
               linetype = "dashed")+
    geom_point(aes(color=genes_of_interest), shape=16, alpha=0.7) +
    geom_point(data=dt[dt$gene %in%genes_interest,], 
               fill="red", shape=21) +
    scale_color_manual(values = c("grey", "red")) + 
    geom_text_repel(data=dt[dt$gene %in%genes_interest,], 
                    aes(label = gene), max.overlaps = 100)+
    ggtitle(test_title[i]) + theme_minimal() + xlim(-max(abs(dt$log2fc)),
                                                    max(abs(dt$log2fc))) +
    theme(legend.position = "none")
  
}

plots

for (i in 1:length(tsv_files)){
  svg(paste0("figures/merged_filtered/volcanoPlot",test_title[i],
             ".svg"), height = 5, width = 5,  family = "sans")
  print(plots[[i]])
  dev.off()
}




## upset plot 
##################

list_all_genes <- lapply(list_dt, function(x) x$gene)
common_genes <- Reduce(intersect, list_all_genes)

# Look where the fold changes are significantly higher (=1)
pander("Look where the fold changes are significantly higher (=1)")

list_genes <- lapply(list_dt, function(x) x$gene[x$decideTests == 1])

for (i in 1:length(list_genes)){
  list_genes[[i]] <- list_genes[[i]][list_genes[[i]] %in% common_genes]
}

names(list_genes) <- test_title

us <- fromList(list_genes)
rownames(us) <- unique(unlist(list_genes))

UpSetR::upset(us,nsets = 5, keep.order = TRUE,sets = names(list_genes))
grid.text("logFC + & q-value<0.05 ",x = 0.75, 
          y=0.95, gp=gpar(fontsize=13))

svg(paste0("figures/merged_filtered/Upsetplot_up.svg"), height = 5, width = 7,  family = "sans")
UpSetR::upset(us,nsets = 5, keep.order = TRUE,sets = names(list_genes))
grid.text("logFC + & q-value<0.05 ",x = 0.75, 
          y=0.95, gp=gpar(fontsize=13))
dev.off()

# Look where the fold changes are significantly lower (=-1)
pander("Look where the fold changes are significantly lower (=-1)")

list_genes <- lapply(list_dt, function(x) x$gene[x$decideTests == -1])

for (i in 1:length(list_genes)){
  list_genes[[i]] <- list_genes[[i]][list_genes[[i]] %in% common_genes]
}

names(list_genes) <- test_title

us <- fromList(list_genes)
rownames(us) <- unique(unlist(list_genes))

UpSetR::upset(us,nsets = 5, keep.order = TRUE,sets = names(list_genes))
grid.text("logFC - & q-value<0.05 ",x = 0.75, 
          y=0.95, gp=gpar(fontsize=13))

svg(paste0("figures/merged_filtered/Upsetplot_down.svg"), height = 6, width = 7,  family = "sans")
UpSetR::upset(us,nsets = 5, keep.order = TRUE,sets = names(list_genes))
grid.text("logFC - & q-value<0.05 ",x = 0.75, 
          y=0.95, gp=gpar(fontsize=13))
dev.off()


##################
##  violin plot ##
##################


adata_X_Macrophage <- adata_X_merged[which(adata_obs_merged$Pheno_name=="Macrophage"),]
adata_obs_Macrophage <- adata_obs_merged[which(adata_obs_merged$Pheno_name=="Macrophage"),]

df <- data.frame(adata_X_Macrophage,Pheno_name=adata_obs_Macrophage$Pheno_name) 

##########################
##    Violins           ## 
##########################

genes_interest <- c("Cd74", "Cxcl2", "H2-Aa", "H2-Ab1", 
                    "H2-D1", "H2-DMa", "H2-DMb1", "H2-Eb1",
                    "H2-K1", "H2-T23")
id <- which(colnames(adata_X_merged) %in% genes_interest)
df <- adata_X_merged[,..id]
df$rowname <- rownames(adata_X_merged)
dff <- adata_obs_merged %>% dplyr::select(c(V1, treatment))
df <- df %>% left_join(dff, by=c("rowname"="V1"))

df <- df %>% pivot_longer(!c(rowname,treatment))
plot <- df %>% ggplot(aes(x = treatment, y = value)) +
  geom_violin() + facet_wrap(~name, scales = "free")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot

svg("figures/merged_filtered/violin_plot.svg", height = 8, width = 10,  family = "sans")
plot
dev.off()

