### run scd_220628 and scd_d660cc data
########################################

import scanpy as sc # https://scanpy.readthedocs.io/en/stable/usage-principles.html
import pandas as pd
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy # https://gseapy.readthedocs.io/en/latest/introduction.html
import random
from collections import Counter
import diffxpy.api as de # https://github.com/theislab/diffxpy
import os
from gseapy.plot import gseaplot
# other modules that are necessary: openpyxl, PhenoGraph

plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['svg.fonttype'] = 'none'

# which data should we use?
tp = "filtered"
expe = "merged"

if not os.path.isdir("figures/"+ expe+"_"+tp):
    os.mkdir("figures/"+ expe+"_"+tp)
    
if not os.path.isdir("results/"+ expe+"_"+tp):
    os.mkdir("results/"+ expe+"_"+tp)

# Make sure environment is always the same
random.seed(770)
np.random.seed(770)

###############
#   Load data
###############

### design ### 
# untreated    Control
# treated      CLD
# Z1	    KO Control
# Z2	    KO CLD
# Z3	    KO GBZ
# Z4	    WT Control
# Z5	    WT GBZ
# Z6	    WT CLD


print("############### Load data ###############")
### untreated
untreated = sc.read_10x_mtx("data/" + "scd_d660cc"+"_" + tp + "/untreated", prefix = "DDI-JZ-g001-"+tp+"-")
type(untreated)
sc.pl.highest_expr_genes(untreated, n_top=20)
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/Top20-most-expressed-genes_untreated.png", dpi=300)
plt.close()

### treated
treated = sc.read_10x_mtx("data/" + "scd_d660cc"+"_" + tp + "/treated", prefix = "DDI-JZ-g002-"+tp+"-")
sc.pl.highest_expr_genes(treated, n_top=20)
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/Top20-most-expressed-genes_treated.png", dpi=300)
plt.close()


### Z1
Z1 = sc.read_10x_mtx("data/" + "scd_220628"+"_" + tp + "/Z1")
type(Z1)
sc.pl.highest_expr_genes(Z1, n_top=20)
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/Top20-most-expressed-genes_Z1.png", dpi=300)
plt.close()

### Z2
Z2 = sc.read_10x_mtx("data/" + "scd_220628"+"_" + tp + "/Z2")
type(Z2)
sc.pl.highest_expr_genes(Z2, n_top=20)
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/Top20-most-expressed-genes_Z2.png", dpi=300)
plt.close()

### Z3
Z3 = sc.read_10x_mtx("data/" + "scd_220628"+"_" + tp + "/Z3")
type(Z3)
sc.pl.highest_expr_genes(Z3, n_top=20)
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/Top20-most-expressed-genes_Z3.png", dpi=300)
plt.close()

### Z4
Z4 = sc.read_10x_mtx("data/" + "scd_220628"+"_" + tp + "/Z4")
type(Z4)
sc.pl.highest_expr_genes(Z4, n_top=20)
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/Top20-most-expressed-genes_Z4.png", dpi=300)
plt.close()

### Z5
Z5 = sc.read_10x_mtx("data/" + "scd_220628"+"_" + tp + "/Z5")
type(Z5)
sc.pl.highest_expr_genes(Z5, n_top=20)
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/Top20-most-expressed-genes_Z5.png", dpi=300)
plt.close()

### Z6
Z6 = sc.read_10x_mtx("data/" + "scd_220628"+"_" + tp + "/Z6")
type(Z6)
sc.pl.highest_expr_genes(Z6, n_top=20)
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/Top20-most-expressed-genes_Z6.png", dpi=300)
plt.close()

# Mark mitochondrial content
untreated.var['mt'] = untreated.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
treated.var['mt'] = treated.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
Z1.var['mt'] = Z1.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
Z2.var['mt'] = Z2.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
Z3.var['mt'] = Z3.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
Z4.var['mt'] = Z4.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
Z5.var['mt'] = Z5.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
Z6.var['mt'] = Z6.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'

###############
#    QC
###############
print("############### QC ###############")
# Largely based on calculateQCMetrics from scater
sc.pp.calculate_qc_metrics(untreated, qc_vars=['mt'], inplace=True)
sc.pp.calculate_qc_metrics(treated, qc_vars=['mt'], inplace=True)
sc.pp.calculate_qc_metrics(Z1, qc_vars=['mt'], inplace=True)
sc.pp.calculate_qc_metrics(Z2, qc_vars=['mt'], inplace=True)
sc.pp.calculate_qc_metrics(Z3, qc_vars=['mt'], inplace=True)
sc.pp.calculate_qc_metrics(Z4, qc_vars=['mt'], inplace=True)
sc.pp.calculate_qc_metrics(Z5, qc_vars=['mt'], inplace=True)
sc.pp.calculate_qc_metrics(Z6, qc_vars=['mt'], inplace=True)

### violin
# untreated
sc.pl.violin(untreated, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("figures/"+expe+"_"+tp+"/untreated_qc.png", dpi=300)
plt.close()
print("untreated", "n_genes_by_counts", np.mean(untreated.obs["n_genes_by_counts"]))
print("untreated", "total_counts", np.mean(untreated.obs["total_counts"]))
print("untreated", "pct_counts_mt", np.mean(untreated.obs["pct_counts_mt"]))

# treated
sc.pl.violin(treated, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("figures/"+expe+"_"+tp+"/untreated_qc.png", dpi=300)
plt.close()
print("treated", "n_genes_by_counts", np.mean(treated.obs["n_genes_by_counts"]))
print("treated", "total_counts", np.mean(treated.obs["total_counts"]))
print("treated", "pct_counts_mt", np.mean(treated.obs["pct_counts_mt"]))

# Z1
sc.pl.violin(Z1, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("figures/"+expe+"_"+tp+"/Z1_qc.png", dpi=300)
plt.close()
print("Z1", "n_genes_by_counts", np.mean(Z1.obs["n_genes_by_counts"]))
print("Z1", "total_counts", np.mean(Z1.obs["total_counts"]))
print("Z1", "pct_counts_mt", np.mean(Z1.obs["pct_counts_mt"]))

# Z2
sc.pl.violin(Z2, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("figures/"+expe+"_"+tp+"/Z2_qc.png", dpi=300)
plt.close()
print("Z2", "n_genes_by_counts", np.mean(Z2.obs["n_genes_by_counts"]))
print("Z2", "total_counts", np.mean(Z2.obs["total_counts"]))
print("Z2", "pct_counts_mt", np.mean(Z2.obs["pct_counts_mt"]))


# Z3
sc.pl.violin(Z3, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("figures/"+expe+"_"+tp+"/Z3_qc.png", dpi=300)
plt.close()
print("Z3", "n_genes_by_counts", np.mean(Z3.obs["n_genes_by_counts"]))
print("Z3", "total_counts", np.mean(Z3.obs["total_counts"]))
print("Z3", "pct_counts_mt", np.mean(Z3.obs["pct_counts_mt"]))

# Z4
sc.pl.violin(Z4, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("figures/"+expe+"_"+tp+"/Z4_qc.png", dpi=300)
plt.close()
print("Z4", "n_genes_by_counts", np.mean(Z4.obs["n_genes_by_counts"]))
print("Z4", "total_counts", np.mean(Z4.obs["total_counts"]))
print("Z4", "pct_counts_mt", np.mean(Z4.obs["pct_counts_mt"]))

# Z5
sc.pl.violin(Z5, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("figures/"+expe+"_"+tp+"/Z5_qc.png", dpi=300)
plt.close()
print("Z5", "n_genes_by_counts", np.mean(Z5.obs["n_genes_by_counts"]))
print("Z5", "total_counts", np.mean(Z5.obs["total_counts"]))
print("Z5", "pct_counts_mt", np.mean(Z5.obs["pct_counts_mt"]))

# Z6
sc.pl.violin(Z1, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("figures/"+expe+"_"+tp+"/Z6_qc.png", dpi=300)
plt.close()
print("Z6", "n_genes_by_counts", np.mean(Z6.obs["n_genes_by_counts"]))
print("Z6", "total_counts", np.mean(Z6.obs["total_counts"]))
print("Z6", "pct_counts_mt", np.mean(Z6.obs["pct_counts_mt"]))


### scatter
# untreated
sc.pl.scatter(untreated, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/"+expe+"_"+tp+"/totalcounts-ngenes_untreated.png", dpi=300)
plt.close()
# treated
sc.pl.scatter(treated, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/"+expe+"_"+tp+"/totalcounts-ngenes_treated.png", dpi=300)
plt.close()
# Z1
sc.pl.scatter(Z1, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/"+expe+"_"+tp+"/totalcounts-ngenes_Z1.png", dpi=300)
plt.close()
# Z2
sc.pl.scatter(Z2, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/"+expe+"_"+tp+"/totalcounts-ngenes_Z2.png", dpi=300)
plt.close()
# Z3
sc.pl.scatter(Z3, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/"+expe+"_"+tp+"/totalcounts-ngenes_Z3.png", dpi=300)
plt.close()
# Z4
sc.pl.scatter(Z4, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/"+expe+"_"+tp+"/totalcounts-ngenes_Z4.png", dpi=300)
plt.close()
# Z5
sc.pl.scatter(Z5, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/"+expe+"_"+tp+"/totalcounts-ngenes_Z5.png", dpi=300)
plt.close()
# Z6
sc.pl.scatter(Z6, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/"+expe+"_"+tp+"/totalcounts-ngenes_Z6.png", dpi=300)
plt.close()


#########################
#  Concatenate the data
#########################
print("############### Concatenate ###############")
adata = untreated.concatenate(treated,Z1,Z2,Z3,Z4,Z5,Z6)
print(adata)

adata = ad.AnnData(X = adata.X, obs = adata.obs, var = adata.var[["gene_ids","mt"]].copy())


##############################
#  Filter genes and cells out
##############################
print("############### Filter genes and cells out ###############")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
sc.pp.filter_cells(adata, min_genes=200) # Filter cell outliers based on counts and numbers of genes expressed
sc.pp.filter_genes(adata, min_cells=10) # Filter genes based on number of cells or counts
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
print(adata)

adata.write_csvs("results/"+ expe+"_"+tp+"/"+expe+"_adata_raw.h5ad",skip_data=False)


##############################
# Normalize per 10k reads + logtransform, reduce variability
##############################
print("############### Normalize, logtransform, reduce variability ###############")
sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=False, inplace=True) # every cell has the same total count after normalization
sc.pp.log1p(adata, base=2)
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt']) # inspired by Seurat’s regressOut function in R
adata = adata[:, [x for x in adata.var_names if not x.startswith("mt-")]].copy()

##############################
# Variable gene detection
##############################
print("############### Variable gene detection ###############")
adata.obs["treatment"] = adata.obs["batch"].map({"0":"untreated","1":"treated",
"2": "KO_Control", "3": "KO_CLD",
"4": "KO_GBZ", "5": "WT_Control", "6": "WT_GBZ", "7": "WT_CLD"})

sc.pp.highly_variable_genes(adata, flavor="seurat")
sc.pl.highly_variable_genes(adata)
plt.savefig("figures/"+expe+"_"+tp+"/highly_variable_prop.png", dpi=300)
plt.close()


##############################
#    PCA
##############################
print("############### PCA ###############")
highly_variable = adata.var.loc[adata.var["highly_variable"]].index
print("# Number of highly variable genes")
print(len(highly_variable))
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True, random_state=770)
adata.obsm['X_pca']
adata.uns['pca']['variance_ratio']
sc.pl.pca(adata, color=["treatment"], cmap="YlOrRd")
plt.tight_layout()
sns.despine()
plt.savefig("figures/"+expe+"_"+tp+"/batch_correction.png", dpi=300)
plt.close()

# Cumulative variance explained:
cml_var_explained = np.cumsum(adata.uns['pca']['variance_ratio'])
x = range(len(adata.uns['pca']['variance_ratio']))
y = cml_var_explained
plt.scatter(x, y, s=4)
plt.xlabel('PC')
plt.ylabel('Cumulative variance explained')
plt.title('Cumulative variance explained by PCs')
plt.savefig("figures/"+expe+"_"+tp+"/components_needed_pervariance.png", dpi=300)
plt.close()

##############################
# data projection on highly variable genes
##############################
print("############### data projection on highly variable genes ###############")
subdata = adata[:, highly_variable].copy()
# PhenoGraph clustering
print("############### PhenoGraph clustering ###############")
sc.external.tl.phenograph(subdata, clustering_algo="leiden", seed=770, k=25, n_iterations=5000)
adata.obs["pheno_leiden"] = subdata.obs["pheno_leiden"]
sc.tl.rank_genes_groups(adata, 'pheno_leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("figures/"+expe+"_"+tp+"/rank_genes_groups.png", dpi=300)
plt.close()


# Data projections
sc.tl.pca(adata, n_comps=26, use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=25, n_pcs=26)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["pheno_leiden"], legend_loc="on data")
sns.despine()
plt.tight_layout()
plt.savefig("figures/"+expe+"_"+tp+"/leiden_umap.png", dpi=300)
plt.close()

##############################
# umap with identified mappings
##############################
print("############### umap with identified mappings ###############")
# Mappings we identified
mapdic = {0: "Mfab5- fibroblast", 1: "Mfab5- fibroblast", 
          2: "Macrophage", 3: "Macrophage",
          4: "duplex", 5: "Tumor cells", 6: "Tumor cells", 7: "Erythrocytes",
          8: "Macrophage", 9: "Macrophage", 10: "Macrophage", 
          11: "CD8",
          12: "Mfab5+ fibroblast", 13: "CD4 FoxP3+",14: "Erythrocytes",
          15: "Mfab5- fibroblast",16: "Tumor cells", 17: "Mfab5+ fibroblast",
          18: "CD4 FoxP3-", 19: "Neutrophil",
          20: "Dendritic cells", 21: "Mast cells",
          22: "Endothelial cells"}

adata.obs["Pheno_name"] = adata.obs["pheno_leiden"].replace(mapdic)
adata.obs["Pheno_name"].fillna("Unknown", inplace=True)

adata.obs.to_csv(os.path.join("results", expe+"_"+tp,expe+"_adata_obs.tsv"), sep="\t")
adata.write_csvs("results/"+ expe+"_"+tp+"/"+expe+"_adata.h5ad",skip_data=False)

sc.pl.umap(adata, color=["Pheno_name"], legend_loc="on data")
plt.tight_layout()
sns.despine()
plt.savefig("figures/"+expe+"_"+tp+"/"+expe+"_UMAP_celltype.png", dpi=300)
plt.close()

for tment in set(adata.obs["treatment"]):
  tindexes = adata.obs.loc[adata.obs["treatment"] == tment].index
  subdata = adata[tindexes, :].copy()
  sc.pl.umap(subdata, color="Pheno_name", legend_loc="on data", legend_fontweight="normal")
  # sc.pl.umap(subdata, color="Pheno_name", legend_loc="on data", palette=mycmap, legend_fontweight="normal")
  plt.tight_layout()
  sns.despine()
  plt.savefig("figures/"+expe+"_"+tp+"/"+expe+"_UMAP_" + tment + "_celltype.png", dpi=300)
  plt.close()

##############################
# Piecharts
##############################
print("############### Piecharts ###############")
countlist = []
for agroup in set(adata.obs["Pheno_name"]):
    cdic = Counter(adata.obs["treatment"].loc[adata.obs["Pheno_name"] == agroup])
    countlist.append([agroup, "untreated", cdic["untreated"]])
    countlist.append([agroup, "treated", cdic["treated"]])
countlist = pd.DataFrame(countlist, columns=["Cell_type", "Treatment", "Count"])
countlist.to_excel("results/"+ expe+"_"+tp+"/"+expe+"_Cell_Counts.xlsx")


##############################
# Differential testing
##############################
# print("############### Differential testing ###############")
# for x in set(adata.obs["Pheno_name"]):
# if not os.path.exists(os.path.join("results", expe+"_"+tp, "Macrophage")):
#     os.mkdir(os.path.join("results",expe+"_"+tp, "Macrophage"))
# 
# for acol in ["Macrophage"]:
#     subindexes = list(adata.obs.loc[adata.obs["Pheno_name"] == acol].index)
#     subdata = adata[subindexes].copy()
#     
# 
#     # treated vs untreated ==============================
#     cond = "treated_untreated"
#     subindexes = list(subdata.obs.loc[np.logical_or(subdata.obs["treatment"] == 'untreated', 
#     subdata.obs["treatment"] == 'treated')].index)
#     subdata1 = subdata[subindexes].copy()
#     subdata1.obs["treatment"] = subdata1.obs["treatment"].cat.remove_unused_categories()
#     new_labels = {'untreated' : 'a_untreated', 'treated' : 'b_treated'}
#     subdata1.obs['treatment'] = subdata1.obs['treatment'].cat.rename_categories(new_labels)
#     np.unique(subdata1.obs["treatment"]) 
#     
#     test = de.test.rank_test(data=subdata1, grouping="treatment", is_logged=True)
#     test.summary()["log2fc"]
#     
#     tester = test.summary().loc[
#         test.summary()["gene"].isin([x for x in test.summary()["gene"] if not x.startswith("mt")])]
# 
#     tester.to_csv(os.path.join("results",expe+"_"+tp,acol, acol+"_"+cond+"_"+"Diffxpy_Statistics.tsv"), sep="\t")
# 
#     # KO_Control vs KO_CLD ==============================
#     cond = "KO_Control_CLD"
#     subindexes = list(subdata.obs.loc[np.logical_or(subdata.obs["treatment"] == 'KO_Control', 
#     subdata.obs["treatment"] == 'KO_CLD')].index)
#     subdata2 = subdata[subindexes].copy()
#     subdata2.obs["treatment"] = subdata2.obs["treatment"].cat.remove_unused_categories()
#     new_labels = {'KO_Control' : 'a_KO_Control', 'KO_CLD' : 'b_KO_CLD'}
#     subdata2.obs['treatment'] = subdata2.obs['treatment'].cat.rename_categories(new_labels)
#     np.unique(subdata2.obs["treatment"]) 
#     
#     test = de.test.rank_test(data=subdata2, grouping="treatment", is_logged=True)
#     test.summary()["log2fc"]
#     
#     tester = test.summary().loc[
#         test.summary()["gene"].isin([x for x in test.summary()["gene"] if not x.startswith("mt")])]
# 
#     tester.to_csv(os.path.join("results",expe+"_"+tp,acol, acol+"_"+cond+"_"+"Diffxpy_Statistics.tsv"), sep="\t")
# 
#     
#     # KO_Control vs KO_GBZ ==============================
#     cond = "KO_Control_GBZ"
#     subindexes = list(subdata.obs.loc[np.logical_or(subdata.obs["treatment"] == 'KO_Control', 
#     subdata.obs["treatment"] == 'KO_GBZ')].index)
#     subdata3 = subdata[subindexes].copy()
#     subdata3.obs["treatment"] = subdata3.obs["treatment"].cat.remove_unused_categories()
#     new_labels = {'KO_Control' : 'a_KO_Control', 'KO_GBZ' : 'b_KO_GBZ'}
#     subdata3.obs['treatment'] = subdata3.obs['treatment'].cat.rename_categories(new_labels)
#     np.unique(subdata3.obs["treatment"]) 
#     
#     test = de.test.rank_test(data=subdata3, grouping="treatment", is_logged=True)
#     test.summary()["log2fc"]
#     
#     tester = test.summary().loc[
#         test.summary()["gene"].isin([x for x in test.summary()["gene"] if not x.startswith("mt")])]
# 
#     tester.to_csv(os.path.join("results",expe+"_"+tp,acol, acol+"_"+cond+"_"+"Diffxpy_Statistics.tsv"), sep="\t")
# 
#     # WT_Control vs WT_CLD ==============================
#     cond = "WT_Control_CLD"
#     subindexes = list(subdata.obs.loc[np.logical_or(subdata.obs["treatment"] == 'WT_Control', 
#     subdata.obs["treatment"] == 'WT_CLD')].index)
#     subdata4 = subdata[subindexes].copy()
#     subdata4.obs["treatment"] = subdata4.obs["treatment"].cat.remove_unused_categories()
#     new_labels = {'WT_Control' : 'a_WT_Control', 'WT_CLD' : 'b_WT_CLD'}
#     subdata4.obs['treatment'] = subdata4.obs['treatment'].cat.rename_categories(new_labels)
#     np.unique(subdata4.obs["treatment"]) 
#     
#     test = de.test.rank_test(data=subdata4, grouping="treatment", is_logged=True)
#     test.summary()["log2fc"]
#     
#     tester = test.summary().loc[
#         test.summary()["gene"].isin([x for x in test.summary()["gene"] if not x.startswith("mt")])]
# 
#     tester.to_csv(os.path.join("results",expe+"_"+tp,acol, acol+"_"+cond+"_"+"Diffxpy_Statistics.tsv"), sep="\t")
# 
#     # WT_Control vs WT_GBZ ==============================
#     cond = "WT_Control_GBZ"
#     subindexes = list(subdata.obs.loc[np.logical_or(subdata.obs["treatment"] == 'WT_Control', 
#     subdata.obs["treatment"] == 'WT_GBZ')].index)
#     subdata5 = subdata[subindexes].copy()
#     subdata5.obs["treatment"] = subdata5.obs["treatment"].cat.remove_unused_categories()
#     new_labels = {'WT_Control' : 'a_WT_Control', 'WT_GBZ' : 'b_WT_GBZ'}
#     subdata5.obs['treatment'] = subdata5.obs['treatment'].cat.rename_categories(new_labels)
#     np.unique(subdata5.obs["treatment"]) 
#     
#     test = de.test.rank_test(data=subdata5, grouping="treatment", is_logged=True)
#     test.summary()["log2fc"]
#     
#     tester = test.summary().loc[
#         test.summary()["gene"].isin([x for x in test.summary()["gene"] if not x.startswith("mt")])]
# 
#     tester.to_csv(os.path.join("results",expe+"_"+tp,acol, acol+"_"+cond+"_"+"Diffxpy_Statistics.tsv"), sep="\t")

    # tester = tester.loc[(tester["gene"].isin(highly_variable))]
    # rnk = tester[["gene", "log2fc"]]
    # 
    # pre_res = gseapy.prerank(rnk=rnk,
    #                          gene_sets="data/Bader_GSEA_GMTs/Mouse_Human_Reactome_February_05_2021_symbol.gmt",
    #                          processes=4, permutation_num=1000,
    #                          outdir=os.path.join("results",expe+"_"+tp, acol, 'GSEA_treatment', 'prerank_report_Reactome'), no_plot=True,
    #                          seed=770)
    # filenam = acol + "_prerank_res" + ".csv"
    # pre_res.res2d.to_csv(os.path.join("results",expe+"_"+tp, acol, 'GSEA_treatment', 'prerank_report_Reactome', filenam))
    # 
    # terms = pre_res.res2d.index
    # for i in range(0, len(terms)):
    #     gseaplot(rank_metric=pre_res.ranking, term=terms[i], **pre_res.results[terms[i]],
    #              ofname=os.path.join("results",expe+"_"+tp, acol, "GSEA_treatment", "prerank_report_Reactome",
    #                                  terms[i].split("%")[0].replace("/", "-") + ".png"))
    #     plt.close()

