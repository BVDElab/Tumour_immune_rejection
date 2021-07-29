import numpy as np
import pandas as pd
import os
import random
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.exceptions import ConvergenceError
from collections import Counter
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['svg.fonttype'] = 'none'


random_seed=770

# Reproducible research 101
random.seed(random_seed)
np.random.seed(random_seed)

# Datasets from Xena
metadata=pd.read_csv(os.path.join("Data", "metadata_tcga_to_xena.txt"), sep="\t", encoding="latin-1")

indata=pd.read_csv("Data/TCGA_survival_data", sep="\t")
indata=indata[["sample", "OS", "OS.time"]].copy()
indata=pd.merge(left=indata, right=metadata, left_on="sample", right_on="sample", how="inner")

# Gene expression Thorsson et al 2018 / GDC
rnaseq=pd.read_csv("Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", sep="\t")
rnaseq.set_index("gene_id", inplace=True)
rnaseq=rnaseq.transpose()
thecols=[x for x in list(rnaseq) if x.find("ADRA")!=-1]
thecols=thecols+["HBA1|3039", "HBA2|3040", "HBB|3043"]

rnaseq.reset_index(drop=False, inplace=True)
rnaseq=rnaseq[thecols+["index"]].copy()
rnaseq["index"]=["-".join(x.split("-")[:-3])[:-1] for x in list(rnaseq["index"])]
indata=pd.merge(left=indata, right=rnaseq, left_on="sample", right_on="index", how="inner")
indata.dropna(subset=["sample"], inplace=True)

indata["sample_type"]=indata["sample_type"].replace({"Additional Metastatic": "Metastatic"})
indata=indata.loc[indata["sample_type"].isin(["Primary Tumor", "Metastatic"])].copy()

counts=Counter(indata["cancer"])
print(counts)

lister=[]
ctypes=counts.keys()

indata=indata[list(indata)[:-1]].copy()

thecols=thecols
indata.dropna(how="any", inplace=True)

for agene in thecols:
    subs=indata.copy()
    subs=subs.loc[~pd.isna(subs[agene])].copy()
    subs[agene]=subs[agene].apply(np.log2)

    for acancer in ctypes:
        subf=indata.copy()
        subf=subf.loc[subf["cancer"]==acancer].copy()
        subf.reset_index(drop=True, inplace=True)

        subf=subf.loc[~pd.isna(subf["OS"])].copy()
        subf=subf.loc[~pd.isna(subf["OS.time"])].copy()


        subf=subf.loc[~pd.isna(subf[agene])].copy()

        if len(subf)>30:
            mymedian=np.median(subf[agene])
            subf["ExpressionState"]=["High" if x>mymedian else "Low" for x in list(subf[agene])]

            low=subf.loc[subf["ExpressionState"]=="Low"]
            high=subf.loc[subf["ExpressionState"]=="High"]

            kmf_low = KaplanMeierFitter()
            kmf_low.fit(durations=low["OS.time"], event_observed=low["OS"], label=r"$\it{"+str(agene).split("|")[0]+"}$ (Low)")
            kmf_high = KaplanMeierFitter()
            kmf_high.fit(durations=high["OS.time"], event_observed=high["OS"], label=r"$\it{"+str(agene).split("|")[0]+"}$ (High)")

            timelow = low["OS.time"]
            timehigh = high["OS.time"]

            eventlow = low["OS"]
            eventhigh = high["OS"]

            logres = logrank_test(timelow, timehigh, eventlow, eventhigh)

            #if logres.p_value<0.05:
            f = plt.figure(figsize=(2.91, 4.10))
            ax=f.add_subplot(111)

            # kmf.plot()
            kmf_low.plot(ax=ax, ci_show=False)
            kmf_high.plot(ax=ax, ci_show=False)

            plt.title(r"$\bf{"+acancer + "}$"+ "\n" + "LogRank pval= " + str(np.round(logres.p_value, 6))+"\nN_patients="+str(len(subf)))
            plt.xlim(6000)
            plt.tight_layout()
            sns.despine()
            plt.savefig("KMcurve_"+agene.split("|")[0]+"_"+acancer+".svg")
            plt.close()