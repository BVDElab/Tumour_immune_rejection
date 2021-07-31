import numpy as np
import pandas as pd
import os
import random
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from collections import Counter

random_seed=770

# Reproducible research 101
random.seed(random_seed)
np.random.seed(random_seed)


# ######################## #
# Processing XCELL data #
# ######################## #


clinical=pd.read_csv("Data/clinical_PANCAN_patient_with_followup.tsv", sep="\t", encoding="cp1252")
metadata=pd.read_csv("Data/Metadata-TCGA-All-18116-Samples.csv")
metadata["case_uuid"]=[x.upper() for x in metadata["case_uuid"]]
metadata["gdc_file_uuid"]=[x.upper() for x in metadata["gdc_file_uuid"]]
metadata=metadata[["gdc_file_uuid", "aliquot_uuid", "primary_site", "sample_type", "investigation", "filename"]].copy()


mdic={}
mdat=open("Data/TCGAutils_tcgamap.txt", "r")
for aline in mdat.read().split("\n"):
    try:
        mdic[aline.split("\t")[1].upper().replace('"','')]=aline.split("\t")[2].replace('"','')
    except IndexError:
        print(aline)
mdat.close()

metadata["aliquot"]=metadata["aliquot_uuid"].map(mdic)
metadata["sample"]=["-".join(x.split("-")[:-3])[:-1] for x in list(metadata["aliquot"])]



# selection frame
infils=pd.read_csv("Data/xCell_TCGA_RSEM.txt", sep="\t")
infils.set_index("Unnamed: 0", inplace=True)
infils=infils.transpose()
infils.reset_index(drop=False, inplace=True)
infils["index"]=[x.replace(".","-") for x in list(infils["index"])]
indata=pd.merge(left=metadata, right=infils, left_on="sample", right_on="index", how="left")
#indata.drop_duplicates(subset=["index"], keep="first", inplace=True)

surv=pd.read_csv("Data/TCGA_survival_data", sep="\t")
surv=surv[["sample", "OS", "OS.time", "PFI", "PFI.time"]]
indata=pd.merge(left=indata, right=surv, left_on="sample", right_on="sample", how="left")

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
indata["investigation"]=[x.replace("TCGA-","") for x in list(indata["investigation"])]

indata["sample_type"]=indata["sample_type"].replace({"Additional Metastatic": "Metastatic"})
indata=indata.loc[indata["sample_type"].isin(["Primary Tumor", "Metastatic"])].copy()

from collections import Counter
counts=Counter(indata["investigation"])
from lifelines.exceptions import ConvergenceError


lister=[]
ctypes=sorted(set(indata["investigation"]))

mergedcols=["ADRA2A|150", "ADRA2B|151", "ADRA2C|152"]
indata["ADRA_SUM"]= indata["ADRA2A|150"]+indata["ADRA2B|151"]+indata["ADRA2C|152"]
indata["ADRA_mean"]= indata["ADRA2A|150"]+indata["ADRA2B|151"]+indata["ADRA2C|152"]
indata["ADRA_mean"]=indata["ADRA_mean"]/3
indata["ADRA_MAX"]=[np.max([x[0], x[1], x[2]]) for x in zip(list(indata["ADRA2A|150"]), list(indata["ADRA2B|151"]), list(indata["ADRA2C|152"]))]


for agene in thecols:
    subs=indata.copy()
    subs=subs.loc[~pd.isna(subs[agene])].copy()
    subs[agene]=subs[agene].apply(np.log2)
    sns.boxplot(x="investigation", y=agene, data=subs, hue="sample_type", order=ctypes,flierprops = dict(markerfacecolor = '0.50', markersize = 2))
    plt.xticks(rotation=90)
    plt.ylabel(r"Log$_2$ gene expression ("+agene+")")
    sns.despine()
    plt.tight_layout()
    plt.savefig("Boxplot_"+agene+".pdf")
    plt.close()

    for acancer in ctypes:
        subf=indata.copy()
        subf=subf.loc[subf["investigation"]==acancer].copy()
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
            kmf_low.fit(durations=low["OS.time"], event_observed=low["OS"], label="Low")
            kmf_high = KaplanMeierFitter()
            kmf_high.fit(durations=high["OS.time"], event_observed=high["OS"], label="High")

            timelow = low["OS.time"]
            timehigh = high["OS.time"]

            eventlow = low["OS"]
            eventhigh = high["OS"]

            logres = logrank_test(timelow, timehigh, eventlow, eventhigh)

            f = plt.figure(figsize=(3, 3))
            ax=f.add_subplot(111)

            # kmf.plot()
            kmf_low.plot(ax=ax, ci_show=False)
            kmf_high.plot(ax=ax, ci_show=False)

            plt.title(agene+" in " + acancer + "\n" + "LogRank pval= " + str(np.round(logres.p_value, 6))+"\nN_patients="+str(len(subf)))
            plt.tight_layout()
            sns.despine()
            plt.savefig("KMcurve_"+agene.split("|")[0]+"_"+acancer+".pdf")
            plt.close()
