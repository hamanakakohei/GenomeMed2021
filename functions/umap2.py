import umap
from scipy.sparse.csgraph import connected_components
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import sys

args = sys.argv
sample_pcs = pd.read_csv(args[1],sep="\t")
sample_capture = pd.read_csv(args[2],sep="\t")
n_pc = args[3]
outputtable = args[4]
outputplot = args[5]

sample_pcs = sample_pcs.reset_index().rename(columns={"index":"sample"})
sample_pcs_cap = pd.merge(sample_pcs.iloc[:,0:(n_pc+1)],sample_capture,left_on="sample",right_on="sample",how="left").fillna({"capture": "Unkown"})
embedding = umap.UMAP(random_state=0).fit_transform(sample_pcs_cap.iloc[:,1:(n_pc+1)].values)

cap_sample_umap1_umap2 = pd.DataFrame({
    "umap1": embedding[:,0],
    "umap2": embedding[:,1],
    "sample": sample_pcs_cap["sample"],
    "capture": sample_pcs_cap["capture"]
})
cap_sample_umap1_umap2.to_csv(outputtable,sep="\t",index=False)

fig = plt.figure(figsize=(40,40))
ax = fig.add_subplot(1,1,1)
capset = set(cap_sample_umap1_umap2["capture"])
for CAP in capset:
    cap_sample_umap1_umap2_CAP = cap_sample_umap1_umap2[cap_sample_umap1_umap2["capture"] == CAP]
    ax.scatter(data=cap_sample_umap1_umap2_CAP, s=5, alpha=0.2, x='umap1', y='umap2', label=CAP)

#plt.colorbar()
ax.legend(loc="right")
ax.set_xlabel("umap1")
ax.set_ylabel("umap2")
plt.legend(fontsize=20)
plt.show()
plt.savefig(outputplot)
