import pandas as pd
import collections
import statsmodels.stats.multitest as multi
import random
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
CluEnsp__path="clusters.proteins.v11.0.human.txt"
Ensp19353__path="ensp19353.txt"
EnspInterest__path="ensp34.txt"
CluBackCount__path="STRINGCluster_count__19353ensp.txt"
Cl54__path="54cl__327GeneQ0.01.lt200.txt"
TotalGene = 19353
TotalCand = 34
CluDesc__path="clusters.info.v11.0.CluDesc.txt"
PERM=1000
val = 4
QThr=0.1
OUTPNG = "string__Q" + str(QThr) + ".34gene.png"
PERMDATA = "NOfQ" + str(QThr) + "_In34gene" + str(PERM) + "perms_In327geneQ0.01.txt"

with open(Cl54__path,"r") as f:
  Cl54__l = [ l.strip("\n") for l in f.readlines() ]

with open(Ensp19353__path,"r") as f:
  Ensp19353__l = [ l.strip("\n") for l in f.readlines() ]


CluDesc = pd.read_csv(CluDesc__path, sep="\t",names=["cl","desc"]).set_index("cl")
CluHitBack = pd.read_csv(CluBackCount__path, sep="\t",names=["cl","CountBack"]).set_index("cl")
CluEnsp = pd.read_csv(CluEnsp__path, sep="\t",names=["mushi","cl","ensp"]).set_index("cl")

NQ = []
for i in range(PERM):
  # for random candidates
  EnspInterest__l = random.sample(Ensp19353__l,TotalCand)
  CluHitCount = pd.DataFrame.from_dict( collections.Counter(CluEnsp[CluEnsp.ensp.isin(EnspInterest__l)].index.to_list()) , orient="index").rename(columns={0:"CountHit"})
  Cl_CountBack_CountHit = pd.merge(CluHitBack,CluHitCount,left_index=True,right_index=True,how="left")
  Cl_CountBack_CountHit = Cl_CountBack_CountHit.fillna(0)
  p__l = []
  for Cl in Cl_CountBack_CountHit.index.to_list():
    CountBack = Cl_CountBack_CountHit.loc[Cl]["CountBack"]
    CountHit  = Cl_CountBack_CountHit.loc[Cl]["CountHit"]
    p__l.append( hypergeom.sf(CountHit-1 ,TotalGene, CountBack, TotalCand) )
  Cl_CountBack_CountHit["p"] = p__l
  Cl_CountBack_CountHit = Cl_CountBack_CountHit[Cl_CountBack_CountHit.CountBack<200]
  Cl_CountBack_CountHit = Cl_CountBack_CountHit[ Cl_CountBack_CountHit.index.isin(Cl54__l) ]
  q__l = multi.multipletests(Cl_CountBack_CountHit["p"].to_list(), alpha=0.01, method='fdr_bh')[1]
  Cl_CountBack_CountHit["q"] = q__l
  Cl_CountBack_CountHit = pd.merge(Cl_CountBack_CountHit, CluDesc,left_index=True,right_index=True,how="inner")
  NQ.append( len( Cl_CountBack_CountHit[Cl_CountBack_CountHit.q < QThr] ) )

write_list(NQ, PERMDATA)

with open(PERMDATA,"r") as f:
  data__l = [ int(l.strip("\n")) for l in f.readlines() ]

Max1 = max(data__l)
Max2 = max( [Max1, val] )
p = len( [ i for i in data__l if i >= val ] ) / PERM
fig = plt.figure(figsize=(6*0.394, 6*0.394),dpi=900)
plt.xlim(0,Max2+1)
plt.hist( data__l ,bins=Max1)
plt.axvline(x=val,linestyle="dashed",color="red")
plt.title("p-value: " + str(p) )
fig.savefig(OUTPNG)
plt.close()
