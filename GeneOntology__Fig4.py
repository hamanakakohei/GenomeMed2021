import pandas as pd
import json
from pandas import json_normalize
import subprocess
import random
import codecs
import statsmodels.stats.multitest as multi
import numpys as np
import matplotlib.pyplot as plt

N=34
JsonFile="tmp.toppfun.random"
Q=0.01
GoFile__328gene="328geneThr20-500fdr1.0.txt"
Entrez19277="enst_entrezid.hg1938.19277.txt"
Entrez52genes="entrez__34gene.txt"

Entrezs__l = pd.read_csv(Entrez19277,sep="\t",header=None)[1].to_list()
Entrez52__l = pd.read_csv(Entrez52genes,sep="\t",header=None)[0].to_list()
with codecs.open(GoFile__328gene,"r","utf-8","ignore") as f:
  Go328__p = pd.read_csv(f,delimiter="\t")

Go328q__p = Go328__p[Go328__p["q-value FDR B&H"] < 0.01]
Mf328q__l = Go328q__p[Go328q__p.Category=="GO: Molecular Function"]["ID"].to_list()
Bp328q__l = Go328q__p[Go328q__p.Category=="GO: Biological Process"]["ID"].to_list()
Cc328q__l = Go328q__p[Go328q__p.Category=="GO: Cellular Component"]["ID"].to_list()

NMf__l = []
NBp__l = []
NCc__l = []
for i in range(1):
  RandGenes = ",".join([ str(i) for i in random.sample(Entrezs__l,N) ])
  Command = 'curl -H \'Content-Type: text/json\' -d \'{"Genes": [' + RandGenes + '],"Categories": [{"Type": "GeneOntologyMolecularFunction", "PValue": 1.0, "MinGenes": 20, "MaxGenes": 500, "MaxResults": 10000, "Correction": "FDR"},{"Type": "GeneOntologyBiologicalProcess","PValue": 1.0,"MinGenes": 20,"MaxGenes": 500,"MaxResults": 10000,"Correction": "FDR"}, {"Type": "GeneOntologyCellularComponent","PValue": 1.0,"MinGenes": 20,"MaxGenes": 500,"MaxResults": 10000,"Correction": "FDR"}] }\' https://toppgene.cchmc.org/API/enrich > tmp.toppfun.random' + str(Q) + '.txt'
  subprocess.call(Command,shell=True)
  with open(JsonFile + str(Q) + ".txt","r") as f:
    j = json.load(f)
    Rand__p = json_normalize(j["Annotations"])
  MfRand__p = Rand__p[Rand__p.Category=="GeneOntologyMolecularFunction"]
  BpRand__p = Rand__p[Rand__p.Category=="GeneOntologyBiologicalProcess"]
  CcRand__p = Rand__p[Rand__p.Category=="GeneOntologyCellularComponent"]
  MfRandP__l = MfRand__p[MfRand__p.ID.isin(Mf328q__l)]["PValue"].to_list()
  BpRandP__l = BpRand__p[BpRand__p.ID.isin(Bp328q__l)]["PValue"].to_list()
  CcRandP__l = CcRand__p[CcRand__p.ID.isin(Cc328q__l)]["PValue"].to_list()
  MfRandPQ__l = multi.multipletests(MfRandP__l, alpha=0.1, method='fdr_bh')
  BpRandPQ__l = multi.multipletests(BpRandP__l, alpha=0.1, method='fdr_bh')
  CcRandPQ__l = multi.multipletests(CcRandP__l, alpha=0.1, method='fdr_bh')
  N_MfRandPQ = len(MfRandPQ__l[1][MfRandPQ__l[1]<Q])
  N_BpRandPQ = len(BpRandPQ__l[1][BpRandPQ__l[1]<Q])
  N_CcRandPQ = len(CcRandPQ__l[1][CcRandPQ__l[1]<Q])
  NMf__l.append(N_MfRandPQ)
  NBp__l.append(N_BpRandPQ)
  NCc__l.append(N_CcRandPQ)

MfRand__p2 = MfRand__p[MfRand__p.ID.isin(Mf328q__l)]
BpRand__p2 = BpRand__p[BpRand__p.ID.isin(Bp328q__l)]
CcRand__p2 = CcRand__p[CcRand__p.ID.isin(Cc328q__l)]
MfRand__p2["QOfGoInSigOneIn328"] = MfRandPQ__l[1]
BpRand__p2["QOfGoInSigOneIn328"] = BpRandPQ__l[1]
CcRand__p2["QOfGoInSigOneIn328"] = CcRandPQ__l[1]
pd.concat([MfRand__p2,BpRand__p2,CcRand__p2]).to_csv("34geneQ" + str(Q) + "0.1In328GeneQ0.01GOs.txt",sep="\t",index=False)

val = 50
TYPE = "bp" #"cc", "mf"
PERM=1000
OUT = TYPE + "Q" + str(Q) + "__34gene.png"
DATA = "NOfQ" + str(Q) + TYPE + "__34gene" + str(PERM) + "perms.txt"
with open(DATA,"r") as f:
  data__l = [ int(l.strip("\n")) for l in f.readlines() ]

Max1 = max(data__l)
Max2 = max( [Max1, val] )
p = len( [ i for i in data__l if i >= val ] ) / PERM
fig = plt.figure(figsize=(6*0.394, 6*0.394),dpi=900)
plt.xlim(0,Max2+1)
plt.hist( data__l ,bins=Max1)
plt.axvline(x=val,linestyle="dashed",color="red")
plt.title("p-value: " + str(p) )
fig.savefig(OUT)
plt.close()

