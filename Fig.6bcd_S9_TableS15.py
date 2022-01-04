import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import feature_column
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn import metrics
from sklearn.neighbors import KernelDensity
import collections
import statistics
import os
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import norm
import scipy.stats as st
import decimal
import seaborn
import math
import copy
import itertools
from google.colab import files

def auc(X, Y):
    return 1/(len(X)*len(Y)) * sum([kernel(x, y) for x in X for y in Y])

def kernel(X, Y):
    return .5 if Y==X else int(Y < X)

def structural_components(X, Y):
    V10 = [1/len(Y) * sum([kernel(x, y) for y in Y]) for x in X]
    V01 = [1/len(X) * sum([kernel(x, y) for x in X]) for y in Y]
    return V10, V01
    
def get_S_entry(V_A, V_B, auc_A, auc_B):
    return 1/(len(V_A)-1) * sum([(a-auc_A)*(b-auc_B) for a,b in zip(V_A, V_B)])

def z_score(var_A, var_B, covar_AB, auc_A, auc_B):
    return (auc_A - auc_B)/((var_A + var_B - 2*covar_AB)**(.5))

def plot_history(histories, key='binary_crossentropy'):
  plt.figure(figsize=(16,10))
  for name, history in histories:
    val = plt.plot(history.epoch, history.history['val_'+key],'--', label=name.title()+' Val')
    plt.plot(history.epoch, history.history[key], color=val[0].get_color(),label=name.title()+' Train')
  plt.xlabel('Epochs')
  plt.ylabel(key.replace('_',' ').title())
  plt.legend()
  plt.xlim([0,max(history.epoch)])

def set_seed(seed=200):
  tf.random.set_seed(seed)
  # for hash seed
  os.environ["PYTHONHASHSEED"] = str(seed)

def df_to_dataset(dataframe, shuffle=True, batch_size=32):
  dataframe = dataframe.copy()
  labels = dataframe.pop('target')
  #labels = pd.concat([dataframe.pop(x) for x in ['target1', 'target2']], 1)
  ds = tf.data.Dataset.from_tensor_slices((dict(dataframe), labels))
  if shuffle:
    ds = ds.shuffle(buffer_size=len(dataframe))
  ds = ds.batch(batch_size)
  return ds

def prep_data(DATA):
  df0 = pd.read_csv(DATA,sep="\t")  #22933*19
  df0["NaCand"]  = [ 1 if i else 0 for i in pd.isna(df0["cand"]).to_list() ]
  df0["NaTrain"] = [ 1 if i else 0 for i in pd.isna(df0["train"]).to_list() ]
  df0["NaCount"] = df0["NaCand"] + df0["NaTrain"]
  ddf = df0.groupby("enst")
  df0 = df0.loc[ddf["NaCount"].idxmin(),:]  #20034*22
  pc = df0.query("my=='yes' & train=='yes'") #328*22
  nc = df0.query("my=='yes' & cand!='yes' & train!='yes'" ); nc = nc[nc["inheritance"].isnull()] #18012*22
 
  df = pd.concat([pc.assign(target=1),nc.assign(target=0)])
  df = df.drop(["ensp","gene","enst380","cand","train","inheritance","my","hc","NaCand","NaTrain","NaCount"],axis=1) #enst #18340*12
  ot__index = df0.index.difference(pc.index).difference(nc.index) #1694
  ot = df0[[ True if i in ot__index.to_list() else False for i in df0.index.to_list() ]] #1694*22
  ot__pc = ot[ot.inheritance.notnull()] #1646*22
  ot__nc = ot[ot.inheritance.isnull()] #48*22
  ot = pd.concat([ot__pc.assign(target=1),ot__nc.assign(target=0)])
  ot = ot.drop(["ensp","enst380","cand","train","inheritance","my","hc","NaCand","NaTrain","NaCount"],axis=1) #enst #1694*13
 
  df = df.dropna(how="any") #12565*12
  df['tsea'] = pd.Categorical(df['tsea'])
  df['tsea'] = df.tsea.cat.codes
  df['string'] = pd.Categorical(df['string'])
  df['string'] = df.string.cat.codes
  df['target'] = pd.Categorical(df['target'])
  df['target'] = df.target.cat.codes
  df = df.replace({"module":{"M2":"other","M3":"other","M5":"other","M6":"other","M8":"other","M9":"other","M10":"other","M11":"other","M12":"other","M14":"other","M15":"other","M16":"other","M17":"other","M18":"other"}})
 
  # fill NA for last 34 HC new gene analysis
  ot_34 = ot[ot.entrez.isin([84078,5789,26528,84444,23396,11198,9646,26038,3736,1995,57555,57599,5978,7514,65986,5702,23031,9703,8314,8019,2011,81603,23468,2932,84146,130507,283450,4440,4691,575,64864,29998,3066,4735])]
  ot_34 = ot_34.fillna({"pLI":np.nanmedian(ot_34.pLI), "oe_lof_upper":np.nanmedian(ot_34.oe_lof_upper), "mis_z":np.nanmedian(ot_34.mis_z),"tsea":np.nanmean(ot_34.tsea), "subre":np.nanmedian(ot_34.subre), "module":"other","go":np.nanmedian(ot_34.go)})
  ot_34['tsea']   = pd.Categorical(ot_34['tsea'])
  ot_34['tsea']   = ot_34.tsea.cat.codes
  ot_34['target'] = pd.Categorical(ot_34['target'])
  ot_34['target'] = ot_34.target.cat.codes
  ot_34['string'] = pd.Categorical(ot_34['string'])
  ot_34['string'] = ot_34.string.cat.codes
  ot_34 = ot_34.replace({"module":{"M2":"other","M3":"other","M5":"other","M6":"other","M8":"other","M9":"other","M10":"other","M11":"other","M12":"other","M14":"other","M15":"other","M16":"other","M17":"other","M18":"other"}})

  ot = ot.dropna(how="any")
  ot['tsea']   = pd.Categorical(ot['tsea'])
  ot['tsea']   = ot.tsea.cat.codes
  ot['target'] = pd.Categorical(ot['target'])
  ot['target'] = ot.target.cat.codes
  ot['string'] = pd.Categorical(ot['string'])
  ot['string'] = ot.string.cat.codes
  ot = ot.replace({"module":{"M2":"other","M3":"other","M5":"other","M6":"other","M8":"other","M9":"other","M10":"other","M11":"other","M12":"other","M14":"other","M15":"other","M16":"other","M17":"other","M18":"other"}})

  # all genes
  all = df0.drop(["ensp","enst380","cand","train","inheritance","my","hc","NaCand","NaTrain","NaCount"],axis=1)
  all = all.fillna({"pLI":np.nanmedian(all.pLI), "oe_lof_upper":np.nanmedian(all.oe_lof_upper), "mis_z":np.nanmedian(all.mis_z),"tsea":np.nanmean(all.tsea), "subre":np.nanmedian(all.subre), "module":"other","go":np.nanmedian(all.go)})
  all['target'] = [1 for i in range(10000)] + [0 for i in range(10034)]
  all['tsea']   = pd.Categorical(all['tsea'])
  all['tsea']   = all.tsea.cat.codes
  all['target'] = pd.Categorical(all['target'])
  all['target'] = all.target.cat.codes
  all['string'] = pd.Categorical(all['string'])
  all['string'] = all.string.cat.codes
  all = all.replace({"module":{"M2":"other","M3":"other","M5":"other","M6":"other","M8":"other","M9":"other","M10":"other","M11":"other","M12":"other","M14":"other","M15":"other","M16":"other","M17":"other","M18":"other"}})

  return df0,df,ot,ot_34,all


def NN(df,ot,ot_34,all):
  trashForNc3, nc2  = train_test_split(df.query("target==0"), test_size=0.08138, random_state=0)
  NoNc2__index = df.index.difference(nc2.index)
  df2 = df[[ True if i in NoNc2__index.to_list() else False for i in df.index.to_list() ]]
  train, test = train_test_split(df2, test_size=0.1, random_state=0) #0.00000001) 

  feature_columns = []
  for header in ["pLI","oe_lof_upper","mis_z","tsea","subre","go","string"]:
    feature_columns.append(feature_column.numeric_column(header))
  module = feature_column.categorical_column_with_vocabulary_list('module', ['M1','M4','M7','M13','other'])
  module_one_hot = feature_column.indicator_column(module)
  feature_columns.append(module_one_hot)
  feature_layer = tf.keras.layers.DenseFeatures(feature_columns)
  batch_size = 5
  train_ds = df_to_dataset(train,shuffle=False, batch_size=batch_size)
  test_ds  = df_to_dataset(test, shuffle=False, batch_size=batch_size) ##########################
  ot_ds    = df_to_dataset(ot, shuffle=False, batch_size=batch_size)
  ot34_ds  = df_to_dataset(ot_34, shuffle=False, batch_size=batch_size)
  nc2_ds   = df_to_dataset(nc2,shuffle=False, batch_size=batch_size)
  all_ds   = df_to_dataset(all,shuffle=False, batch_size=batch_size)

  model = tf.keras.Sequential([ feature_layer, layers.Dense(128, activation='relu'), layers.Dense(128, activation='relu'), layers.Dense(1, activation='sigmoid') ])
  model.compile(optimizer='adam',loss='binary_crossentropy',metrics=['accuracy'])
  model_hist = model.fit(train_ds, epochs=5)

  train["nn"] = [ i[0] for i in model.predict(train_ds) ]
  test["nn"]  = [ i[0] for i in model.predict(test_ds) ] #################################
  ot["nn"]    = [ i[0] for i in model.predict(ot_ds) ]
  ot_34["nn"] = [ i[0] for i in model.predict(ot34_ds) ]
  nc2["nn"]   = [ i[0] for i in model.predict(nc2_ds) ]
  all["nn"]   = [ i[0] for i in model.predict(all_ds) ]

  ot_ddg = ot[ot.entrez.isin([473,23065,10236,79947,9204,26523,4520,23499,23499,11146,1301,476,8364,1141,103,777,57165,10000,56896,6654,25885,3899,5903,2736,785,90,64135,6328,4036,2335,10058,5077,9901,5894,5915,93,11280,776,285195,7545,11342,2122,1857,1173,2257,3516,2895,5530,10586,25836,348980,57688,3796,6228,7025,2201,5307,51780,5159,23291,7280,347733,221692,6659,8364,8364,8364,8364,8364,8364,8364,8364,8364,8364,203068,860,56479,9096,8936,116150,2697,7432,60,816,7532,9863,54809,2783,6608,2146,23221,1135,2138,23414,5885,2131,23513,23189,23586,80036,4915,9568,6709,7248,57582,22884,11155,9211,5832,3265,10522,9024,3746,203859,5080,7490,4076,6506,2132,50801,9379,6712,4041,8726,10413,2893,867,54538,84623,894,1822,4040,7846,6602,57609,6334,6601,23592,59341,23316,26960,1282,5587,27133,7253,8812,23241,161725,6263,5888,53944,53944,8766,8912,23162,112476,112755,5336,1013,10381,10381,5048,5430,6844,4628,4626,4621,93649,9611,84282,6928,7703,7067,7283,2670,5718,2186,60,5881,23347,63895,10939,1000,23327,1938,5605,10382,9230,22983,112939,4854,6324,7040,1760,2906,50944,8200,2036,57167,1137,4661,150094,8216,23384,9681,4627,57502,412,4281,9758,5277,8905,4810,6611,80311,1968,5422,170302,11141,2710,1756,7102,10159,4128,5127,6853,24140,10682,10084,11152,28952,389856,23708,3028,23133,2245,81887,4983,1947,1896,24137,1741,54413,8473,51132,84061,1349,50945,254065,7552,27286,3188,5354,2182,5063,79868,186,7319,65109,2892,4952,51114,3547,2239,2719,10479,2273,9459,6658,2334,3423,10046,6535,10134,215,3897,3054,2316,2664,116442,1193])]

  X_plot = np.arange(0,1.0,0.01).reshape(-1,1)
  X = np.array(nc2.nn).reshape(-1,1)
  kdeNc2 = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(X) #0.015
  X = np.array(ot_ddg.nn).reshape(-1,1)
  kdeDdg = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(X)  #0.02

  ot_34["prior"] = (52-380*0.05)/52
  ot_34["preodds"] = ot_34["prior"]/(1-ot_34["prior"])
  ot_34["nnlr"] = np.exp(kdeDdg.score_samples(ot_34.nn.values.reshape(-1,1))) / np.exp(kdeNc2.score_samples(ot_34.nn.values.reshape(-1,1)))
  ot_34["posodds"] = ot_34["preodds"]*ot_34["nnlr"]
  ot_34["poste"] = ot_34["posodds"]/(1+ot_34["posodds"])
  ot_34["nnrank"] = ot_34.nn.rank()
  ot_34["nn2"] = [ decimal.Decimal(i)+0 for i in ot_34.nn ]
  ot_34["xlab"] = ot_34.gene + " (" + ot_34.nn2.astype(str) + ")"
  ot_34 = ot_34.fillna({"poste":1})

  return train,test,ot,ot_34,nc2,all,ot_ddg



##set_seed(0)
DATA="tableS15.txt"
aucs_cds = []
for i in range(1):
  df0,df,ot,ot_34,all = prep_data(DATA)
  train,test,ot,ot_34,nc2,all,ot_ddg = NN(df,ot,ot_34,all)
  pc3 = ot_ddg.copy()
  pc3["target"] = 1
  nc3 = pd.merge(df0[["enst","gene"]],nc2)
  nc3["target"] = 0
  pc3_nc3 = pd.concat([pc3,nc3])
  auc = roc_auc_score(pc3_nc3['target'],pc3_nc3['nn'])
  aucs_cds = aucs_cds + [auc]


plt.rcParams["font.size"] = 8
fig = plt.figure(figsize=(4.5*0.394, 4.5*0.394),dpi=900)
ax = fig.add_subplot(1, 1, 1)
ax.violinplot([
  train.query("target==0").nn.to_list(),
  train.query("target==1").nn.to_list(),
  test.query("target==0").nn.to_list(),
  test.query("target==1").nn.to_list()],widths=1, showmeans=False, showmedians=True, showextrema=False)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['NC1', 'PC1', 'NC2', 'PC2'])
fig.subplots_adjust(left=0.2)
plt.show()
plt.close()
fig.savefig("SubsetModel.png")


fig = plt.figure()
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax1.hist(nc2.nn.to_list(),   density=True, bins=15) #
ax2.hist(ot_ddg.nn.to_list(),density=True, bins=15)
ax1.plot(X_plot[:, 0], np.exp(kdeNc2.score_samples(X_plot)),color="orange")
ax2.plot(X_plot[:, 0], np.exp(kdeDdg.score_samples(X_plot)),color="orange")


gene_nnmean = gene_nns.median(axis=1)
gene_nnmean.name = "nn"
gene_postemean = gene_postes.median(axis=1)
gene_postemean.name = "poste"
gene_nn_poste = pd.merge(gene_nnmean,gene_postemean,left_index=True,right_index=True).reset_index()
gene_nn_poste["xlab"] = gene_nn_poste.gene + " (" + gene_nn_poste.nn.round(3).astype(str) + ")"
gene_nn_poste["nnrank"] = gene_nn_poste.nn.rank()
gene_nn_poste.to_csv("tableS12median.txt",sep="\t")

decimal.getcontext().prec = 2
plt.rcParams["font.size"] = 8

WIDTH = 15
HEIGHT = 5
fig = plt.figure(figsize=(WIDTH*0.394, HEIGHT*0.394),dpi=900)
ax = fig.add_subplot(1, 1, 1)
plt.scatter(gene_nn_poste.nnrank, gene_nn_poste.poste)
plt.xticks(gene_nn_poste.nnrank.to_list(), gene_nn_poste.xlab.to_list(),rotation=270) #,fontsize = 7
plt.hlines(0.9, 0, 35, colors='black', linestyle='dashed', linewidth=1)
plt.xlim(0,35)
plt.ylim(-0.05,1.05)
plt.subplots_adjust(bottom=0.5)
plt.show()
plt.close()
fig.savefig("PP" + str(WIDTH) + "-" + str(HEIGHT) + ".median.png")


# vs RVIS, GDI, HGC

pc12__276gene = pd.merge(df0[["enst","gene"]],pd.concat([train,test]).query('target==1')).rename({"gene":"Gene"},axis=1)["Gene"].to_list()
GDI_path = 'GDI_full.txt'
RVIS_path = 'RVIS_Unpublished_ExACv2_March2017.txt' 
HGCS_path = 'pc3_nc3__hgcs_additional26.txt'
Gene_GDI = pd.read_table(GDI_path)[['Gene','GDI']]
Gene_RVIS = pd.read_table(RVIS_path).rename({"CCDSr20":"Gene", "RVIS[pop_maf_0.05%(any)]":"RVIS"},axis=1)[["Gene","RVIS"]] #RVIS v4, constructed on ExAC v2
Gene_HGC = pd.read_table(HGCS_path)[["Target","Distance","Rank"]].dropna(how="any").rename({"Target":"Gene"},axis=1).drop_duplicates()
Gene_Rank = Gene_HGC.loc[Gene_HGC.groupby('Gene')['Rank'].idxmin(),:][["Gene","Rank"]]
Gene_Distance = Gene_HGC.loc[Gene_HGC.groupby('Gene')['Distance'].idxmin(),:][["Gene","Distance"]]
Gene_predictors = df0[['gene', 'tsea', 'pLI', 'oe_lof_upper', 'mis_z', 'subre','module', 'go', 'string']].rename({"gene":"Gene"},axis=1).drop_duplicates(subset="Gene")
Gene_predictors = Gene_predictors.replace({"module":{"M1":1,"M2":0,"M3":0,"M4":1,"M5":0,"M6":0,"M7":1,"M8":0,"M9":0,"M10":0,"M11":0,"M12":0,"M13":1,"M14":0,"M15":0,"M16":0,"M17":0,"M18":0,"other":0}})
Gene_predictors['string'] = pd.Categorical(Gene_predictors['string'])
Gene_predictors['string'] = Gene_predictors.string.cat.codes

pc3 = ot_ddg
pc3["target"] = 1
nc3 = pd.merge(df0[["enst","gene"]],nc2)
nc3["target"] = 0
Gene_target_nn_GDI_RVIS_Rank_Distance =pd.merge(pd.merge(pd.merge(pd.merge(pd.concat([pc3,nc3]).rename({"gene":"Gene"},axis=1)[["Gene","target","nn"]], Gene_GDI,how="left"),Gene_RVIS,how="left"),Gene_Rank,how="left"),Gene_Distance,how="left")
pc3nc3__1247gene = Gene_target_nn_GDI_RVIS_Rank_Distance["Gene"].to_list()
AllComb = list(itertools.product(pc12__276gene,pc3nc3__1247gene))
Core_Gene_Target = pd.DataFrame({"Core_Gene":[i[0] for i in AllComb],"Target":[i[1] for i in AllComb]})
Gene_Distance_Rank__median = pd.merge(Core_Gene_Target,(pd.read_table(HGCS_path)[["Core_Gene","Target","Distance","Rank"]]),how="left").fillna({"Distance":92.23896, "Rank":16723}).rename({"Target":"Gene"},axis=1)[["Gene","Distance","Rank"]].groupby('Gene').median().rename({"Distance":"Distance_median","Rank":"Rank_median"},axis=1).reset_index()
Gene_target_nn_GDI_RVIS_Rank_Distance =pd.merge(Gene_target_nn_GDI_RVIS_Rank_Distance, Gene_Distance_Rank__median,how="left")
Gene_target_nn_GDI_RVIS_Rank_Distance_predictors = pd.merge(Gene_target_nn_GDI_RVIS_Rank_Distance,Gene_predictors,how="left")

VAR = "pLI"         ; auc_pli = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "oe_lof_upper"; auc_loe = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "mis_z"       ; auc_mis = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "tsea"        ; auc_tse = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "subre"       ; auc_sub = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "module"      ; auc_mod = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "go"          ; auc_goo = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "string"      ; auc_str = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "RVIS"        ; auc_rvi = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "GDI"         ; auc_gdi = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])
VAR = "HGC"         ; auc_dis = roc_auc_score(Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)["target"], Gene_target_nn_GDI_RVIS_Rank_Distance_predictors.dropna(subset=[VAR],axis=0)[VAR])

WIDTH = 5
HEIGHT = 5
fig = plt.figure(figsize=(WIDTH*0.394,HEIGHT*0.394),dpi=900)
plt.violinplot([aucs_cds],widths=1,showmeans=False,showextrema=False,showmedians=False,vert=True)
plt.scatter([1,2,3,4,5,6,7,8,9,10,11,12],[np.median(aucs_cds),auc_pli,1-auc_loe,auc_mis,auc_tse,auc_sub,auc_mod,auc_goo,auc_str,1-auc_rvi,1-auc_gdi,1-auc_ran],c='red',s=3,alpha=1)
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12],["NN","pLI","LOEUF","mis-z","TSEA","Subregion","Module","GO","STRING","RVIS","GDI","HGC"])
plt.xticks(rotation=-90)
plt.yticks([0.5,0.6,0.7,0.8,0.9],["0.5","0.6","0.7","0.8","0.9"])
plt.ylabel("AUC for PC3 and NC3")
plt.ylim(0.5,0.9)
plt.subplots_adjust(bottom=0.5)
plt.show()
plt.close()
fig.savefig("auc" + str(WIDTH) + "-" + str(HEIGHT) + ".png")
