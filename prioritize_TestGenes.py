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
import statistics
import os
import argparse

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

#parser = argparse.ArgumentParser(description='aaa')    
#parser.add_argument('arg1', default='TrainingGenes.tsv', help='aaa')
#parser.add_argument('arg2', default='TestGenes.tsv', help='aaa')
#parser.add_argument('arg3', default='TestGenesScores.tsv', help='aaa')
#args = parser.parse_args()   
#TrainingGenes = args.arg1
#TestGenes = args.arg2
#TestGenesScores = args.arg3
TrainingGenes = 'TrainingGenes.tsv'
TestGenes = 'TestGenes.tsv'
TestGenesScores = 'TestGenesScores.tsv'

df = pd.read_csv(TrainingGenes,sep="\t")
ot_34 = pd.read_csv(TestGenes,sep="\t") 

df = df.dropna(how="any")
df['tsea'] = pd.Categorical(df['tsea'])
df['tsea'] = df.tsea.cat.codes
df['string'] = pd.Categorical(df['string'])
df['string'] = df.string.cat.codes
df['target'] = pd.Categorical(df['target'])
df['target'] = df.target.cat.codes
df = df.replace({"module":{"M2":"other","M3":"other","M5":"other","M6":"other","M8":"other","M9":"other","M10":"other","M11":"other","M12":"other","M14":"other","M15":"other","M16":"other","M17":"other","M18":"other"}})
 
# fill NA for last 34 HC new gene analysis
ot_34 = ot_34.fillna({"pLI":np.nanmedian(ot_34.pLI), "oe_lof_upper":np.nanmedian(ot_34.oe_lof_upper), "mis_z":np.nanmedian(ot_34.mis_z),"tsea":np.nanmean(ot_34.tsea), "subre":np.nanmedian(ot_34.subre), "module":"other","go":np.nanmedian(ot_34.go)})
ot_34['tsea']   = pd.Categorical(ot_34['tsea'])
ot_34['tsea']   = ot_34.tsea.cat.codes
ot_34['target'] = pd.Categorical(ot_34['target'])
ot_34['target'] = ot_34.target.cat.codes
ot_34['string'] = pd.Categorical(ot_34['string'])
ot_34['string'] = ot_34.string.cat.codes
ot_34 = ot_34.replace({"module":{"M2":"other","M3":"other","M5":"other","M6":"other","M8":"other","M9":"other","M10":"other","M11":"other","M12":"other","M14":"other","M15":"other","M16":"other","M17":"other","M18":"other"}})

trash, nc3  = train_test_split(df.query("target==0"), test_size=0.08138, random_state=0)
NoNc3__index = df.index.difference(nc3.index)
df2 = df[[ True if i in NoNc3__index.to_list() else False for i in df.index.to_list() ]]
train, test = train_test_split(df2, test_size=0.00000001, random_state=0) #0.00000001

feature_columns = []
for header in ["pLI","oe_lof_upper","mis_z","tsea","subre","go","string"]:
  feature_columns.append(feature_column.numeric_column(header))
module = feature_column.categorical_column_with_vocabulary_list('module', ['M1','M4','M7','M13','other'])
module_one_hot = feature_column.indicator_column(module)
feature_columns.append(module_one_hot)
feature_layer = tf.keras.layers.DenseFeatures(feature_columns)
batch_size = 5
train_ds = df_to_dataset(train,shuffle=False, batch_size=batch_size)
ot34_ds  = df_to_dataset(ot_34, shuffle=False, batch_size=batch_size)

model = tf.keras.Sequential([ feature_layer, layers.Dense(128, activation='relu'), layers.Dense(128, activation='relu'), layers.Dense(1, activation='sigmoid') ])
model.compile(optimizer='adam',loss='binary_crossentropy',metrics=['accuracy'])
model.fit(train_ds, epochs=5)

ot_34["nn"] = [ i[0] for i in model.predict(ot34_ds) ]
ot_34.to_csv(TestGenesScores)
