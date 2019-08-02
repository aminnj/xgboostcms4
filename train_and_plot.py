import os
import sys
import pickle

import numpy as np
import uproot
from tqdm import tqdm
import xgboost as xgb

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score

from utils import write_json, json_to_cfunc

def load_data(fname,feature_names):
    f = uproot.open(fname)
    t = f["t"]
    arrs = t.arrays(t.keys())
    y_data = (arrs["truth"] > 0)
    x_data = np.column_stack([arrs[name] for name in feature_names])
    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.25, random_state=42)
    return x_train,x_test,y_train,y_test

feature_names = [
        "miniiso",
        "ptratio",
        "ptrel",
        "nmiss",
        "hovere",
        ]
x_train,x_test,y_train,y_test = load_data("output_tt.root",feature_names)

dtrain = xgb.DMatrix( x_train, label=y_train)
dtest = xgb.DMatrix( x_test, label=y_test)
evallist  = [(dtrain,'train'), (dtest,'eval')]

num_round = 200

param = {}
param['objective'] = 'binary:logistic'
param['eta'] = 0.3
param['max_depth'] = 3
param['silent'] = 1
param['nthread'] = 15
param['eval_metric'] = "auc"
param['subsample'] = 0.6
param['alpha'] = 8.0
param['gamma'] = 2.0
param['lambda'] = 1.0
param['min_child_weight'] = 1.0
param['colsample_bytree'] = 1.0

sumw_pos = np.abs(dtrain.get_label()==1).sum()
sumw_neg = np.abs(dtrain.get_label()==0).sum()
param["scale_pos_weight"] = sumw_neg/sumw_pos
pklname = "test.pkl"
if not os.path.exists(pklname):
    os.nice(10)
    bst = xgb.train( param.items(), dtrain, num_round, evallist, early_stopping_rounds=10 )
    pickle.dump(bst,open(pklname,"wb"))
    write_json("model.json", bst, feature_names)
    json_to_cfunc("model.json", fname_out="func.h")
else:
    print "[!] Found pickle file...loading"
    bst = pickle.load(open(pklname,"rb"))

y_pred_train = bst.predict(dtrain)
y_pred_test = bst.predict(dtest)

# PLOTTING

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
preds_bkg_test = y_pred_test[y_test==0]
preds_sig_test = y_pred_test[y_test==1]
preds_bkg_train = y_pred_train[y_train==0]
preds_sig_train = y_pred_train[y_train==1]

density = True
bins = np.linspace(0.0,1,50)
edges = 0.5*(bins[:-1]+bins[1:])

from scipy import stats
pks_bkg = stats.ks_2samp(preds_bkg_train,preds_bkg_test)[1]
pks_sig = stats.ks_2samp(preds_sig_train,preds_sig_test)[1]

counts_train_bkg,_,_ = ax.hist(preds_bkg_train, bins=bins,histtype="stepfilled",alpha=0.45, normed=density, label="bkg, train",color=["b"])
counts_train_sig,_,_ = ax.hist(preds_sig_train, bins=bins,histtype="stepfilled",alpha=0.45, normed=density, label="sig, train",color=["r"])
counts_test_bkg,_,_ = ax.hist(preds_bkg_test, bins=bins,histtype="step",alpha=1.0, normed=density, label="bkg, test (KS prob = {:.2f})".format(pks_bkg),color=["b"], lw=1.5, linestyle="solid")
counts_test_sig,_,_= ax.hist(preds_sig_test, bins=bins,histtype="step",alpha=1.0, normed=density, label="sig, test (KS prob = {:.2f})".format(pks_sig),color=["r"], lw=1.5, linestyle="solid")

ax.set_yscale("log")
ax.legend()

ax.set_ylim([0.01,ax.get_ylim()[1]])
fig.set_tight_layout(True)
fig.savefig("disc.pdf")
os.system("which ic && ic disc.pdf")

fpr_test,tpr_test,_ = roc_curve(y_test, y_pred_test)
fpr_train,tpr_train,_ = roc_curve(y_train, y_pred_train)
auc_test = roc_auc_score(y_test, y_pred_test)
auc_train = roc_auc_score(y_train, y_pred_train)
fig, ax = plt.subplots()
ax.plot(fpr_test,tpr_test, label="test AUC = {:.3f}".format(auc_test))
ax.plot(fpr_train,tpr_train, label="train AUC = {:.3f}".format(auc_train))
ax.set_xlabel("bkg eff")
ax.set_ylabel("sig eff")
ax.set_xlim([0.001,0.8])
ax.set_ylim([0.6,1.0])
ax.set_title("ROC curves")
ax.legend()
fig.set_tight_layout(True)
ax.set_yscale("log")
ax.set_xscale("log")
fig.savefig("roc.pdf")
os.system("which ic && ic roc.pdf")
