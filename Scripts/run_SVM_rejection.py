import os
from sys import argv
from pathlib import Path

import numpy as np
import pandas as pd
import time as tm
from sklearn.svm import LinearSVC
import rpy2.robjects as robjects
from sklearn.calibration import CalibratedClassifierCV


def run_SVM_rejection(DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath = "", NumGenes = 0, Threshold = 0.7):

    # read the Rdata file
    robjects.r['load'](CV_RDataPath)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1
    test_ind = np.array(robjects.r['Test_Idx'])
    train_ind = np.array(robjects.r['Train_Idx'])

    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)

    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]

    # read the feature file
    if (NumGenes > 0):
        features = pd.read_csv(GeneOrderPath,header=0,index_col=None, sep=',')

    # normalize data
    data = np.log1p(data)

    Classifier = LinearSVC()
    clf = CalibratedClassifierCV(Classifier)

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        if (NumGenes > 0):
            feat_to_use = features.iloc[0:NumGenes,i]
            train = train.iloc[:,feat_to_use]
            test = test.iloc[:,feat_to_use]

        start=tm.time()
        clf.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = clf.predict(test)
        prob = np.max(clf.predict_proba(test), axis = 1)
        unlabeled = np.where(prob < Threshold)
        predicted[unlabeled] = 'Unknown'
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    OutputDir = Path(OutputDir)
    truelab.to_csv(str(OutputDir / Path("SVM_rejection_true.csv")),
                   index = False)
    pred.to_csv(str(OutputDir / Path("SVM_rejection_pred.csv")),
                index = False)
    tr_time.to_csv(str(OutputDir / Path("SVM_rejection_training_time.csv")),

                   index = False)
    ts_time.to_csv(str(OutputDir / Path("SVM_rejection_test_time.csv")),
                   index = False)

run_SVM_rejection(argv[1], argv[2], argv[3], argv[4], argv[5], int(argv[6]))
