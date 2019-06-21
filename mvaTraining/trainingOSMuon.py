#!/usr/bin/env python

from ROOT import TMVA, TFile, TTree, TCut
from subprocess import call
from os.path import isfile
import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'

import numpy as np
import keras
import h5py

from keras.models import Sequential
from keras.layers import Dense, Activation, Conv2D, MaxPooling2D, Flatten, Dropout
from keras.regularizers import l2
from keras.optimizers import SGD, Adam
from keras.callbacks import ModelCheckpoint

##### FUNCTIONS

def getKerasModel(inputDim, modelName, nLayers = 3, layerSize = 200, dropValue = 0.2):
    model = Sequential()
    model.add(Dense(layerSize, activation='relu', kernel_initializer='normal', input_dim=inputDim))
    if dropValue != 0:
        model.add(Dropout(dropValue))

    for i in range(1, nLayers):
        model.add(Dense(layerSize, activation='relu', kernel_initializer='normal'))
        if dropValue != 0:
            model.add(Dropout(dropValue))

    # Anything below this point should not be changed in order for the network to work with TMVA, exception: the optimizer
    model.add(Dense(2, activation='softmax'))

    opt = Adam(lr=0.001)
    model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])
    model.save(modelName)
    model.summary()
    return modelName


# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

# Load data
file = '../ntuples/ntuBsDG0MC2018.root'

data = TFile.Open(file)

tree = data.Get('PDsecondTree')

# Prepare factory
name = 'OsMuonHLTJpsiMu'

outputName = 'TMVA' + name + '.root'

output = TFile.Open(outputName, 'RECREATE')

factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification')

dataloader = TMVA.DataLoader('dataset')

# variable list
varList = [
    ('muoPt', 'F')
    ,('muoEta', 'F')
    ,('muoDxy', 'F')
    ,('muoExy', 'F')
    ,('muoDz', 'F')
    ,('muoEz', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDrB', 'F')
    ,('muoPFIso', 'F')
    ,('muoConePt', 'F')
    ,('muoConePtRel', 'F')
    ,('muoConeDr', 'F')
    ,('muoConeEnergyRatio', 'F')
    ,('muoConeQ', 'F')
    ,('muoConeCF', 'F')
    ,('muoConeNCH', 'I')
    # ,('muoJetDFprob', 'F')
    ]

varListClean = [
    ('muoPt', 'F')
    ,('muoEta', 'F')
    ,('muoDxy', 'F')
    ,('muoExy', 'F')
    ,('muoDz', 'F')
    ,('muoEz', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDrB', 'F')
    ,('muoPFIso', 'F')
    ,('muoConeCleanPt', 'F')
    ,('muoConeCleanPtRel', 'F')
    ,('muoConeCleanDr', 'F')
    ,('muoConeCleanEnergyRatio', 'F')
    ,('muoConeCleanQ', 'F')
    ]

varList = varListClean

# automatic variable counting and adding
nVars = 0
for var in varList:
    dataloader.AddVariable( var[0], var[1] )
    nVars += 1

# prepare dataloader
# Event wise selection
cut = 'osMuon==1'

cutSgn = cut + '&&osMuonTag==1' #Correcly tagged events selection i.e. sign(charge lepton) -> correct flavour 
cutBkg = cut + '&&osMuonTag==0' #Uncorrecly tagged events

# same tree, add selection later
dataloader.AddSignalTree(tree)
dataloader.AddBackgroundTree(tree)

# evtWeight variable in the ntuple, address simulation bias
dataloader.SetWeightExpression( 'evtWeight' );

# nBkg = 169226
# nSgn = 393356

nBkg = '126920'
nSgn = '295017'
nBkgTest = '42306'
nSgnTest = '98339'

dataloaderOpt = 'nTrain_Signal=' + nSgn + ':nTrain_Background=' + nBkg + ':nTest_Signal=' + nSgnTest + ':nTest_Background=' + nBkgTest
dataloaderOpt += ':SplitMode=Random:NormMode=NumEvents:V:'

dataloader.PrepareTrainingAndTestTree(TCut(cutSgn), TCut(cutBkg), dataloaderOpt)

# Create Keras Model
nLayers = 3
layerSize = 200
dropValue = 0.4

modelName = getKerasModel(nVars, 'model' + name + '.h5', nLayers, layerSize, dropValue)
# modelName = 'TrainedModel_DNNOsMuonHLTJpsiMu.h5'
# Book methods
dnnOptions = '!H:!V:NumEpochs=50:TriesEarlyStopping=10:BatchSize=1024:ValidationSize=33%:SaveBestOnly=True'
dnnOptions = dnnOptions + ':Tensorboard=./logs:FilenameModel=' + modelName
# 

# Preprocessing string creator, loop was for selection of which variable to apply gaussianification
preprocessingOptions = ':VarTransform=N,G,N'

dnnName = 'DNN' + name

factory.BookMethod(dataloader, TMVA.Types.kPyKeras, dnnName, dnnOptions + preprocessingOptions)

# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
