#here's all of the imports, etc. not all of these are needed, but i don't
#want to sift through them all to find the ones i need to keep

import os
import shutil
import deepchem as dc
import numpy as np
import pandas as pd
from deepchem.feat import RdkitGridFeaturizer
from deepchem.models import AtomicConvModel
from deepchem.utils import download_url, load_from_disk
from deepchem.utils.vina_utils import prepare_inputs
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem
from simtk.openmm.app import PDBFile
from sklearn.ensemble import RandomForestRegressor

from deepchem.utils.vina_utils import prepare_inputs
from deepchem.utils import download_url, load_from_disk
from deepchem.models import AtomicConvModel
from deepchem.feat import RdkitGridFeaturizer
from sklearn.ensemble import RandomForestRegressor
from deepchem.utils.evaluate import Evaluator

#Find the sanitized molecules file and store all pairs to a dictionary
#the file should exist in the current directory, but further down
protLigDict = {}
progDirec = os.getcwd()
for eachfile in os.listdir(progDirec):
    if eachfile == 'PDBFiles':
        os.chdir(progDirec+'\PDBFiles')
        pdbDirec = os.getcwd()
    else:
        continue

for eachfile in os.listdir(pdbDirec):
    if eachfile == 'SanitizedFiles':
        os.chdir(pdbDirec+'\SanitizedFiles')
        sanitDirec = os.getcwd()
    else:
        continue

#now that I found the files, I will parse through them and match the ligand protein
#combinations before storing them in a dictionary
for eachfile in os.listdir(sanitDirec):
    if 'ligand_%s.pdb' % eachfile[:4] in os.listdir(sanitDirec):
        protLigDict[eachfile[:4]] = ['ligand_%s.pdb' % eachfile[:4], eachfile]
    else:
        continue

for keys in protLigDict:
  print(keys, ':', protLigDict[keys])

#save the data to text documents, there's too much to put into
#a dictionary
#start by making a protein ligand file
protLigFile = open('proteinLigandSMILES.txt', 'a')

for prot in protLigDict:
    #start with making smiles data of proteins
    molProt = Chem.rdmolfiles.MolFromPDBFile(protLigDict[prot][1])
    protSmiles = Chem.rdmolfiles.MolToSmiles(molProt)
    
    #grab the dicts items and pull/ make SMILES data for ligands
    protLigItem = protLigDict[prot]
    molLig = Chem.rdmolfiles.MolFromPDBFile(protLigItem[0])
    ligSmiles = Chem.rdmolfiles.MolToSmiles(molLig)

    #store each to the created appended file
    protLigFile.write(protSmiles+','+ligSmiles+'\n')
    
#close the file for sanitary sake
protLigFile.close()

#here are all of our featurizers. I'll push this down as coding demands adding in
#new lines as needed

#rdkitdescriptor --> smiles iterable, circularfingerprint --> smiles iterable, MACCSKeysFingerprint --> rdkitmol object
protLigFile = open('proteinLigandSMILES.txt', 'r')
#ligFeats = open('FeaturizedLigands.txt','a')
for eachline in protLigFile:
    molLig = Chem.rdmolfiles.MolFromSmiles(eachline(-1))
    rd_featurizer = dc.feat.RDKitDescriptors(size=2048)
    rd_features = rd_featurizer.featurize(molLig)
   
    cf_featurizer = dc.feat.CircularFingerprint(size=2048)
    cf_features = cf_featurizer.featurize(molLig)

    mk_featurizer = dc.feat.MACCSKeysFingerprint(size=2048)
    mk_features = mk_featurizer.featurize(molLig)

    pc_featurizer = dc.feat.PubChemFingerprint(size=2048)
    pc_features = pc_featurizer.featurize(molLig)

'''
for prot in protLigDict:
    
    rd_featurizer = dc.feat.RdkitGridFeaturizer(size=2048)
    rd_features = rd_featurizer.featurize(protLigDict[prot])
   
    cp_featurizer = dc.feat.ContactCircularFingerprint(size=2048)
    cp_features = cp_featurizer.featurize(protLigDict[prot])

    ac_featurizer = dc.feat.AtomicConvFeaturizer(size=2048)
    ac_features = ac_featurizer.featurize(protLigDict[prot])
'''