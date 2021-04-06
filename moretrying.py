import os
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

#if no such file exists, make the file. This should make a ne file in the 
#currently active directory
if not os.path.exists(os.path.dirname(os.path.realpath(__file__)) +'\PDBFiles'):
    #start by making a new file in the current directory
    print('Creating file directory...')
    os.makedirs(os.path.dirname(os.path.realpath(__file__)) +'\PDBFiles')
else:  
    #change to the correct directory
    os.chdir(os.path.dirname(os.path.realpath(__file__)) +'\PDBFiles')

#Now Check for the dataset being downloaded
if not os.path.exists(os.path.dirname(os.path.realpath(__file__)) +'\pdbbind_v2019_other_PL.tar.gz'):
  wholeSet = input('Would you Like to download whole pdbbind Dataset? \n(Y/N)\n')
  if wholeSet.casefold() == 'y':
    #Now download the pdbbind compressed dataset to target location
    print(os.path.dirname(os.path.realpath(__file__)) +'\pdbbind_v2019_other_PL.tar.gz')
    print('Missing dataset. Downloading file. This might take a minute...')
    download_url("http://www.pdbbind-cn.org/download/pdbbind_v2019_other_PL.tar.gz")
    print('File Downloaded!')
  else:
    pass

if not os.path.exists(os.path.dirname(os.path.realpath(__file__)) +'\pdbbind_core_df.csv.gz'):
  coreSet = input('Would you Like to download core pdbbind Dataset? \n(Y/N)\n')
  if coreSet.casefold() == 'y':
    print('Downloading file...')
    download_url("https://s3-us-west-1.amazonaws.com/deepchem.io/datasets/pdbbind_core_df.csv.gz", dest_dir = os.getcwd())
    print('File downloaded!')
  else:
    pass

#lets look for active directories
direc = os.getcwd() #current working directory
fileDict = {}
counter = 1
for filDir in os.listdir(direc):
  fileDict[str(counter)] = filDir
  counter += 1

#move to another directory depending on the dataset
for keys in fileDict:
  print(keys, ':', fileDict[keys])
dataFile = input('Which Dataset are you working with?\n')
print(fileDict[dataFile])
os.chdir(fileDict[dataFile])

#now grab the name of each file in the directory
direc = os.getcwd()
for pdbid in os.listdir(direc):
    fixer = PDBFixer(url='https://files.rcsb.org/download/%s.pdb' % (pdbid))
    PDBFile.writeFile(fixer.topology, fixer.positions, open('%s.pdb' % (pdbid), 'w'))