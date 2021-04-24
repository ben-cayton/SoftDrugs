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

#if no such file exists, make the file. This should make a ne file in the 
#currently active directory
if not os.path.exists(os.path.dirname(os.path.realpath(__file__)) +'\PDBFiles'):
    #start by making a new file in the current directory
    print('Creating file directory...')
    os.makedirs(os.path.dirname(os.path.realpath(__file__)) +'\PDBFiles')
else:  
    #change to the correct directory
    os.chdir(os.path.dirname(os.path.realpath(__file__)) +'\PDBFiles')

#get the pdb files directory for later use
pdbDirec = os.getcwd()

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

#now grab the name of each file in the directory and
#save the current working directory before i mess with it!
direc = os.getcwd()
for pdbid in os.listdir(direc):
  #print(pdbid)
  if not direc == os.getcwd():
    os.chdir(direc)
  else:
    pass

  p, m = None, None
  try:
    if not os.path.exists(os.path.dirname(os.path.realpath(__file__)) + '\%s' % (pdbid)):
      fixer = PDBFixer(url='https://files.rcsb.org/download/%s.pdb' % (pdbid))
      PDBFile.writeFile(fixer.topology, fixer.positions, open('%s.pdb' % (pdbid), 'w'))
    #write all of the ligand files as well
    #WE WANT SMILES STUFF BEFORE  WE SAVE THE FILE
    os.chdir(direc+'/%s' % (pdbid))
    ligand = None
    for eachfile in os.listdir(os.getcwd()):
      if eachfile == ('%s_ligand.mol2' % (pdbid)):
        molfile = Chem.rdmolfiles.MolFromMol2File(eachfile)
        ligand = Chem.rdmolfiles.MolToSmiles(molfile)
        #Chem.rdmolfiles.MolToPDBFile(molfile, 'ligand_%s.pdb' % (pdbid))
        #shutil.move(os.getcwd()+'\ligand_%s.pdb' % pdbid, os.path.dirname(os.getcwd()))
      else:
        continue
    #return to the previous directory
    os.chdir(direc)

    if not os.path.exists(os.path.dirname(os.path.realpath(__file__)) + '\Processedfiles.txt'):
      newfile = open('Processedfiles.txt', 'w')
      newfile.close()

    #now run the full fixing step
    try:
      newFile = open('Processedfiles.txt', 'r')
      if str(pdbid) not in newFile.read():
        newFile.close()
        print("Sanitizing %s" %(pdbid))
        proFile = open('Processedfiles.txt', 'a')
        proFile.write(str(pdbid+'\n'))
        proFile.close()
        p, m = prepare_inputs('%s.pdb' % (pdbid), ligand)
        newFile.close()
      else:
        continue
    except:
      print('%s failed sanitization' % (pdbid))
  except:
    continue
    
  if p and m:  # protein and molecule are readable by RDKit
    #Start in the directory of pdb files
    os.chdir(pdbDirec)

    #make a new file called 'sanitized files"and switch to it
    if not os.path.exists(os.path.dirname(os.path.realpath(__file__)) +'\SanitizedFiles'):
      os.makedirs(os.path.dirname(os.path.realpath(__file__)) +'\SanitizedFiles')
    else:
      pass
    os.chdir(os.getcwd()+'\SanitizedFiles')

    #check if there is already a sanitized file present. if not, make one. if so, skip it
    if not (os.path.exists(os.path.dirname(os.path.realpath(__file__)) + '\%s.pdb' % (pdbid)) and os.path.exists(os.path.dirname(os.path.realpath(__file__)) + '\ligand_%s.pdb' % (pdbid))):
      Chem.rdmolfiles.MolToPDBFile(p, '%s.pdb' % (pdbid))
      Chem.rdmolfiles.MolToPDBFile(m, 'ligand_%s.pdb' % (pdbid))
