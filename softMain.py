#for simplicities sake, I'm going to be using the exact data set
#and code used in tutorial 13, which contains the pdbbind data

#start with common imports
import os
import numpy as np
import pandas as pd
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
import deepchem as dc
from simtk.openmm.app import PDBFile
from pdbfixer import PDBFixer

from deepchem.utils.vina_utils import prepare_inputs
from deepchem.utils import download_url, load_from_disk
from deepchem.models import AtomicConvModel
from deepchem.feat import RdkitGridFeaturizer

#slight changes are being made to follow my normal coding structure
#create a data directory that will contain our file
dataDir = dc.utils.get_data_dir()
datasetFile = os.path.join(dataDir, "pdbbind_core_df.csv.gz")

#if the file location doesnt exist, the program will find and
#download the file from designated url
if not os.path.exists(datasetFile):
    print('File does not exist. Downloading file...')
    download_url("https://s3-us-west-1.amazonaws.com/deepchem.io/datasets/pdbbind_core_df.csv.gz")
    print('File downloaded...')

#place the data file into a program "understood" format
rawDataset = load_from_disk(datasetFile)
rawDataset = rawDataset[['pdb_id', 'smiles', 'label']]

#visualize a little bit for sanity's sake
rawDataset.head()

#pull the data set into pieces that allow us to manipulate it
pdbid = rawDataset['pdb_id'].iloc[25]
ligand = rawDataset['smiles'].iloc[25]

#now apply a fixer to the pdbbind data. sets a fixer and writes to pdb file
fixer = PDBFixer(pdbid=pdbid)
PDBFile.writeFile(fixer.topology, fixer.positions, open('%s.pdb' % (pdbid), 'w'))

#some generic stuff from the pdbbind fixer. still working on understanding it,
# but i might just accept the outputs without understanding it
p, m = None, None
# fix protein, optimize ligand geometry, and sanitize molecules
try:
    p, m = prepare_inputs('%s.pdb' % (pdbid), ligand)
except:
    print('%s failed PDB fixing' % (pdbid)) 

if p and m:  # protein and molecule are readable by RDKit
    print(pdbid, p.GetNumAtoms())
    Chem.rdmolfiles.MolToPDBFile(p, '%s.pdb' % (pdbid))
    Chem.rdmolfiles.MolToPDBFile(m, 'ligand_%s.pdb' % (pdbid))

#Here we split the data into two variables, proteins and ligands, for input into featurizer
pdbids = rawDataset['pdb_id'].values
ligand_smiles = rawDataset['smiles'].values

proteins = [f for f in os.listdir('.') if len(f) == 8 and f.endswith('.pdb')]
ligands = [f for f in os.listdir('.') if f.startswith('ligand') and f.endswith('.pdb')]

small_dataset = rawDataset[rawDataset['pdb_id'].isin(pdbids)]
labels = small_dataset.label

#Now begin Featurizing
#first we'll start with the RDKIT grid featurizer
fp_featurizer = dc.feat.RdkitGridFeaturizer(size=2048)
fp_features = fp_featurizer.featurize(zip(ligands, proteins))

#now initialize the Weave Featurizer
wv_featurizer = dc.feat.WeaveFeaturizer()
wv_features = wv_featurizer.featurize(zip(ligands,proteins))

#Activation of the machine learning algorithim
dataset = dc.data.NumpyDataset(X=features, y=labels, ids=pdbids)
train_dataset, test_dataset = dc.splits.RandomSplitter().train_test_split(dataset, seed=42)
