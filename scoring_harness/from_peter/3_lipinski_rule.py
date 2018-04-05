#!/usr/bin/python

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    17.12.12    
#
#   Calculate molecular properties via RDKit
#
###########################################################################

import os,sys

msg = '''\n\t{0}\n\t\t[library of ligands smi|sdf]
\t\t[Output .csv/.xlsx prefix]
\n\te.g.>\t x.py mol.smi mol_result\n'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

import rdkit
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rd

from rdkit_open import *
from CommonUtility import *

##########################################################################

def main(in_file, output):

  Cmpds  = {}
  InMols = rdkit_open([in_file])
  print('\n # Number of input molecule: {0}'.format(len(InMols)))
  for mol in InMols:
    m = {}

    name = mol.GetProp('_Name').split()[0]
    
    m['Name'] = name
    m['Formula'] = rd.CalcMolFormula(mol)
    m['SMILES'] = Chem.MolToSmiles(mol)

    m['MW']   = rd._CalcMolWt(mol)               # Molecular Weight
    m['logP'] = rd.CalcCrippenDescriptors(mol)[0]  # Partition coefficient
    m['HDon'] = rd.CalcNumLipinskiHBD(mol)      # Lipinski Hbond donor
    m['HAcc'] = rd.CalcNumLipinskiHBA(mol)      # Lipinski Hbond acceptor
    m['TPSA'] = rd.CalcTPSA(mol)                # Topological polar surface area

    m['Rotat'] = rd.CalcNumRotatableBonds(mol, strict=True) # Rotatable bond
    m['MolRef'] = rd.CalcCrippenDescriptors(mol)[1]         # Molar refractivity
    m['AliRing'] = rd.CalcNumAliphaticRings(mol)        # Aliphatic ring number
    m['AroRing'] = rd.CalcNumAromaticRings(mol)         # Aromatic ring number
#    m['Stereo'] = rd.CalcNumAtomStereoCenters(mol)      # Stereo center number
#    m['UnspStereo'] = rd.CalcNumUnspecifiedAtomStereoCenters(mol)  # unspecified stereo

    m['SMILES'] = Chem.MolToSmiles(mol, 
                    isomericSmiles=True, allHsExplicit=False)
    Cmpds[name] = m

  ####################################

  df = pd.DataFrame.from_dict(Cmpds, orient='index')
  df.index.name = 'Name'

  # Columns of data to print out
  Columns = [ 'Formula',
              'MW',    'logP',   'HDon',    'HAcc',    'TPSA',
              'Rotat', 'MolRef', 'AliRing', 'AroRing', 
              #'Stereo', 'UnspStereo', 
              'SMILES', ]
  reorder = df[Columns]

  # Output to CSV
  reorder.to_csv( output+'.csv', sep=',', na_rep='NA', encoding='utf-8',
                  float_format='%.5f', header=True )

  # Output to Excel
  reorder.to_excel( output+'.xlsx', header=True, na_rep='NA' )


##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2] )
