#!/bin/python

import rdkit
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rd


def main():
    sub_df = pd.read_csv("submissions_final_result.csv")

    cmp_ds  = []

    for _,row in sub_df.iterrows():
        cmp_dict = {}
        mol = Chem.MolFromSmiles(row['smiles_string'])
        cmp_dict['submission_id'] = row['submission_id']
        cmp_dict['smiles_string'] = row['smiles_string']
        
        # Lipinski's rule
        cmp_dict['h_bond_donor'] = rd.CalcNumLipinskiHBD(mol)     # Lipinski Hbond donor
        cmp_dict['h_bond_acceptor'] = rd.CalcNumLipinskiHBA(mol)  # Lipinski Hbond acceptor
        cmp_dict['moluclar_mass'] = rd._CalcMolWt(mol)            # Molecular Weight
        cmp_dict['log_p'] = rd.CalcCrippenDescriptors(mol)[0]     # Partition coefficient
        
        # Topological polar surface area
        cmp_dict['topological_polar_surface_area'] = rd.CalcTPSA(mol)

        cmp_ds.append(cmp_dict)

    result = pd.merge(sub_df,pd.DataFrame(cmp_ds),on=['submission_id','smiles_string'])
    result.to_csv("lipinski_psa_result.csv",index=False,encoding='utf-8')

if __name__ == '__main__':
    main()