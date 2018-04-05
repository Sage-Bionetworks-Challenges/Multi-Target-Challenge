#!/bin/python

import pandas
from rdkit import Chem
from rdkit import DataStructs

from rdkit.Chem import AllChem

def get_mol(smiles_string):
    return(Chem.MolFromSmiles(smiles_string))

def get_fingerprint(mol):
    return(AllChem.GetMorganFingerprintAsBitVect(mol,2,2048))

def compare_similarity(fp1, fp2):
    return(DataStructs.FingerprintSimilarity(fp1,fp2))

def _get_result(smiles, ref_df):
    result = {}
    mol = get_mol(smiles)
    fp = get_fingerprint(mol)

    for _, row in ref_df.iterrows():
        ref_mol = get_mol(row['canonical_smiles'])
        ref_fp = get_fingerprint(ref_mol)
        coeff = compare_similarity(fp,ref_fp)
        result[row['molecule_chembl_id']] = coeff
    return result

def get_result(sub_df, ref_df, target_name):
    print("calculating...")
    result = []
    for _, row in sub_df.iterrows():
        smiles_str = row['smiles_string']
        row_cal_dict = _get_result(smiles_str,ref_df)
        row_dict = {}
        row_dict[target_name] = all(v < 0.4 for v in row_cal_dict.values())
        row_dict['submission_id'] = row['submission_id']
        row_dict['smiles_string'] = row['smiles_string']
        result.append(row_dict)
    return result

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Calculate tanimoto")
    parser.add_argument('submission_file', help='submission file, txt/tsv format')
    parser.add_argument('reference_file', help='reference_file, txt/tsv format')
    parser.add_argument('target_name', help='target_name, use CHEMBL ID')

    args = parser.parse_args()

    submission_file = args.submission_file
    reference_file = args.reference_file
    target_name = args.target_name
    
    print("target: "+target_name)

    sub_df = pandas.read_csv(submission_file)
    
    ref_df = pandas.read_table(reference_file)
    ref_df = ref_df[['molecule_chembl_id','canonical_smiles']]
    ref_df = ref_df.drop_duplicates()

    result_ls = get_result(sub_df,ref_df,target_name)
    result_df = pandas.merge(sub_df,pandas.DataFrame(result_ls),on=['submission_id','smiles_string'])
    result_df.to_csv(target_name+"_tanimoto_result.csv",index=False)

if __name__ == '__main__':
    main()
        
