#!/bin/python

import pandas as pd

def check_lipinski(row):
    score = 0
    # Lipinski Hbond donor - no more than 5
    score = score + 1 if row['h_bond_donor'] <= 5 else score      
    # Lipinski Hbond acceptor - no more than 10
    score = score + 1 if row['h_bond_acceptor'] <= 10 else score
    # Molecular weight - less than 500
    score = score + 1 if row['moluclar_mass'] < 500 else score
    # Partition coefficient - not greater than 5
    score = score + 1 if row['log_p'] <= 5 else score

    # need to satify 3/4 of the rules
    return(True if score >= 3 else False)

def check_polar_surface_area(num):
    return(True if num < 75 else False)

# lipinski and psa
df_lipinksi_psa = pd.read_csv("lipinski_psa_result.csv")
df_lipinksi_psa['pass_lipinski'] = df_lipinksi_psa.apply(check_lipinski, axis=1)
df_lipinksi_psa['pass_polar_surface_area'] = df_lipinksi_psa['topological_polar_surface_area'].apply(lambda x: check_polar_surface_area(x))
df_lipinksi_psa.to_csv("lipinski_psa_result_final.csv",index=False,encoding="utf-8")

# tanimoto
df_tanimoto = pd.read_csv("tanimoto_result.csv")
target_info = pd.read_csv("chembl/chembl_ids_info.csv")

all_ids = list(set(target_info['chembl_id']))

probelm_group = target_info.groupby('problem')
p1_df = probelm_group.get_group(1)
p1_ids = list(p1_df['chembl_id'])
p1_ids_required = list(p1_df['chembl_id'][p1_df['required']])
p2_df = probelm_group.get_group(2)
p2_ids = list(p2_df['chembl_id'])
p2_ids_required = list(p2_df['chembl_id'][p2_df['required']]) 

result = []
for _, row in df_tanimoto.iterrows():
    row_dict = {}
    row_dict['pass_for_all'] = all(row[all_ids])
    if(row['problem'] == 1):
        row_dict['pass_for_all'] = all(row[p1_ids])
        row_dict['pass_for_required'] = all(row[p1_ids_required])
    else:
        row_dict['pass_for_all'] = all(row[p2_ids])
        row_dict['pass_for_required'] = all(row[p2_ids_required])
    row_dict['submission_id'] = row['submission_id']
    row_dict['smiles_string'] = row['smiles_string']
    result.append(row_dict)

result_df = pd.merge(df_tanimoto,pd.DataFrame(result),on=['submission_id','smiles_string'])
result_df.to_csv("tanimoto_result_final.csv",index=False,encoding="utf-8")
