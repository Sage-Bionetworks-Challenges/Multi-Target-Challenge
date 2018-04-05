import pandas as pd
from chembl_webresource_client.new_client import new_client

target = new_client.target
activity = new_client.activity

target_info = pd.read_csv('../chembl/chembl_ids_info.csv')

for _, row in target_info.iterrows():
    print(row['chembl_id'])
    compound_query = activity.filter(target_chembl_id=row['chembl_id'],standard_type__iregex='(IC50|Kd)')

    compound_info = pd.DataFrame(list(compound_query.all()))
    compound_info[['standard_value']] = compound_info[['standard_value']].astype(float)
    if(row['bind']):
        compound_info = compound_info.loc[(compound_info['standard_value']<10000) & compound_info['standard_relation'].isin(['<','<=','='])]
    else:
        compound_info = compound_info.loc[(compound_info['standard_value']>=10000) & compound_info['standard_relation'].isin(['=','>','>='])]
    
    cpd_file_name = "../chembl/"+row['chembl_id']+"_ligands.tsv"
    compound_info.to_csv(cpd_file_name,sep='\t',index=False,encoding='utf-8')