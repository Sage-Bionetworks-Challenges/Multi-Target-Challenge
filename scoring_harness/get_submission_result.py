import pandas as pd
import yaml

s = pd.read_csv("submissions_final.csv")
submission_ids = list(s['submissionId'])

score_result = []

for s_id in submission_ids:
    print(s_id)
    file_name = 'files/submission_'+str(s_id)+".yaml"
    print(file_name)
    with open(file_name,"r") as f:
        submission = yaml.load(f)
    for i in xrange(1,3,1):
        submission_p = submission["Problem "+str(i)]
        for index, value in enumerate(submission_p):
            soln = value.keys()[0]
            soln_num = soln.split(' ')[1]
            smiles_str = value[soln]['SMILES string']
            score_result.append({"submission_id":s_id,"problem":i,"solution":soln_num,"smiles_string":smiles_str})

submission_result = pd.DataFrame(score_result)
submission_result = submission_result[['submission_id', 'problem', 'solution', 'smiles_string']]

submission_result.to_csv("submission_final_result.csv",index=False)