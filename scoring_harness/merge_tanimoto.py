#!/bin/python

import pandas as pd
import glob

def main():
    path =r'.'
    allFiles = glob.glob(path + "/*_tanimoto_result.csv")
    dfList = []
    for file_ in allFiles:
        df = pd.read_csv(file_,index_col=None)
        dfList.append(df)
    frame = reduce(lambda x, y: pd.merge(x,y,on =['submission_id','smiles_string','problem','solution']), dfList)

    frame.to_csv('tanimoto_result.csv',index=False)

if __name__ == '__main__':
    main()