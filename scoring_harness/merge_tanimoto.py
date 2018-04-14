#!/bin/python

import pandas as pd
import glob

def merge_files(path,pattern):
    allFiles = glob.glob(path + pattern)
    dfList = []
    for file_ in allFiles:
        df = pd.read_csv(file_,index_col=None)
        dfList.append(df)
    frame = reduce(lambda x, y: pd.merge(x,y,on =['submission_id','smiles_string','problem','solution']), dfList)
    return frame

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Merge files")
    parser.add_argument('path', help='path to files for merging')
    parser.add_argument('pattern', help='files pattern')
    parser.add_argument('output', help='output file; CSV format')

    args = parser.parse_args()

    frame = merge_files(args.path,args.pattern)
    frame.to_csv(args.output,index=False)

if __name__ == '__main__':
    main()