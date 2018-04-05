#!/bin/bash

for filename in chembl/*_ligands.tsv; do
    echo "reading file - $filename"
    new_file=${filename/_ligands.tsv/}
    new_file=${new_file/chembl\//}
    python calculate_tanimoto.py submissions_final_result.csv $filename $new_file
done

echo "merging files..."
python merge_tanimoto.py

echo "completed."