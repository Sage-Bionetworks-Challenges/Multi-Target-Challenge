## Open files

import re,gc
from CommonUtility import file_handle
from rdkit import Chem

def rdkit_open(File_Tulup):
  List = []
  for f in File_Tulup:
    handle = file_handle(f)
    if re.search(r'.sdf', f):
      Mol = [x for x in Chem.ForwardSDMolSupplier(handle, removeHs=False)
             if x is not None]
    if re.search(r'.smi', f):
      Mol = [x for x in Chem.SmilesMolSupplier(f, titleLine=True,
                                                  delimiter=' |\t|,')
             if x is not None]
    if re.search(r'.mol2',f):
      Mol = [x for x in Chem.MolFromMol2File(f, removeHs=False)
             if x is not None]
    print "# Found mol in {0}: {1}".format(f, len(Mol))
    for mol in Mol: List.append(mol)
 
  gc.collect()
  return List

