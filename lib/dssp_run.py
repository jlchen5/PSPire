#!/usr/bin/env python
# Date: 2023-07-17
# Author: Shuang Hou
# Contact: houshuang@tongji.edu.cn


import os
import pandas as pd
from lib.Utility import *
from Bio.PDB.DSSP import make_dssp_dict

# ------------------------------------------------
# get secondary structure from dssp result
# ------------------------------------------------
def dssp_run(uids, files, tmpDir, resume,ofile):
    dic = {'H': 1, 'G': 1, 'I': 1, 'B': 2, 'E': 2, 'T': 0, 'S': 0, '-': 0}
    rows = list()
    destDir = tmpDir + '/dssp_result'
    make_dir(destDir,resume)
    for uid, file in zip(uids, files):
        dsspFile = destDir + '/' + os.path.basename(file).replace('.pdb', '.out')
        if not os.path.isfile(dsspFile):
            run_cmd('mkdssp -i '+file+' -o '+dsspFile)
        else:
            info(f'{file} DSSP result already exists.')
        out = make_dssp_dict(dsspFile)
        dssp = list(out)[0]
        out = [dssp[i][1] for i in dssp.keys()]
        binary = list(pd.Series(out).map(dic))
        structures = ''.join(map(str, binary))
        asa = [str(dssp[i][2]) for i in dssp.keys()]
        aas = [dssp[i][0] for i in dssp.keys()]

        rows.append([uid, structures,','.join(aas),','.join(asa)])

    df = pd.DataFrame(rows, columns=['uniprot_id', 'secondary_structure_binary','aas','dssp_asa_score'])
    df = df.set_index('uniprot_id')
    df.to_pickle(ofile)
    # df.to_csv(destDir + '/dssp_result.csv')
    info('DSSP finished!')
    return df
