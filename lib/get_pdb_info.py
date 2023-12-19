#!/usr/bin/env python
# Date: 2023-07-17
# Author: Shuang Hou
# Contact: houshuang@tongji.edu.cn


from distutils.log import warn
import itertools
import re
import os
import sys

import numpy as np
import pandas as pd

from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1

from lib.Utility import *

sys.path.append('../software')


# ------------------------------------------------
# get IDRs from AlphaFold PDB file or Mobidb-lite
# ------------------------------------------------

def get_seq_and_name(file):
    pdbparser = PDBParser()
    structure = pdbparser.get_structure('aa', file)
    chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
    seq = chains['A']
    return seq,len(seq)

def getIdr_pos(structure, threshold):
    '''position is 1-based'''
    idr = list()
    num = 0
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                num += 1
                if residue.has_id("CA"):
                    ca = residue["CA"]
                    if ca.get_bfactor() < threshold:
                        idr.append(num)
    return idr


def intervals_extract(iterable):
    '''Convert list of sequential number into intervals'''
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]


def matMorphology(seq, rmax=3):
    # One, two and three residues-long ID stretches flanked on both sides by one,
    # two and three residue-long ordered stretches are converted to order and vice- versa.

    # Disorder expansion
    seq = rmax*"D"+seq+rmax*"D"

    for r in range(1, rmax + 1):
        pattern = 'D'*r + 'S'*r + 'D'*r
        new_pattern = 'D'*r + 'D'*r + 'D'*r

        for i in range(0, r + 1):
            seq = seq.replace(pattern, new_pattern)

    # Disorder contraction
    seq = rmax*"S"+seq[rmax:-rmax]+rmax*"S"

    for r in range(1, rmax + 1):
        pattern = 'S'*r + 'D'*r + 'S'*r
        new_pattern = 'S'*r + 'S'*r + 'S'*r

        for i in range(0, r + 1):
            seq = seq.replace(pattern, new_pattern)

    return seq[rmax:-rmax]


def mobidb_run(seq, idr, length, cutoff=10):
    '''Modified algorithm of MobiDB_10-lite to get long disordered regions'''
    state = ['D' if i in idr else 'S' for i in range(1, length+1)]
    state = ''.join(state)

    consensus = matMorphology(state, 3)

    # Structured stretches of up to 10 consecutive residues are then converted to ID
    # if they are flanked by two disordered regions of at least 20 residues
    flag = True
    while flag:
        m = re.search("D{21,}S{1,10}D{21,}", consensus)
        if m:
            matchLength = m.end(0) - m.start(0)
            consensus = consensus[:m.start(
                0)] + "D" * matchLength + consensus[m.end(0):]
        else:
            flag = False

    position = [i for i in np.arange(1, length+1) if consensus[i-1] == 'D']

    idr_intervals = list(intervals_extract(position))
    idr_pos = list()
    # idr_seq = list()
    for i, j in idr_intervals:
        if j-i >= cutoff-1:
            # idr_seq.append(seq[(i-1):j])
            for pos in range(i, j+1):
                idr_pos.append(pos)

    idr_binary = np.zeros(length)
    idr_binary[list(np.array(idr_pos)-1)] = 1
    idr_binary = [str(int(i)) for i in idr_binary]

    # return ''.join(idr_binary), ','.join(idr_seq)
    return ''.join(idr_binary)


def filter_idr(idr,length,seq,secondary_structure_binary):
    idr_binary_ori = np.zeros(length)
    idr_binary_ori[list(np.array(idr)-1)]=1
    
    # idr - sheet - helix
    idr_binary = np.where((idr_binary_ori==1) & (secondary_structure_binary==1),0,idr_binary_ori)
    
    position = [i for i in np.arange(1,length+1) if idr_binary[i-1]==1]
    
    # idr_binary_filter_mobidb_10 = mobidb_run(seq,position,length,10)
    idr_binary_filter_mobidb = mobidb_run(seq,position,length,20)
    
    return idr_binary_filter_mobidb


def get_pdb_info(mobidb_ofile, dsspFile, uids, files, threshold, tmpDir, mobidb):
    if mobidb:
        if os.path.isfile(mobidb_ofile):
            tmp = pd.read_csv(mobidb_ofile,header=None)
            tmp = tmp[0].str.split('\t',expand=True)
            mobidb_data = tmp.loc[tmp[3].isin([None]),:].copy()
        else:
            raise_error(f'The mobidb-lite file is empty: {mobidb_ofile}')
    else:
        if os.path.isfile(dsspFile):
            dssp = pd.read_pickle(dsspFile)
        else:
            raise_error(f'The dssp file is empty: {dsspFile}')
    
    rows = list()
    newFiles = files.copy()
    newUids = uids.copy()
    for uid, file in zip(uids, files):
        seq, length = get_seq_and_name(file)

        if mobidb:
            if os.path.isfile(mobidb_ofile):
                tmp = mobidb_data.loc[mobidb_data[0].str.contains('|'+uid+'|',regex=False),:].copy()
                if tmp.empty:
                    idr_length_filter_mobidb = 0
                    idr_percentage_filter_mobidb = 0
                    idr_binary_filter_mobidb = ''.join(['0']*length)
                else:
                    idr_intervals = [(int(row[1]),int(row[2])) for _,row in tmp.iterrows()]
                    idr_pos = list()
                    for i, j in idr_intervals:
                        for pos in range(i, j+1):
                            idr_pos.append(pos)
                    idr_binary = np.zeros(length)
                    idr_binary[list(np.array(idr_pos)-1)] = 1
                    idr_binary = [str(int(i)) for i in idr_binary]
                    
                    idr_length_filter_mobidb = len(idr_pos)
                    idr_percentage_filter_mobidb = round(idr_length_filter_mobidb/length,3)
                    idr_binary_filter_mobidb = ''.join(idr_binary)
            else:
                idr_length_filter_mobidb = 0
                idr_percentage_filter_mobidb = 0
                idr_binary_filter_mobidb = ''.join(['0']*length)
        else:
            structure_binary = dssp.loc[uid, 'secondary_structure_binary'].replace('2', '1')
            structure_binary = np.array(list(map(int, list(structure_binary))))
            if len(structure_binary) != length:
                warn(f'SeqIO package can not get true number of amino acids from {file}: (1)seqIO length: {length}; (2)dssp length: {len(structure_binary)}. Skip this file.')
                newFiles.remove(file)
                newUids.remove(uid)
                continue
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('_', file)
            idr = getIdr_pos(structure, threshold)
            idr_binary_filter_mobidb = filter_idr(idr,length,seq,structure_binary)
            idr_length_filter_mobidb = idr_binary_filter_mobidb.count('1')
            idr_percentage_filter_mobidb = round(idr_length_filter_mobidb/length,3)

        rows.append([uid, seq, length, idr_length_filter_mobidb,
                    idr_percentage_filter_mobidb, idr_binary_filter_mobidb])

    if len(newFiles)==0:
        warn('There is no valid id to analyze.')
        sys.exit(1)

    df = pd.DataFrame(rows, columns=['uniprot_id', 'sequence', 'length',
                                     'idr_length_filter_mobidb', 'idr_percentage_filter_mobidb', 'idr_binary_filter_mobidb'])
    df = df.set_index('uniprot_id')
    df.to_pickle(tmpDir + '/pdbInfo.pkl')
    info('Get PDB info finished!')
    return df, newFiles, newUids
