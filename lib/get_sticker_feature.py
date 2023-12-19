#!/usr/bin/env python
# Date: 2023-07-17
# Author: Shuang Hou
# Contact: houshuang@tongji.edu.cn

"""
Description: calculate sticker features
"""

import numpy as np
from lib.Utility import *
import scipy.cluster.hierarchy as hcluster
from pymol import cmd
from Bio.PDB import PDBParser


def cal_sticker_feature(numbers,numbers_neg,numbers_ssup,pdbfile,ssup_length):
    thresh = 14
    cluster_method = 'centroid'
    cmd.load(pdbfile,'pdb')
    cmd.select("aa", "/pdb//A/*/CA")
    score = list()
    for num in numbers_ssup:
        cmd.select("bb", "aa within "+str(thresh)+" of /pdb//A/"+str(num)+"/CA")
        ids = []
        cmd.iterate("bb", lambda atom: ids.append(atom.resi))
        ids = [int(i) for i in ids]
        ids = set(ids) & set(numbers_ssup)
        pos_num = len(set(numbers) & set(ids))
        neg_num = len(set(numbers_neg) & set(ids))
        score.append(pos_num - neg_num)
    
    cmd.delete('pdb')
    
    score = np.array(score)
    if np.sum((score>=3) | (score<=-3))==0:
        return 0,0,'',''
    
    parser = PDBParser()
    structure = parser.get_structure('test', pdbfile)
    model = structure[0]
    chain = model['A']
    aas = list(chain)
    
    numbers_ssup = np.array(numbers_ssup)
    filter_numbers_ssup = numbers_ssup[(score>=3) | (score<=-3)]
    
    coords = [aas[i-1]['CA'].coord for i in filter_numbers_ssup]
    coords = np.vstack(coords)
    
    if np.sum((score>=3) | (score<=-3))==1:
        clusters = np.array([1])
    else:
        clusters = hcluster.fclusterdata(coords, int(thresh), criterion="distance",metric='euclidean',method=cluster_method)
    
    filter_score = score[(score>=3) | (score<=-3)]
        
    result = list()
    for i in range(1,np.max(clusters)+1):
        scores_tmp = filter_score[clusters==i]
        if np.sum(scores_tmp>0) > np.sum(scores_tmp<0):
            charge_type = 'pos'
        elif np.sum(scores_tmp>0) < np.sum(scores_tmp<0):
            charge_type = 'neg'
        else:
            charge_type = 'both'
        
        if charge_type == 'pos':
            tmp = np.max(scores_tmp)
        elif charge_type == 'neg':
            tmp = np.min(scores_tmp)
        else:
            tmp = 0
            
        result.append(tmp)
        
    result = np.array(result)
    pos_num = np.sum(result>0)
    neg_num = np.sum(result<0)
    
    cluster_positions = [tuple(np.array(filter_numbers_ssup)[clusters==i]) for i in range(1,np.max(clusters)+1)]
    cluster_scores = [tuple(filter_score[clusters==i]) for i in range(1,np.max(clusters)+1)]
    pos_positions = list()
    neg_positions = list()
    cluster_positions = cluster_positions
    for num,scores_tmp in enumerate(cluster_scores):
        if len(scores_tmp)==1:
            if int(scores_tmp[0])>0:
                pos_positions.append(cluster_positions[num])
            else:
                neg_positions.append(cluster_positions[num])
        else:
            scores_tmp = np.array(scores_tmp)
            if np.sum(scores_tmp>0) > np.sum(scores_tmp<0):
                pos_positions.append(cluster_positions[num])
            elif np.sum(scores_tmp>0) < np.sum(scores_tmp<0):
                neg_positions.append(cluster_positions[num])
    
    min_sticker_num = min(pos_num,neg_num)
    sticker_num = len(result)
    
    min_sticker_num_percent = min_sticker_num/ssup_length
    sticker_num_percent = sticker_num/ssup_length
    
    return min_sticker_num_percent,sticker_num_percent,str(pos_positions),str(neg_positions)

