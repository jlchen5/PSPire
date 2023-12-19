#!/usr/bin/env python
# Date: 2023-07-17
# Author: Shuang Hou
# Contact: houshuang@tongji.edu.cn


"""
Description: calculate final features
"""

import os
import itertools

import pandas as pd
import numpy as np

from lib.Utility import *
from lib.get_sticker_feature import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# ------------------------------------
# static parameters
# ------------------------------------
RESIDUES = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

residue_key = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 
               'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE', 
               'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 
               'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER', 
               'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}

residue_key2 = {residue_key[i]:i for i in residue_key.keys()}

dssp_standard_asa = {"A": 106.0, "R": 248.0, "N": 157.0, "D": 163.0,
                "C": 135.0, "Q": 198.0, "E": 194.0, "G": 84.0,
                "H": 184.0, "I": 169.0, "L": 164.0, "K": 205.0,
                "M": 188.0, "F": 197.0, "P": 136.0, "S": 130.0,
                "T": 142.0, "W": 227.0, "Y": 222.0, "V": 142.0}

# AAindex: GRAR740102
Polarity = {'A':8.1, 'L':4.9, 'R':10.5, 'K':11.3,
            'N':11.6, 'M':5.7, 'D':13.0, 'F':5.2,
            'C':5.5, 'P':8.0, 'Q':10.5, 'S':9.2,
            'E':12.3, 'T':8.6, 'G':9.0, 'W':5.4,
            'H':10.4, 'Y':6.2, 'I':5.2, 'V':5.9}

psaia_col = ['chain_id','ch_total_ASA','ch_b_bone_ASA','ch_s_chain_ASA','ch_polar_ASA','ch_n_polar_ASA',
            'res_id','res_name','total_ASA','b_bone_ASA','s_chain_ASA','polar_ASA','n_polar_ASA',
            'total_RASA','b_bone_RASA','s_chain_RASA','polar_RASA','n_polar_RASA','average_DPX','s_avg_DPX',
            's_ch_avg_DPX','s_ch_s_avg_DPX','max_DPX','min_DPX','average_CX','s_avg_CX','s_ch_avg_CX',
            's_ch_s_avg_CX','max_CX','min_CX','Hydrophobicity']
psaia_col_keep = ['res_name','total_RASA']

## columns
idr_length = 'filter_mobidb'
idr_length_col = 'idr_length_'+str(idr_length)
idr_percentage_col = 'idr_percentage_'+str(idr_length)

group_percentage_col = ['Xle','Aliphaticity','Hydrophobicity','alpha_helix','Charged']
group_percentage_col.extend(['fraction_'+res for res in RESIDUES])
group_percentage_col.extend(['IEP','molecular_weight','gravy'])
group_percentage_col.extend(['Asx','Glx','Pos_charge','Neg_charge','Aromatic','Small','Hydrophilic','beta_turn','beta_sheet','aromaticity'])

columns = ['uniprot_id','ssup_regions','min_sticker_num_percent','sticker_num_percent','pos_sticker_regions','neg_sticker_regions']

for struc in ['idr','ssup']:
    for col in group_percentage_col:
        columns.append('group_'+col+'_'+struc)

for struc in ['idr','ssup']:
    for col in ['hydropathy','Polarity']:
        columns.extend([col + '_'+struc])

# ------------------------------------
# functions
# ------------------------------------
def get_KD_original():
    '''Function which returns the original KD hydropathy lookup table'''
    return  {'ILE': 4.5,                
             'VAL': 4.2,
             'LEU': 3.8,
             'PHE': 2.8,
             'CYS': 2.5,
             'MET': 1.9,
             'ALA': 1.8,
             'GLY': -0.4,
             'THR': -0.7,
             'SER': -0.8,
             'TRP': -0.9,
             'TYR': -1.3,
             'PRO': -1.6,
             'HIS': -3.2,
             'GLU': -3.5,
             'GLN': -3.5,
             'ASP': -3.5,
             'ASN': -3.5,
             'LYS': -3.9,
             'ARG': -4.5}

def get_KD_shifted():
    """ 
    Function which returns the shifted KD hydropathy lookup table (such that
    it runs from 0 to 9 instead of -4.5 to 4.5)
    """
    
    original  = get_KD_original()
    
    shifted = {}
    for i in original:
        shifted[i] = original[i]+4.5

    return shifted

def get_KD_uversky():
    """
    Returns a 0-to-1 normalized KD scale

    """
    shifted = get_KD_shifted()
    
    uversky = {}
    for i in shifted:
        uversky[i] = shifted[i]/9.0

    return uversky

def intervals_extract(iterable):
    '''Convert list of sequential number into intervals'''
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]
        

def intervals_format(intervals):
    if len(intervals)==0:
        formated_intervals = ''
    else:
        formated_intervals = [str(i)+'-'+str(j) if i!=j else str(i) for i,j in intervals]
        formated_intervals = ','.join(formated_intervals)
    return formated_intervals

def cal_group_percentage(seq):
    df = dict()
    for res in RESIDUES:
        res_count =seq.count(res)
        seq_len = len(seq)
        df['fraction_'+res] = res_count/seq_len
    
    seqanalysis = ProteinAnalysis(seq)
    IEP = seqanalysis.isoelectric_point()
    molecular_weight = seqanalysis.molecular_weight()
    gravy = seqanalysis.gravy()

    Xle =df['fraction_I'] + df['fraction_L']
    Aliphaticity=df['fraction_V'] + df['fraction_I'] + df['fraction_L']+ df['fraction_M']
    Hydrophobicity = (df['fraction_V'] + df['fraction_I'] + df['fraction_L']
                                            + df['fraction_F'] + df['fraction_W'] + df['fraction_Y'] 
                                            + df['fraction_M'])
    alpha_helix=(df['fraction_V'] + df['fraction_I'] + df['fraction_Y'] + df['fraction_F']
                      + df['fraction_W'] + df['fraction_L'])
    Charged = df['fraction_K'] + df['fraction_R'] + df['fraction_D'] + df['fraction_E']

    Asx=df['fraction_D'] + df['fraction_N']
    Glx=df['fraction_E'] + df['fraction_Q']
    Pos_charge=df['fraction_K'] + df['fraction_R'] + df['fraction_H']
    Neg_charge=df['fraction_D'] + df['fraction_E']
    Aromatic=df['fraction_F'] + df['fraction_W'] + df['fraction_Y'] + df['fraction_H']
    Small=df['fraction_P'] + df['fraction_G'] + df['fraction_A'] + df['fraction_S']
    Hydrophilic=(df['fraction_S'] + df['fraction_T'] + df['fraction_H'] + 
                                df['fraction_N'] + df['fraction_Q'] + df['fraction_E'] +
                                df['fraction_D'] + df['fraction_K'] + df['fraction_R'])
    beta_turn=df['fraction_N'] + df['fraction_P'] + df['fraction_G'] + df['fraction_S']
    beta_sheet=df['fraction_E'] + df['fraction_M'] + df['fraction_A'] + df['fraction_L']
    # Calculates the aromaticity value of a protein according to Lobry, 1994. 
    # It is simply the relative frequency of Phe+Trp+Tyr.
    aromaticity=df['fraction_F'] + df['fraction_W'] + df['fraction_Y']
    
    out = [Xle,Aliphaticity,Hydrophobicity,alpha_helix,Charged]
    out.extend([df['fraction_'+res] for res in RESIDUES])
    out.extend([IEP,molecular_weight,gravy])
    out.extend([Asx,Glx,Pos_charge,Neg_charge,Aromatic,Small,Hydrophilic,beta_turn,beta_sheet,aromaticity])
    
    return out


def calculate_features(tmpDir, cutoff, files, uids, pdbInfo, ignore,dsspFile):
    normalizedKD  = get_KD_uversky()

    ### load data
    if pdbInfo.empty:
        raise_error('Please specify pdb information data.')
    
    psaiaDir = tmpDir + '/PSAIA_result/'
    if not os.path.isdir(psaiaDir):
        if os.path.isfile(dsspFile):
            dssp = pd.read_pickle(dsspFile)
        else:
            raise_error(f'The dssp file is empty: {dsspFile}')
        
    ### calculate features
    final_df = list()
    for uid, file in zip(uids, files):
        info('UniProt ID: ' + uid)
        info('PDB file: ' + file)
        pdbName = os.path.basename(file).replace('.pdb','')
        
        if os.path.isdir(psaiaDir):
            psaia_file = run_cmd('ls '+psaiaDir+pdbName+'_20*_unbound.tbl', 'result')
            if not psaia_file:
                raise_error(f'{file} has no PSAIA result.')
            info('PSAIA result file: ' + psaia_file)
            
            data = pd.read_csv(psaia_file,skiprows=8,sep='\s+',header=None)
            data.columns=psaia_col
            data = data[psaia_col_keep]
            
            data['rsa'] = (data['total_RASA']/100).values
            data['aa_short'] = data['res_name'].map(residue_key2).tolist()
        else:
            data = pd.DataFrame()
            data['aa_short'] = list(dssp.loc[uid,'aas'].split(','))
            data['res_name'] = data['aa_short'].map(residue_key).tolist()
            data['standard_asa']=data['aa_short'].map(dssp_standard_asa).values
            data['asa'] = list(map(float,dssp.loc[uid,'dssp_asa_score'].split(',')))
            data['rsa'] = np.round(data['asa']/data['standard_asa'],3)
        
        data['Polarity'] = data['aa_short'].map(Polarity).values
        data['HP'] = data['res_name'].map(normalizedKD).values
        data['pos_charge'] = np.where(data['aa_short'].isin(['K','R']),1,0)
        data['neg_charge'] = np.where(data['aa_short'].isin(['D','E']),1,0)
        data['alipha'] = np.where(data['aa_short'].isin(['V','L','I','M']),1,0)
        
        idr_binary_name = 'idr_binary_'+idr_length
        idr_binary = list(map(int,list(pdbInfo.loc[uid, idr_binary_name])))
        data[idr_binary_name] = idr_binary
        data['non_'+idr_binary_name] = np.logical_not(data[idr_binary_name]).astype('int')
        
        i = cutoff
        superficial = '_'.join(['superficial',str(i)])
        data[superficial] = np.where(data['rsa']>i/100,1,0)
        idr = 'idr_binary_'+idr_length
        non_idr = 'non_idr_binary_'+idr_length
        ssup = '_'.join(['ssup',str(i),idr_length])
        data[ssup] = data[superficial] & data[non_idr]
        
        tmp = [uid]
        
        ## sticker features
        binary = data[ssup].tolist()
        numbers_ssup = [kk+1 for kk in range(len(binary)) if binary[kk]==1]

        ssup_length = len(numbers_ssup)
        if ssup_length<=1:
            tmp.extend(['',0,0,'',''])
        else:
            out = list(intervals_extract(numbers_ssup))
            cluster_regions = intervals_format(out)
            tmp.extend([cluster_regions])
            binary = (data['pos_charge'] & data[ssup]).tolist()
            numbers = [kk+1 for kk in range(len(binary)) if binary[kk]==1]
            binary = (data['neg_charge'] & data[ssup]).tolist()
            numbers_neg = [kk+1 for kk in range(len(binary)) if binary[kk]==1]
            min_sticker_num_percent,sticker_num_percent,pos_cluster_regions,neg_cluster_regions = cal_sticker_feature(numbers,numbers_neg,numbers_ssup,file,ssup_length)
            tmp.extend([min_sticker_num_percent,sticker_num_percent,pos_cluster_regions,neg_cluster_regions])
            
        # group percentage
        for struc in [idr,ssup]:
            if ignore and struc.startswith('idr'):
                result = [0]*len(group_percentage_col)
            else:
                seq = ''.join(list(data['res_name'][data[struc]==1].map(residue_key2)))
                if len(seq)==0:
                    result = [0]*len(group_percentage_col)
                else:
                    result = cal_group_percentage(seq)
            tmp.extend(result)
            
        # type2
        for struc in [idr,ssup]:
            for rcol in ['HP','Polarity']:
                if ignore and struc.startswith('idr'):
                    result = 0
                else:
                    if struc == 'PRT':
                        transData = data[rcol]
                    else:
                        transData = data[rcol][data[struc]==1]
                    if len(transData)==0:
                        result = 0
                    else:
                        result = np.mean(transData)
                tmp.append(round(result,3))
        
        final_df.append(tmp)

    ## merge result
    final_df = pd.DataFrame(final_df, columns=columns)
    if ignore:
        final_df[idr_length_col] = 0
        final_df[idr_percentage_col] = 0
    else:
        final_df = pd.merge(final_df, pdbInfo[[idr_length_col,idr_percentage_col]], 
                        how='inner',left_on='uniprot_id',right_index=True)
    final_df.to_pickle(tmpDir + '/final_features.pkl')

    return final_df

