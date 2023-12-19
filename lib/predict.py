#!/usr/bin/env python
# Date: 2023-07-17
# Author: Shuang Hou
# Contact: houshuang@tongji.edu.cn


import os
import numpy as np
import pandas as pd
import pickle
import bz2
from lib.Utility import *
    
# ------------------------------------
# functions
# ------------------------------------
def decompress_data(featureFile):
    features = bz2.BZ2File(featureFile, 'rb')
    features = pickle.load(features)
    return features


def compress_data(outfile,data):
    with bz2.BZ2File(outfile, 'w') as f:
        pickle.dump(data, f)
        
        
def preprocess_data(data, scaler):
    df = data.copy()
    df = df.fillna(0)
    df = df._get_numeric_data()
    columns = df.columns
    df = scaler.transform(df) 
    df = pd.DataFrame(df, columns=columns)
    return df
    

def predict(dataDir,pdata_ori):
    pd.set_option('mode.chained_assignment', None)
    
    fea1 = ['sticker_num_percent', 'min_sticker_num_percent', 'idr_length_filter_mobidb', 'idr_percentage_filter_mobidb', 
            'group_beta_sheet_idr', 'group_fraction_P_idr', 'group_fraction_Q_idr', 'group_Small_idr', 'group_fraction_G_idr', 
            'group_fraction_W_idr', 'hydropathy_idr', 'group_Asx_idr', 'group_fraction_E_idr', 'group_fraction_H_idr', 
            'group_fraction_V_idr', 'Polarity_idr', 'group_fraction_M_idr', 'group_fraction_T_idr', 'group_alpha_helix_idr', 
            'group_fraction_N_idr', 'group_Glx_idr', 'group_fraction_Y_idr', 'group_Pos_charge_idr', 'group_fraction_F_idr', 
            'group_Neg_charge_idr', 'group_fraction_I_idr', 'group_Xle_idr', 'group_Aliphaticity_idr', 'group_IEP_idr', 
            'group_fraction_S_idr', 'group_aromaticity_idr', 'group_fraction_R_idr', 'group_fraction_C_idr', 'group_fraction_K_idr', 
            'group_Hydrophilic_idr', 'group_Hydrophobicity_idr', 'group_beta_turn_idr', 'group_molecular_weight_idr', 'group_fraction_A_idr', 
            'group_Aromatic_idr', 'group_fraction_D_idr', 'group_fraction_L_idr', 'group_Charged_idr', 
            'group_beta_sheet_ssup', 'group_fraction_P_ssup', 'group_fraction_Q_ssup', 'group_Small_ssup', 'group_fraction_G_ssup', 
            'group_fraction_W_ssup', 'hydropathy_ssup', 'group_Asx_ssup', 'group_fraction_E_ssup', 'group_fraction_H_ssup', 
            'group_fraction_V_ssup', 'Polarity_ssup', 'group_fraction_M_ssup', 'group_fraction_T_ssup', 'group_alpha_helix_ssup', 
            'group_fraction_N_ssup', 'group_Glx_ssup', 'group_fraction_Y_ssup', 'group_Pos_charge_ssup', 'group_fraction_F_ssup', 
            'group_Neg_charge_ssup', 'group_fraction_I_ssup', 'group_Xle_ssup', 'group_Aliphaticity_ssup', 'group_IEP_ssup', 
            'group_fraction_S_ssup', 'group_aromaticity_ssup', 'group_fraction_R_ssup', 'group_fraction_C_ssup', 'group_fraction_K_ssup', 
            'group_Hydrophilic_ssup', 'group_Hydrophobicity_ssup', 'group_beta_turn_ssup', 'group_molecular_weight_ssup', 'group_fraction_A_ssup', 
            'group_Aromatic_ssup', 'group_fraction_D_ssup', 'group_fraction_L_ssup', 'group_Charged_ssup']

    fea2 = fea1 + ['phos_PRT']
    
    type = 'phos' if 'phos_PRT' in pdata_ori.columns else 'noPhos'
    features = fea2 if type=='phos' else fea1
    
    modelFile = dataDir + 'human_phos_scaler.sav' if type=='phos' else dataDir + 'human_noPhos_scaler.sav'
    if not os.path.isfile(modelFile):
        raise_error(f'{modelFile} does not exists. Please download the file from https://github.com/TongjiZhanglab/PSPire and put the file under the data folder of the PSPire installation path.')
    else:
        scaler = decompress_data(modelFile)
    
    score = pdata_ori[['uniprot_id']].copy()
    pdata = pdata_ori[features].copy()
    pdata = preprocess_data(pdata,scaler)
    
    ## predict
    modelFile = dataDir + 'human_phos_model.sav' if type=='phos' else dataDir + 'human_noPhos_model.sav'
    if not os.path.isfile(modelFile):
        raise_error(f'{modelFile} does not exists. Please download the file from https://github.com/TongjiZhanglab/PSPire and put the file under the data folder of the PSPire installation path.')
    else:
        modlist = decompress_data(modelFile)
    
    i = 0
    for model in modlist:
        probability = model.predict_proba(pdata)[:, 1]
        score['probability_'+str(i)] = probability
        i += 1
    
    score = score.set_index('uniprot_id')
    outData = pd.DataFrame(np.mean(score, axis=1), columns=['score'])
    outData.reset_index(inplace=True)
    return outData
