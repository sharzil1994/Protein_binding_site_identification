# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 06:36:25 2022

@author: PC
"""

import pickle
import numpy as np

# from tqdm import tqdm
# from Bio import SeqIO
import os
import pandas as pd

import subprocess


def site_predictor(input_ID):
    path='./Feature/'
     
    with open(f'{path}/data_seq/data_ID.dat','rb') as FH:
        data0 =pickle.load(FH)
        
    found=0
    i_0=0
    for index_0 in data0[:,0]:  
        if input_ID.lower()==index_0.lower():
            found=1
            seq=data0[i_0,1]
            break
        i_0=i_0+1
    
    if found==1:
        pass
    else :
        print('ID not present in dataset')
        return
        
        
    
    tEMP_0=np.load(f'{path}pssm/{input_ID}.npy')
    tEMP_1=np.load(f'{path}dssp/{input_ID}.npy')
    tEMP_2=np.load(f'{path}hmm/{input_ID}.npy')
    
    
    
    def windowing(features,seq,w_size,start,stop):   
    
        # all_features=features
        seq_len=len(seq)
        
        
        # a=all_features
        fea_len=np.shape(features)[1]
        # print(fea_len)
    
        finalout1=np.zeros([seq_len,w_size,fea_len],'float')
        # finalout1=np.zeros([to,w_size,69],'float')
        
        l=0
        for j in range( 0, seq_len):
            
            for k in range(0,w_size):
                
                k1=int(j+k-((w_size-1)/2))
                
                if k1<0 or k1 > seq_len-1:
                    pass
                else:
                    finalout1[l,k,:]=features[k1,start:stop]
    
            l=l+1    
        return finalout1
    
    import tensorflow as tf
    
    w_size=19
    f1=windowing(tEMP_0,seq,w_size,0,20)
    f2=windowing(tEMP_1,seq,w_size,0,14)
    f3=windowing(tEMP_2,seq,w_size,0,20)
    new_model = tf.keras.models.load_model('./model/model_19.h5')
    pred_y = new_model.predict([f1,f2,f3])
    
    for i in range(len(pred_y)):
        if pred_y[i] < 0.5:
            pred_y[i] = 0;
        else:
            pred_y[i] = 1;
    
    pred_y=np.squeeze(pred_y)        
    pred_y=pred_y.astype('int')
    
    out=''
    for index_0 in pred_y:
        out=out+str(index_0)        
    
    print (seq)
    print (out)
    
    return
    # return seq , out


# input_ID='3zeuD'
input_ID=input('give protien id and chain : example 3zeuD\n')
site_predictor(input_ID)
    
    
    
    
    
    
    
    