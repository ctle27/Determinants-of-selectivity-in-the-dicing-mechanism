# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 22:22:19 2022

@author: congt
"""

ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
tri_nu_comb = []
for nu1 in ['A','T','G','C']:
    for nu2 in ['A','T','G','C']:
        for nu3 in ['A','T','G','C']:
            tri_nu_comb.append(nu1+nu2+nu3)
seq = ''
seq_list = []  
for i in range(1,6):
     for comb1 in tri_nu_comb:
         for comb2 in tri_nu_comb:
             seq = ini_seq[:i+18] + comb1 + ini_seq[i+21:67] + ini_seq[67:81-i] + comb2 + ini_seq[84-i:]
             if i == 5:
                 seq = seq[:22] + 'C' + seq[23:79] + 'G' + seq[80:]
             if seq not in seq_list:
                 seq_list.append(seq)
                 
for i in range(1,7):
    for comb1 in tri_nu_comb:
        for comb2 in tri_nu_comb:
            seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:21] + ini_seq[21:81] + ini_seq[81:87-i] + comb2 + ini_seq[90-i:]
            if seq not in seq_list:
                seq_list.append(seq)
                
var_list = []
shRNA_list = []
truncated_shRNA_list = []

import pandas as pd
df = pd.DataFrame()
for i, seq in enumerate(seq_list):
    var_list.append('>Variant_'+str(i+1).zfill(5)) #zfill outputs 00001 instead of 1
    shRNA_list.append(seq[:35] + seq[67:])
    truncated_shRNA_list.append(seq[13:35] + seq[67:89])
df['Variant'] = var_list
df['shRNA_list'] = shRNA_list
df['truncated_shRNA_list'] = truncated_shRNA_list

f = open('DICER-delHelicase-lib-reference.fa','w+')
for i,var in enumerate(var_list):
  f.write(var+'\n')
  f.write(shRNA_list[i]+'\n')
f.close()

df.to_csv('DICER-delHelicase-lib-reference.bed', sep='\t', index=False, header=False)
