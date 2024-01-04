# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 20:08:02 2023

@author: congt
"""

import pandas as pd
df1 = pd.read_csv('count/count_crossing_miR_filtered_DicerWT-rep1.bed', sep='\t')
df2 = pd.read_csv('count/count_crossing_miR_filtered_DicerWT-rep2.bed', sep='\t')
df3 = pd.read_csv('count/count_crossing_miR_filtered_DicerWT-rep3.bed', sep='\t')
df4 = pd.read_csv('count/count_crossing_miR_filtered_DicerdelHelicase-rep1.bed', sep='\t')
df5 = pd.read_csv('count/count_crossing_miR_filtered_DicerdelHelicase-rep2.bed', sep='\t')
df6 = pd.read_csv('count/count_crossing_miR_filtered_DicerdelHelicase-rep3.bed', sep='\t')

def calculating_5p_percentage(df_input, sample):
    df = df_input.copy()
    df = df[df['miRNA_strand'] == '3p']
    df['Count_5p'] = df.groupby(['miRNA', '5p', 'hairpin'])['Count'].transform('sum')
    df['Count_miRNA'] = df.groupby(['miRNA', 'hairpin'])['Count'].transform('sum')
    df.drop_duplicates(subset=['miRNA', '5p', 'hairpin'], keep='first', inplace=True)
    df = df[df['Count_miRNA'] >= 30] #set minimum read count
    df['5p_percentage'] = df['Count_5p'] / df['Count_miRNA']
    df.drop(['3p', 'Count', '3p_modification', 'Strand', 'Sequence', 'miRNA_strand', 'Count_5p', 'Count_miRNA'], axis=1, inplace=True)
    df.rename(columns={'Count_5p':'Count_5p_'+sample, 'Count_miRNA': 'Count_miRNA_'+sample, '5p_percentage':'5p_percentage_'+sample}, inplace=True)
    '''
    ID: miRNA+hairpin, to select ID that appears in all 3 repeats
    '''
    df['ID'] = df['miRNA'] + '-' + df['hairpin']
    df.reset_index(inplace=True, drop=True)
    return (df)
df1 = calculating_5p_percentage(df1, 'DicerWT_rep1')
df2 = calculating_5p_percentage(df2, 'DicerWT_rep2')
df3 = calculating_5p_percentage(df3, 'DicerWT_rep3')
df4 = calculating_5p_percentage(df4, 'DicerdelHelicase_rep1')
df5 = calculating_5p_percentage(df5, 'DicerdelHelicase_rep2')
df6 = calculating_5p_percentage(df6, 'DicerdelHelicase_rep3')

#%%merge samples
'''
select miRNAs that appear in all 3 dataframes
'''
def merge_repeat(df_input1, df_input2, df_input3, sample):
    df1 = df_input1.copy()
    df2 = df_input2.copy()
    df3 = df_input3.copy()
    
    df1_unique = set(df1['ID'])
    df2_unique = set(df2['ID'])
    df3_unique = set(df3['ID'])
    common_elements = df1_unique & df2_unique & df3_unique
    
    df1 = df1[df1['ID'].isin(common_elements)]
    df2 = df2[df2['ID'].isin(common_elements)]
    df3 = df3[df3['ID'].isin(common_elements)]
    
    df = df1.merge(df2, on=['miRNA', '5p', 'hairpin', 'hairpin_sequence', 'checking_5p', 'ID'], how='outer')
    df = df.merge(df3, on=['miRNA', '5p', 'hairpin', 'hairpin_sequence', 'checking_5p', 'ID'], how='outer')
    
    del df['ID']
    df.reset_index(inplace=True, drop=True)
    df.fillna(0, inplace=True)
    df['Mean_percentage'] = (df['5p_percentage_'+sample+'_rep1'] + df['5p_percentage_'+sample+'_rep2'] + df['5p_percentage_'+sample+'_rep3']) / 3
    
    return (df)

df_wt = merge_repeat(df1, df2, df3, 'DicerWT')
df_delHel = merge_repeat(df4, df5, df6, 'DicerdelHelicase')

#pearson correlation analysis
import numpy as np
data = df_delHel
sample = 'DicerdelHelicase'
data1 = data['5p_percentage_'+sample+'_rep1'].to_numpy()
data2 = data['5p_percentage_'+sample+'_rep2'].to_numpy()
r = np.corrcoef(data1, data2)
print (r)


#%%calculate variation of 5p end
'''
1. identify the most abundant isomir in DicerWT
2. need to check the same isomir in DicerdelHelicase
3. calculate the variation
for each isomir of the miRNA:
absolute distance (nt) x abundancy
sum variation of all isomir is the variation of that miRNA
'''

#1. identify the most abundant isomir in DicerWT
def most_abundant_isomir(df_input):
    df = df_input.copy()
    df.sort_values(['hairpin', 'miRNA', 'Mean_percentage'], ascending=[True, True, False], inplace=True)
    df.drop_duplicates(subset=['hairpin', 'miRNA'], keep='first', inplace=True)
    df = df[['miRNA', 'hairpin', '5p']]
    df.rename(columns={'5p':'Major_5p_DICERWT'}, inplace=True)
    return (df)
df = most_abundant_isomir(df_wt)

#2,3.calculate variation
def calculate_variation(df_input1, df_input2):
    df1 = df_input1.copy()
    df2 = df_input2.copy()
    
    df1 = df1.merge(df2, on=['miRNA', 'hairpin'], how='left')
    
    df1.dropna(inplace=True) #to remove miRNA that do not found in DicerWT sample
    
    df1['Local_variation'] = abs(df1['5p'] - df1['Major_5p_DICERWT']) * df1['Mean_percentage']
    df1['Global_variation'] = df1.groupby(['miRNA', 'hairpin'])['Local_variation'].transform('sum')
    
    df1.sort_values(['hairpin', 'miRNA', 'Mean_percentage'], ascending=[True, True, False], inplace=True)
    #only keep rows of the most abundant isomir. In DicerdelHelicase, it might be different from Major_5p_DICERWT column
    df1.drop_duplicates(subset=['miRNA', 'hairpin'], keep='first', inplace=True)
    
    df1 = df1[['miRNA', '5p', 'Major_5p_DICERWT', 'hairpin', 'hairpin_sequence', 'checking_5p', 'Global_variation']]
    
    df1.reset_index(inplace=True, drop=True)
    
    return (df1)

df_wt = calculate_variation(df_wt, df)
df_delHel = calculate_variation(df_delHel, df)

#%%now classify miRNA into 2-nt and non-2nt cleavage
df = pd.read_csv(path+'hsa_miRNA_sequence_and_structure_miRBAse_v22.1.bed', sep='\t')
df = df[df['miRNA_strand'] == '3p']
def find_5p_position(row):
    return row['hairpin_sequence'].find(row['miRNA_sequence'])
df['5p_position'] = df.apply(find_5p_position, axis=1)

df = df[['miRNA', 'hairpin', 'hairpin_sequence', 'miRNA_sequence', '5p_position', 'dot_struct', 'forgi_struct', 'new_define_struct2']]

df_wt = df_wt.merge(df, on=['miRNA', 'hairpin', 'hairpin_sequence'], how='left')
df_delHel = df_delHel.merge(df, on=['miRNA', 'hairpin', 'hairpin_sequence'], how='left')

'''
1. need to use most abundant isomir in DICER-WT to check if it is 2-nt from loop --> observe the change of that isomir in DICER-WT and DICERdelHelicase
2. most abundant isomir in DICERdelHelicase might be different from major_5p_DICERWT isomir.
'''
#1. checking if 2-nt from loop using DICER-WT sample

def classify_isomir(df_input):
    df = df_input.copy()
    for i,seq in enumerate(df['hairpin_sequence']):
        position = df['checking_5p'][i]
        forgi_struct = df['forgi_struct'][i]
        new_df_struct = df['new_define_struct2'][i]
        forgi_split = forgi_struct.split(' ')
        
        if new_df_struct[position - 1] == 'M' and new_df_struct[position - 2] == 'M':
            forgi_minus_one = int(forgi_split[position - 2 + 1])
            if new_df_struct[position - 3] != 'M' or new_df_struct[forgi_minus_one] != 'M':
                df.loc[i,'2-nt (from DICER-WT isomir)'] = 'Yes'
            else:
                df.loc[i,'2-nt (from DICER-WT isomir)'] = 'No'
        else:
            df.loc[i,'2-nt (from DICER-WT isomir)'] = 'No'    
    
    df = df[['miRNA','hairpin','Major_5p_DICERWT','2-nt (from DICER-WT isomir)']]
    return (df)
df = classify_isomir(df_wt)

def assign_feature(df_input1, df_input2, sample):
    df1 = df_input1.copy()
    df2 = df_input2.copy()
    
    df1 = df1.merge(df2, on=['miRNA','hairpin','Major_5p_DICERWT'], how='left')
    
    df1['ID'] = df1['miRNA'] + '-' + df1['hairpin']
    df1['Sample'] = [sample] * len(df1)
    
    return (df1)

df_wt = assign_feature(df_wt, df, 'DicerWT')
df_delHel = assign_feature(df_delHel, df, 'DicerdelHelicase')
    
dfwt_unique = set(df_wt['ID'])
dfdelHel_unique = set(df_delHel['ID'])

common_elements = dfwt_unique & dfdelHel_unique

df_wt = df_wt[df_wt['ID'].isin(common_elements)]
df_delHel = df_delHel[df_delHel['ID'].isin(common_elements)]
    
#%%plotting

df_merge = pd.concat([df_wt, df_delHel], axis=0, ignore_index=True)

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

ax = plt.figure(figsize=(3,3))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True
my_color = sns.color_palette("Blues", as_cmap=True)

ax = plt.figure(figsize=(3,3))


ax = sns.boxplot(data=df_merge, x='2-nt (from DICER-WT isomir)', y='Global_variation', palette=['limegreen','dodgerblue'], hue='Sample',
                 showfliers=False,linewidth=1,zorder=10,hue_order=['DicerWT','DicerdelHelicase'],order=['Yes','No'])
ax = sns.stripplot(data=df_merge, x='2-nt (from DICER-WT isomir)', y='Global_variation', palette=['limegreen','deepskyblue'], edgecolor=".2", 
                   hue='Sample',dodge=True,zorder=1,size=2,hue_order=['DicerWT','DicerdelHelicase'],order=['Yes','No'])

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)

ax.get_legend().remove()

plt.ylim(-0.05,2.05)
plt.yticks([0,0.5,1,1.5,2])

plt.ylabel('')
plt.xlabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)

plt.savefig('boxplot_global_variation_DICERWT_DICERdelHelicase_rescue.png', dpi=150, bbox_inches='tight')
plt.show()

#%%check p-val and Cohen's d
import pandas as pd
import scipy.stats as stats
from scipy.stats import wilcoxon
import numpy as np
def cohen_d(group1, group2):
    mean_diff = np.mean(group1) - np.mean(group2)
    pooled_std = np.sqrt((np.std(group1, ddof=1) ** 2 + np.std(group2, ddof=1) ** 2) / 2)
    return mean_diff / pooled_std

for order in ['Yes','No']:
    df_check = df_merge[df_merge['2-nt (from DICER-WT isomir)'] == order]
    list_wt = df_check[df_check['Sample'] == 'DicerWT']['Global_variation'].tolist()
    list_delHel = df_check[df_check['Sample'] == 'DicerdelHelicase']['Global_variation'].tolist()
    result = stats.ttest_ind(list_wt, list_delHel)
    print (result)
    print (cohen_d(list_wt, list_delHel))
    
    statistic, p_value = wilcoxon(list_wt, list_delHel)
    
    print ('Wilcoxon test ')
    print (statistic)
    print (p_value)



    
    


































