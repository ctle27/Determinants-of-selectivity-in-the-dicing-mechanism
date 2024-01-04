
#%% control sample
path = 'D:/HKUST_Research/Dicer-loop-counting-redefinition/ngs_analysis/'
import pandas as pd
df_rep1 = pd.read_csv(path+'rep1/Dicer-del-Helicase-control-raw-count.bed',sep='\t',names=['Variant','shRNA_sequence','Raw_count'])
df_rep2 = pd.read_csv(path+'rep2/Control-raw-count-Rep2.bed',sep='\t',names=['Variant','Raw_count'])


df_rep1['RPM_control_rep1'] = df_rep1['Raw_count'] / df_rep1['Raw_count'].sum() * 1000000
df_rep2['RPM_control_rep2'] = df_rep2['Raw_count'] / df_rep2['Raw_count'].sum() * 1000000
del df_rep1['Raw_count']
del df_rep2['Raw_count']

df_ctrl = df_rep1.merge(df_rep2, on=['Variant'], how='inner')
df_ctrl = df_ctrl.merge(df, on=['shRNA_sequence'], how='inner')

df_check = df_ctrl.copy()
df_check['count'] = df_check.groupby(['concrete_struct'])['shRNA_sequence'].transform('count')
df_check.drop_duplicates(subset=['concrete_struct'], keep='first', inplace=True)

#%%cleavage samples
path = 'D:/HKUST_Research/Dicer-loop-counting-redefinition/ngs_analysis/'
df_dc_rep1 = pd.read_csv(path+'rep1/3.Dicerdelhelicase-DC-score.bed',sep='\t')
df_sc_rep1 = pd.read_csv(path+'rep1/3.Dicerdelhelicase-SC-score.bed',sep='\t')

df_dc_rep2 = pd.read_csv(path+'rep2/3.Dicerdelhelicase-DC-score.bed',sep='\t')
df_sc_rep2 = pd.read_csv(path+'rep2/3.Dicerdelhelicase-SC-score.bed',sep='\t')

def pre_process(df_input1,type_clv1,df_input2,type_clv2,rep):
    df1 = df_input1.copy()
    df1['Type-general'] = [type_clv1] * len(df1.index) #add one column to show this is SC or DC cleavage product
    df1['Sum_clv_site_count_'+rep] = df1.groupby(['Variant','Start','End'])['Count'].transform('sum')
    df1.sort_values(['Variant','Start','End','Count'],ascending=[True,True,True,False],inplace=True)
    df1.drop_duplicates(subset=['Variant','Start','End'],keep='first',inplace=True)
    df1.drop(['Sequence','Relative_score','Count','RPM'],axis=1,inplace=True)
    
    df2 = df_input2.copy()
    df2['Type-general'] = [type_clv2] * len(df2.index) #add one column to show this is SC or DC cleavage product
    df2['Sum_clv_site_count_'+rep] = df2.groupby(['Variant','Start','End'])['Count'].transform('sum')
    df2.sort_values(['Variant','Start','End','Count'],ascending=[True,True,True,False],inplace=True)
    df2.drop_duplicates(subset=['Variant','Start','End'],keep='first',inplace=True)
    df2.drop(['Sequence','Relative_score','Count','RPM'],axis=1,inplace=True)
    
    #merge df_dc and df_sc vertically and merge with df_control to get information of new define structure and structure id in df_ctrl.
    df = pd.concat([df1,df2], ignore_index=True)
    return df

df_delhel_rep1 = pre_process(df_dc_rep1, 'DC',df_sc_rep1, 'SC','rep1')
df_delhel_rep2 = pre_process(df_dc_rep2, 'DC',df_sc_rep2, 'SC','rep2')

df_delhel = df_delhel_rep1.merge(df_delhel_rep2, on=['Variant','shRNA_sequence','Start','End','Cleavage-type','Cleavage-site','Type-general'], how='outer')
df_delhel.fillna(0.1, inplace=True)

df_delhel = df_ctrl.merge(df_delhel, on=['Variant', 'shRNA_sequence'], how='outer')

#%%since some structures contain bulges --> need to reassign cleavage sites
def re_assign_clv_site(df_input):
    df = df_input.copy()
    df.dropna(inplace=True)
    df.reset_index(inplace=True, drop=True)
    for i,strc in enumerate(df['new_define_struct2']): #struct2 converts AAA-AA in struct1 to ASS-SS
        x = int(df['Start'][i])
        y = int(df['End'][i]) + 32 #32N barcode
        if x < 6:
            bp_5p = str(0) #all 3p SC will be assigned to 0-y, doesnt matter 5p starts from 0,1,2,3,4 or 5
        elif x >= 6:
            if strc[x-1] == 'M' or strc[x-1] == 'S':
                bp_5p = str(strc[:x].count('M') + strc[:x].count('S') + strc[:x].count('F'))
            elif strc[x-1] != 'M' and strc[x-1] != 'S':
                segment_5p = strc[:x].replace('S','M') #replace S with M for convenient of finding the last M/S
                pos_5p = segment_5p.rfind('M')
                bp_5p = str(strc[:x].count('M') + strc[:x].count('S') + strc[:x].count('F')) + strc[pos_5p+1:x-1]
        if y in range(100,105):
            bp_3p = str(0) #all 5p SC will be assigned to x-0, doesnt matter 3p ends at 68,69,70,71 or 72
        if y not in range(100,105):
            if strc[y] == 'M' or strc[y] == 'S':
                bp_3p = str(strc[y:].count('M') + strc[y:].count('S') + strc[y:].count('T'))
            elif strc[y] != 'M' and strc[y] != 'S':
                segment_3p = strc[y:102].replace('S','M') #replace S with M for convenient of finding the first M/S
                pos_3p = segment_3p.find('M')
                bp_3p = str(strc[y:].count('M') + strc[y:].count('S') + + strc[y:].count('T')) + strc[y:y+pos_3p]
        '''
        #assign annotation for 1 of 15 clv site: DC-2nt and SC from 19-23.
        other cases: DC-1nt, DC-3nt, cut at B and A to 'other. results stored in 5p-3p-alternative
        '''
        combine = bp_5p + '-' +  bp_3p
        df.loc[i,'5p-3p'] = combine
        if bp_5p == bp_3p:
            df.loc[i,'5p-3p-alternative'] = combine
        if bp_5p != bp_3p:
            if 'A' not in combine and 'B' not in combine:
                if bp_5p == '0' or bp_3p == '0':
                    df.loc[i,'5p-3p-alternative'] = combine
                if bp_5p != '0' and bp_3p != '0':
                    df.loc[i,'5p-3p-alternative'] = 'other'
        if bp_5p != bp_3p:
            if 'A' in combine or 'B' in combine:
                df.loc[i,'5p-3p-alternative'] = 'other'
        if bp_5p == '0':
            df.loc[i,'SC_on'] = '3p'
        elif bp_3p == '0':
            df.loc[i,'SC_on'] = '5p'
        else:
            df.loc[i,'SC_on'] = 'None'
                
    # remove SC and DC at position 18. because in initial filter, only select
    '''
    in initial filter, only select cleavage site from 19-23
    now, due to bulge, need to re-assign cleavage site. some cleavage site may change from 19 or 20 to 18
    it does not mean that other variants were not cleaved at position 18, however, the reads were filtered out
    therefore, need to remove clv site at position 18 in this step
    '''
    df = df[~df['5p-3p'].isin(['18-18','0-18','18-0'])]
    '''
    #re-calculate RPM
    '''
    df['Sum_clv_type_rep1'] = df.groupby(['Type-general'])['Sum_clv_site_count_rep1'].transform('sum')
    df['Sum_clv_type_rep2'] = df.groupby(['Type-general'])['Sum_clv_site_count_rep2'].transform('sum')
    
    df['RPM_CP_rep1'] = df['Sum_clv_site_count_rep1'] / df['Sum_clv_type_rep1'] * 1000000
    df['RPM_CP_rep2'] = df['Sum_clv_site_count_rep2'] / df['Sum_clv_type_rep2'] * 1000000

    df.drop(['Sum_clv_site_count_rep1','Sum_clv_site_count_rep2','Cleavage-type','Cleavage-site'
             ,'Sum_clv_type_rep1','Sum_clv_type_rep2'], axis=1, inplace=True)
    return df

reassigned_df_delhel = re_assign_clv_site(df_delhel)


#%%now calculate cleavage efficiency, accuracy
def calculating_accuracy_efficiency(df_input):
    df = df_input.copy()
    import numpy as np
    for rep in ['_rep1','_rep2']:
        df['Sum_clv_site_RPM'+rep] = df.groupby(['Variant','5p-3p','Type-general'])['RPM_CP'+rep].transform('sum') #--> give positional clv efficiency
        df['Sum_alternative_clv_site_RPM'+rep] = df.groupby(['Variant','5p-3p-alternative','Type-general'])['RPM_CP'+rep].transform('sum') # merge all 'other' using this command
        df['Sum_variant_RPM_in_each_type'+rep] = df.groupby(['Variant','Type-general'])['RPM_CP'+rep].transform('sum')
        
        df.drop_duplicates(subset=['Variant','5p-3p','Type-general'], keep='first', inplace=True) 
        df.sort_values(['Variant'], ascending=True, inplace=True)
        df.reset_index(inplace=True)
        df.drop(['index','RPM_CP'+rep], axis=1, inplace=True)
     
        df['Positional_efficiency'+rep] = np.log2(df['Sum_clv_site_RPM'+rep]+0.1) - np.log2(df['RPM_control'+rep]+0.1)
        df['Positional_efficiency_of_alternative_5p_3p'+rep] = np.log2(df['Sum_alternative_clv_site_RPM'+rep]+0.1) - np.log2(df['RPM_control'+rep]+0.1)
        
        df['Type_general_efficiency'+rep] = np.log2(df['Sum_variant_RPM_in_each_type'+rep]+0.1) - np.log2(df['RPM_control'+rep]+0.1)
       
        df['Cleavage_accuracy'+rep] = df['Sum_clv_site_RPM'+rep] / df['Sum_variant_RPM_in_each_type'+rep]
        df['Cleavage_accuracy_of_alternative_5p_3p'+rep] = df['Sum_alternative_clv_site_RPM'+rep] / df['Sum_variant_RPM_in_each_type'+rep]
    
    df['Mean_Positional_efficiency'] = (df['Positional_efficiency_rep1']+df['Positional_efficiency_rep2']) / 2
    df['Mean_Positional_efficiency_of_alternative_5p_3p'] = (df['Positional_efficiency_of_alternative_5p_3p_rep1']+df['Positional_efficiency_of_alternative_5p_3p_rep2']) / 2
    df['Mean_Type_general_efficiency'] = (df['Type_general_efficiency_rep1']+df['Type_general_efficiency_rep2']) / 2
    df['Mean_Cleavage_accuracy'] = (df['Cleavage_accuracy_rep1']+df['Cleavage_accuracy_rep2']) / 2
    df['Mean_Cleavage_accuracy'] = (df['Cleavage_accuracy_rep1']+df['Cleavage_accuracy_rep2']) / 2
    df['Mean_Positional_efficiency'] = (df['Positional_efficiency_rep1']+df['Positional_efficiency_rep2']) / 2
    df['Mean_Cleavage_accuracy_of_alternative_5p_3p'] = (df['Cleavage_accuracy_of_alternative_5p_3p_rep1']+
                                                         df['Cleavage_accuracy_of_alternative_5p_3p_rep2']) / 2
    df.drop(['Start','End','Sum_alternative_clv_site_RPM_rep1','Sum_alternative_clv_site_RPM_rep2'], axis=1, inplace=True)  
    return df

processed_df_delhel = calculating_accuracy_efficiency(reassigned_df_delhel)
