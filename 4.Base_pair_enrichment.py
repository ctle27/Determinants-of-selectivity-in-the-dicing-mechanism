#%%%   DRES
'''
plotting for DRES
'''

import pandas as pd

df = pd.read_csv('list_of_333_DRES.txt', sep='\t')

df.fillna('N-N', inplace=True)

for i,pos__4 in enumerate(df['Position-4']):
    pos__3 = df['Position-3'][i]
    pos__2 = df['Position-2'][i]
    pos__1 = df['Position-1'][i]
    pos_0 = df['Position0'][i]
    pos_1 = df['Position1'][i]
    
    df.loc[i, '5p-arm'] = pos__4[0] + pos__3[0] + pos__2[0] + pos__1[0] + pos_0[0] + pos_1[0]
    df.loc[i, '3p-arm'] = pos__4[2] + pos__3[2] + pos__2[2] + pos__1[2] + pos_0[2] + pos_1[2]

import logomaker as lm
import matplotlib.pyplot as plt

def draw_weblogo (list_sequence, fig_name, save_fig = 'no'):
    # counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat['correct_index'] = counts_mat.index.map(lambda x: x+1)
    counts_mat = counts_mat.set_index('correct_index')
    
    '''
    column N represents any nt.
    divide the frequency of N at each position by 4 and add it to each of the column of frequency of A, T, G, C
    '''
    counts_mat['A'] += counts_mat['N'] / 4
    counts_mat['C'] += counts_mat['N'] / 4
    counts_mat['G'] += counts_mat['N'] / 4
    counts_mat['U'] += counts_mat['N'] / 4
    
    del counts_mat['N']
    
    crp_logo  = lm.Logo(counts_mat, 
                        figsize = [0.4*len(counts_mat), 1.2], 
                        color_scheme = {'A': '#60d394', 'C': '#f08080',  'G': '#3da4dc', 'U': 'lightgrey', 'T': '#3da4dc'},
                        font_name='Arial Rounded MT Bold', zorder = 3)
    
    for _, spine in crp_logo.ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_color('black')
    
    plt.yticks([0,1,2], fontsize = 0, color = 'white')
    plt.xticks([1,2,3,4,5,6], fontsize = 0, color = 'white')
    plt.tick_params(axis='y', width = 1, length=6)
    plt.tick_params(axis='x', width = 1, length=6)
    plt.grid(axis = 'both', color = 'lightgrey', linestyle = '--', linewidth = 0.5)
    
    if save_fig != 'no':
        # plt.title(save_fig)
        # plt.show()
        plt.savefig(f'{fig_name}.png'.format(save_fig), bbox_inches="tight", dpi =1000)
    else:
        plt.show()
    return()

draw_weblogo(df['5p-arm'].tolist(),  'DRES_5p_arm', save_fig = 'yes')
draw_weblogo(df['3p-arm'].tolist(),  'DRES_3p_arm',  save_fig = 'yes')

#%%%   YCR
'''
plotting for YCR
'''
import pandas as pd

df = pd.read_csv('YCR_scoring.txt', sep='\t')
df = df.head(35) #top YCR motifs


import logomaker as lm
import matplotlib.pyplot as plt

def draw_weblogo (list_sequence, fig_name, save_fig = 'no'):
    # counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat['correct_index'] = counts_mat.index.map(lambda x: x+1)
    counts_mat = counts_mat.set_index('correct_index')
    '''
    remove position -3 (Y) from analysis
    '''
    counts_mat.drop([1], axis=0, inplace=True)
    crp_logo  = lm.Logo(counts_mat, 
                        figsize = [0.4*len(counts_mat), 1.2], 
                        color_scheme = {'A': '#60d394', 'C': '#f08080',  'G': '#3da4dc', 'U': 'lightgrey', 'T': '#3da4dc'},
                        font_name='Arial Rounded MT Bold', zorder = 3)
    
    for _, spine in crp_logo.ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_color('black')
    
    plt.yticks([0,1,2], fontsize = 0, color = 'white')
    plt.xticks([2,3], fontsize = 0, color = 'white')
    plt.tick_params(axis='y', width = 1, length=6)
    plt.tick_params(axis='x', width = 1, length=6)
    plt.grid(axis = 'both', color = 'lightgrey', linestyle = '--', linewidth = 0.5)
    
    if save_fig != 'no':
        # plt.title(save_fig)
        # plt.show()
        plt.savefig(f'{fig_name}.png'.format(save_fig), bbox_inches="tight", dpi =1000)
    else:
        plt.show()
    print (counts_mat)
    return()

draw_weblogo(df["5'-arm (from 5' to 3')"].tolist(),  'YCR_5p_arm', save_fig = 'yes')
draw_weblogo(df["3'-arm (from 3' to 5')"].tolist(),  'YCR_3p_arm',  save_fig = 'yes')





























