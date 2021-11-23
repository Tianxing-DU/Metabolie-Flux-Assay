# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12

# Import
import pandas as pd, matplotlib.pyplot as plt, seaborn as sns, numpy as np, pingouin as pg
from sklearn import preprocessing
from library import *

# Set up infos
# +++ sorting of chain orders

# +++ extra plotting orders
order = ['Kontrolle', 'DMSO', 'Resveratrol', 'Genistein', 'Xanthohumol', 'EGCG', 'Bezafibrat', 'Etomoxir']

# 0. DataFrame Processing:
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Datentabelle_humane Adipocyten 20110221_2014-01-27 corri.xls'
sheet_name = 'Behandlung'
human_df = pd.read_excel(file_path, sheet_name=sheet_name)  # open the correct sheet
human_df.columns = human_df.iloc[0]  # rename the columns
human_df = human_df.drop(index=0).reset_index(drop=True)  # delete the duplicate column names
yaxislabel = 'Acylcarnitine products \n (normalized data on a log scale)'
xaxislabel = 'Oleic acid (µM)'
# +++ reorganise data
human_df['Date'] = human_df['Sample ID'].map(lambda x: x.split('_')[0])
human_df['Cell'] = human_df['Sample ID'].map(lambda x: x.split('_')[1])
human_df['Unknown'] = human_df['Sample ID'].map(lambda x: x.split('_')[2])
human_df['Time'] = human_df['Sample ID'].map(lambda x: x.split('_')[3])
human_df['Repeat'] = human_df['Sample ID'].map(lambda x: x.split('_')[4])
human_df['Charge'] = human_df['Sample ID'].map(lambda x: x.split('_')[5])
human_df['Treatment'] = human_df['Sample ID'].map(lambda x: x.split('_')[6])

# 1. creat simplified and purified df for plotting
main_df = human_df
main_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
main_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
main_df.replace(',', '.', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '100-C0:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = main_df[['Treatment', 'Charge']]
# extract the real interested columns for plotting
# for x in dict_LOD:
#     Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN
# 问题在于计算p-value的时候不能使用这个
# 3. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # for preparation of normalisation by CSA
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back

# 4. merge the int and str df
simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1)

# 5. set new colums with treatments-charge from df
melt_df = pd.melt(simplified_df, id_vars=['Treatment', 'Charge'], var_name='C-Chains',
                  value_name='Normalised Value (mean)')
pair_df = pd.DataFrame(melt_df['Treatment'].str.cat(melt_df['Charge'], sep='-')).rename(
    columns={'Treatment': 'Pair'})  # pair 'Treatment' and 'Charge
final_df = pd.concat([melt_df, pair_df], axis=1).loc[:, ['C-Chains', 'Normalised Value (mean)', 'Pair']]
# all pairs / infos for plotting are included

# 6. calculation of all P-Value
final_df['Subject'] = np.array(range(0, 9504))
# 这里我一般会然程序跑一次，然后报错，就知道这里的range是多少了
pairwise_test_df = pg.pairwise_ttests(dv='Normalised Value (mean)', between=['Pair'], within=['C-Chains'],
                                      padjust='fdr_bh', data=final_df, subject='Subject', alpha=0.05)
# 这里是用的是bh adjusted p-value
a = pairwise_test_df.loc[pairwise_test_df['C-Chains'] != '-'].reset_index(drop=True)  # simplify rows
b = pd.DataFrame(a['A'].str.cat(a['B'], sep=' vs. ')).rename(columns={'A': 'pair'})  # pair 'Treatment'
# 这里str的用法，以及把两个columns合在一起，中间用vs隔开，值得学习

c = a.loc[:, ['C-Chains', 'p-corr', 'p-unc']]  # extract only p-value and BH-correction adjusted p-value
d = pd.concat([b, c], axis=1)

# +++ pval belongs to Xanthohumol vs. DMSO
pval_Xanthohumol_oleic_yes_df = d.loc[d['pair'] == 'DMSO-Oleic+ vs. Xanthohumol-Oleic+']
pval_EGCG_oleic_yes_df = d.loc[d['pair'] == 'EGCG-Oleic+ vs. Kontrolle-Oleic+']
pval_Genistein_oleic_yes_df = d.loc[d['pair'] == 'DMSO-Oleic+ vs. Genistein-Oleic+']
pval_Resveratrol_oleic_yes_df = d.loc[d['pair'] == 'DMSO-Oleic+ vs. Resveratrol-Oleic+']
pval_Kontrolle_oleic_yes_df = d.loc[d['pair'] == 'DMSO-Oleic+ vs. Kontrolle-Oleic+']


# +++ corr Xanthohumol vs. DMSO
"""Xanthohumol_DMSO_with_oleic"""
# !!! 这里曾解决了一个最重要的问题， min()empty, 方法是去掉reverse 函数，而采用[::-1]，即 list_custom[::-1]
# !!! 这里曾尝试过使用subplot但是没有成功，主要是特别复杂的facid 不能（或者我不知道）装入subplot当中
# --- 这里记录一下 subplot我的尝试痕迹
# --- +++ 构建图框 三排，每排一个图，这里并不是把下面三行代码都一次性写下来，而是在每一个子图前写
# plt.subplot(311, sharey=pval_no_ax)
# plt.subplot(312, sharey=pval_no_ax)
# plt.subplot(313, sharey=pval_no_ax)

# --- +++ 画图（这里没有具体的过程，但是和上面普通作图没有区别）

# --- +++ adjust the plot 可以使用的一些指令
# plt.gcf().set_size_inches(5.5,6) #这个不记得有没有用了
# plt.subplots_adjust(left=0.125,
#                     bottom=0.1,
#                     right=0.9,
#                     top=0.9,
#                     wspace=0.2,
#                     hspace=0.35)
# # 一些其他可以调节的部分constrained_layout=True
# plt.tight_layout()  # Important! so that the plot can have a better distribution

plt.figure(figsize=(18.5,5))
Xanthohumol_pval_yes_ax = sns.pointplot(data=pval_Xanthohumol_oleic_yes_df, x='C-Chains', y='p-unc',
                            order=list_custom_r, palette='plasma_r', scale=.5, alpha=.6)
plt.xticks(rotation=90)
Xanthohumol_pval_yes_ax.set(yscale="log")
plt.axhline(0.05, color='red', alpha=.8, linewidth=.8, label='Adjusted p-value= 0.05')
Xanthohumol_pval_yes_ax.legend(frameon=False, fontsize='large', title_fontsize='large')
plt.title('Xanthohumol', fontsize='large', fontweight='bold')
plt.legend(title='p-value',loc=0)
plt.tight_layout()
save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210920_7h/XN_pval_unc.pdf'
plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()


"""EGCG_Kontrolle_with_oleic"""
plt.figure(figsize=(18.5,5))
EGCG_pval_yes_ax = sns.pointplot(data=pval_EGCG_oleic_yes_df, x='C-Chains', y='p-unc',
                            order=list_custom_r, palette='plasma_r', scale=.5, alpha=.6)
plt.xticks(rotation=90)
EGCG_pval_yes_ax.set(yscale="log")
plt.axhline(0.05, color='red', alpha=.8, linewidth=.8, label='Adjusted p-value= 0.05')
EGCG_pval_yes_ax.legend(frameon=False, fontsize='large', title_fontsize='large')
plt.title('EGCG', fontsize='large', fontweight='bold')
plt.legend(title='p-value',loc=0)
plt.tight_layout()
save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210920_7h/EGCG_pval_unc.pdf'
plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

"""Resveratrol_DMSO_with_oleic"""
plt.figure(figsize=(18.5,5))
Resveratrol_pval_yes_ax = sns.pointplot(data=pval_Resveratrol_oleic_yes_df, x='C-Chains', y='p-unc',
                            order=list_custom_r, palette='plasma_r', scale=.5, alpha=.6)
plt.xticks(rotation=90)
Resveratrol_pval_yes_ax.set(yscale="log")
plt.axhline(0.05, color='red', alpha=.8, linewidth=.8, label='Adjusted p-value= 0.05')
Resveratrol_pval_yes_ax.legend(frameon=False, fontsize='large', title_fontsize='large')
plt.title('Resveratol', fontsize='large', fontweight='bold')
plt.legend(title='p-value',loc=0)
plt.tight_layout()
save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210920_7h/Resveratrol_pval_unc.pdf'
plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

"""Genistein_DMSO_with_oleic"""
plt.figure(figsize=(18.5,5))
Genistein_pval_yes_ax = sns.pointplot(data=pval_Genistein_oleic_yes_df, x='C-Chains', y='p-unc',
                            order=list_custom_r, palette='plasma_r', scale=.5, alpha=.6)
plt.xticks(rotation=90)
Genistein_pval_yes_ax.set(yscale="log")
plt.axhline(0.05, color='red', alpha=.8, linewidth=.8, label='Adjusted p-value= 0.05')
Genistein_pval_yes_ax.legend(frameon=False, fontsize='large', title_fontsize='large')
plt.title('Genistein', fontsize='large', fontweight='bold')
plt.legend(title='p-value',loc=0)
plt.tight_layout()
save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210920_7h/Genistein_pval_unc.pdf'
plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

"""Kontrolle_DMSO_with_oleic"""
plt.figure(figsize=(18.5,5))
Kontrolle_pval_yes_ax = sns.pointplot(data=pval_Kontrolle_oleic_yes_df, x='C-Chains', y='p-unc',
                            order=list_custom_r, palette='plasma_r', scale=.5, alpha=.6)
plt.xticks(rotation=90)
Kontrolle_pval_yes_ax.set(yscale="log")
plt.axhline(0.05, color='red', alpha=.8, linewidth=.8, label='Adjusted p-value= 0.05')
Kontrolle_pval_yes_ax.legend(frameon=False, fontsize='large', title_fontsize='large')
plt.title('Kontrolle', fontsize='large', fontweight='bold')
plt.legend(title='p-value',loc=0)
plt.tight_layout()
save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210920_7h/Kontrolle_pval_unc.pdf'
plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
