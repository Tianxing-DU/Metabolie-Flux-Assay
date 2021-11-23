"""如果只看XN, Genistein和dmso，两组之间对比(数据需要提前提取一次)"""

import pandas as pd, matplotlib.pyplot as plt, seaborn as sns, numpy as np
from sklearn import preprocessing
from library import *

# 1. input:
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Behandlung_3T3-L1 2011-12-20_2014-01-27_corr.xls'
main_df = pd.read_excel(file_path)

# 2. creat simplified and purified df
main_df.replace(' ', '', regex=True, inplace=True)
main_df.replace('µM', '', regex=True, inplace=True)
main_df.replace(',', '.', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '100-C0:0': '229-C18:0-DC']
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)
Extract_str_df = main_df[['Treatment', 'Charge', 'CSA']]

# for x in dict_LOD:
#     Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN  # remove of values under LOD

Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # normalisation based on CSA

x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)

# Extract_int_df_log = Extract_int_df_div.apply(np.log)  # logarithmising the data for better visualisation
# Extract_int_df_log.replace([np.inf, -np.inf], np.nan, inplace=True)
# log_df = Extract_int_df_log.dropna(axis=1)
# log使用的问题就是之后不能使用heatmap的形式，dropnan过后又会导致我的数据不够全面

str_df = main_df[['Treatment', 'Charge']]
simplified_df = pd.concat([str_df, Extract_int_df_normalised], axis=1)  # without normalisation of CSA

# 3. select positive charged data
simplified_oleic_yes_df = simplified_df.loc[simplified_df['Charge'] == 'Oleic+']
simplified_oleic_no_df = simplified_df.loc[simplified_df['Charge'] == 'Oleic-']
pca_p0_df = simplified_oleic_yes_df.drop(columns='Charge').reset_index(drop=True)

# 4. select wished treatments for comparison in cluster (for single comparison)
concentrated0_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'DMSO']
concentrated1_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Resveratrol']
concentrated2_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Xanthohumol']
concentrated3_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'EGCG']
concentrated4_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Genistein']
concentrated5_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Bezafibrat']
concentrated6_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Etomoxir']
concentrated7_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'K']
con_list = [concentrated1_df, concentrated2_df, concentrated3_df, concentrated4_df, concentrated5_df,
            concentrated6_df]
pca_p0_df = pd.concat(con_list, axis=0).reset_index(drop=True)

# 5. clustermap
cluster_df = pca_p0_df.iloc[:, 1:45]
sns.set_theme(style='ticks', font='Avenir', font_scale=.6)
index = pca_p0_df.loc[:, ['Treatment']]
# 这种方式改index值得学习
cluster_df.index = index
cm = sns.clustermap(cluster_df, yticklabels=True, xticklabels=True, dendrogram_ratio=0.2, figsize=(12, 17),
                    colors_ratio=0.0, cmap='CMRmap_r')
# 在最开始调figuresize使用plt是不起作用的，必须要在sns这里内部调整
# 这里发现，如果数据不使用log那么数据的分离就不好.但我觉得使用自定义那么也是ok的
# 调整颜色要用cmap！！ cbar，palette都不行（可选用颜色：gnuplot2_r,CMRmap_r
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210921_6h/clustermap_mouse_drop_DMSO_K.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
