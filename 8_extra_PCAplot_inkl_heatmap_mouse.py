import pandas as pd, matplotlib.pyplot as plt, seaborn as sns, numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
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

for x in dict_LOD:
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN
# remove of values under LOD
# 无论这里产不产生NaN值，后面的log都会出现NaN值，而减少我的样品数量
# 然而使用普通的标准化，又不能让数据分离

Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # normalisation based on CSA

# 用了log后，数据已经被很好的分离，就没有在使用这里的min max normalisation
# x = Extract_int_df_div.values  # returns a numpy array
# min_max_scaler = preprocessing.MinMaxScaler()
# x_scaled = min_max_scaler.fit_transform(x)
# Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns).dropna(axis=1)

Extract_int_df_log = Extract_int_df_div.apply(np.log)  # logarithmising the data for better visualisation
Extract_int_df_log.replace([np.inf, -np.inf], np.nan, inplace=True)
log_df = Extract_int_df_log.dropna(axis=1)

str_df = main_df[['Treatment', 'Charge']]
simplified_df = pd.concat([str_df, log_df], axis=1)  # without normalisation of CSA

# 3. select positive charged data
simplified_oleic_yes_df = simplified_df.loc[simplified_df['Charge'] == 'Oleic+']
simplified_oleic_no_df = simplified_df.loc[simplified_df['Charge'] == 'Oleic-']
pca_p0_df = simplified_oleic_yes_df.drop(columns='Charge').reset_index(drop=True)

# 4. select wished treatments for comparison in PCA
concentrated0_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'DMSO']
concentrated1_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Resveratrol']
concentrated2_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Xanthohumol']
concentrated3_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'EGCG']
concentrated4_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Genistein']
concentrated5_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Bezafibrat']
concentrated6_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'Etomoxir']
concentrated7_df = pca_p0_df.loc[pca_p0_df['Treatment'] == 'K']
con_list = [concentrated0_df, concentrated2_df, concentrated5_df,
            concentrated6_df, concentrated7_df]
pca_p0_df = pd.concat(con_list, axis=0).reset_index(drop=True)

# 5. standardise the data
x = pca_p0_df.loc[:, 'C2:0': 'C12:1'].values  # Separating out the features
x = preprocessing.scale(x). # 使用中心标准化，即将变量都转化成z分数的形式，避免量纲问题对压缩造成影响
# 还可以使用StandardScale，也是在preprocessing里面

y = pca_p0_df.loc[:, ['Treatment']]  # Separating out the target

# 6. question component numbers of PCA（也就是估计多少个PC比较好）
pca = PCA(n_components=17)  # 1. 选择我想要多少个维度（如果只是测试那么就是所有的维度）
pca.fit(x)  # fit一次原始的数据??
principalComponents_17 = pca.fit_transform(x)  # 使用fit_transform 再来一次

# +++ 累积解释变异程度 (2 PC are enough)
plt.plot(np.cumsum(pca.explained_variance_ratio_), linewidth=3, alpha=.6)
plt.xlabel('Components / Factors')
plt.ylabel('Cumulative explained variance / Eigenvalue');
plt.grid(True)
plt.show()
pca_17 = pca.explained_variance_ratio_
pca_17_reshape = pca_17.reshape(-1, 1)  # -1表示任意行数，1表示1列

# 7.Plotting
# +++ preparation for generate the real PCA model (2 PC)
pca_real = PCA(n_components=2)  # 直接与变量个数相同的主成分
pca_real.fit(x)
principalComponents_2 = pca_real.fit_transform(x)  # fit_transform 表示将生成降维后的数据
pca_2 = pca_real.explained_variance_ratio_  # 这一步对于我观察我的PCA分离度很有帮助

# # 查看规模差别
# print("原始数据集规模:   ", x.shape)
# print("降维后的数据集规模:", pca_real_df.shape) #  这里可以看到数据被将为压缩了 
# 注意这里pca_real 要dataframe化 

# +++ for potential heatmap may needed
pca_real_df = pd.DataFrame(principalComponents_2)
pca_real_df.columns = ['pca_1', 'pca_2']
pca_real_df.index = y

plt.figure(figsize=(10, 12))
sns.set_theme(style='whitegrid', font='Avenir')
heatmap_ax = sns.heatmap(pca_real_df, linewidth=.5, linecolor='w', xticklabels=True,
                         cmap='twilight_shifted',
                         cbar_kws={'location': "right", 'shrink': 1, "fraction": .1, 'aspect': 20})
heatmap_ax.yaxis.set_tick_params(labelsize=10)
heatmap_ax.xaxis.set_tick_params(labelsize=10)
plt.title('Heatmap', fontsize=16)
plt.show()

# +++ for plotting PCA model
principalDf = pd.DataFrame(data=principalComponents_2, columns=['principal component 1', 'principal component 2'])
finalDf = pd.concat([principalDf, y], axis=1)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1)
# 和seaborn 不一样的作图思维以及作图方式，先构化框架
ax.set_xlabel(f"Principal Component 1 ({format(pca_2[0] * 100, '.2f')}%)", fontsize=15)  # 这里只取小数点后2位
ax.set_ylabel(f"Principal Component 2 ({format(pca_2[1] * 100, '.2f')}%)", fontsize=15)
ax.set_title('2 Component PCA', fontsize=16)

# Treatments = [ 'Xanthohumol','EGCG', 'Genistein',"Resveratrol" ]
# colors = [ 'y','m','C2','C9']
# Treatments = [ 'DMSO', 'K', 'Genistein', 'Etomoxir','Xanthohumol', 'EGCG', 'Bezafibrat', ]
# colors = ['w', 'k', 'c','m','y','r','b']
Treatments = ['DMSO', 'K', 'Etomoxir', 'Xanthohumol', 'Bezafibrat', ]
colors = ['w', 'C7', 'r', 'g', 'b']
for treatment, color in zip(Treatments, colors):
    indicesToKeep = finalDf['Treatment'] == treatment
    # 选择某个label下的数据进行绘制
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c=color, edgecolor='k', alpha=.85, linewidth=1
               , s=90)
# ax.grid()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210810_6h/NuG0_PCA.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
