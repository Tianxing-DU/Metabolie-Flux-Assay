# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12

# 0. set up
from library_new_version import *
order = ['H2O', 'DMSO', 'Resveratrol', 'Genistein', 'Xanthohumol', 'EGCG', 'Bezafibrat', 'Etomoxir']
order_mainbody = ['H2O', 'DMSO', 'EGCG', 'Resveratrol', 'Genistein', 'Xanthohumol']
order_supplyment = ['H2O', 'DMSO', 'Bezafibrat', 'Etomoxir']

# 1. data input:
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Behandlung_3T3-L1 2011-12-20_2014-01-27_corr.xls'
main_df = pd.read_excel(file_path)

# 2. creat simplified and purified df for plotting
main_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
main_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
main_df.replace(',', '.', regex=True, inplace=True)
main_df.replace('K', 'H2O', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '100-C0:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = main_df[['Treatment', 'Charge']]
# extract the real interested columns for plotting

for x in dict_LOD:
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN

# 3. processing data
# +++ normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # for preparation of normalisation by CSA
# reasons here without normalisation is demonstrated later downside

# +++ merge the int and str df
simplified_df = pd.concat([Extract_str_df, Extract_int_df_div], axis=1)

# +++ set new colums with treatments-charge from df
melt_df = pd.melt(simplified_df, id_vars=['Treatment', 'Charge'], var_name='C-Chains',value_name='Normalised Value (mean)')
pair_df = pd.DataFrame(melt_df['Treatment'].str.cat(melt_df['Charge'], sep='&')).rename(columns={'Treatment': 'Pair'})  # pair 'Treatment' and 'Charge
final_df = pd.concat([melt_df, pair_df], axis=1).loc[:, ['C-Chains', 'Normalised Value (mean)', 'Pair']]
# all pairs / infos for plotting are included

# +++ groupby
# +++ calculate mean in set
grouped = final_df['Normalised Value (mean)'].groupby([final_df['Pair'], final_df['C-Chains']]).mean()
# some programing thoughts for groupby here:
# 'Normalised Value (mean)'是我传入的第一个数据，或者说我想要分析的数据。这里也可以不用写，因为就这一个剩下
# 如果想要一次传入多个数组，后面的就是多个column的引用作为key索引
# 如果直接groupby，他不会'智慧的'按照词条进行默认的索引分组

grouped_df = pd.DataFrame(grouped)
# some programing thoughts: 不需要使用这个把index变成列，下面的代码reset的时候，直接就分配到了列。
# alternative：grouped_df['Treatment'] = grouped_df.index
grouped_df.reset_index(inplace=True)
# dropnan removes here C0:0 (which we are not interested)

grouped_df['Treatment'] = grouped_df['Pair'].map(lambda x: x.split('&')[0])
grouped_df['Charge'] = grouped_df['Pair'].map(lambda x: x.split('&')[1])
grouped_df.drop('Pair', axis=1, inplace=True)

# +++ select oleic+
heat_oleic_yes_df = grouped_df.loc[grouped_df['Charge'] == 'Oleic+']  # for the moment, only oleic+ are interested
heat_oleic_no_df = grouped_df.loc[grouped_df['Charge'] == 'Oleic-']

# +++ prepare pivot table for heatmap
heat_table = pd.pivot_table(heat_oleic_yes_df, values='Normalised Value (mean)', index=['Treatment'], columns=['C-Chains'])
heat_table_recolumn = heat_table[list_custom]
heat_table_organised = heat_table_recolumn.reindex(order)
# sort index and columns

# 4. visualisation for comparison
# +++ separated data processing to calculate relative values
# effect due to treatments is compared with their corresponding ctrl.
# The relative values are presented through doing division

# etomoxir is controlled by K/H2O, so an extra division process is needed
# attention to the dividing consequences
heat_table_organised.loc['EGCG'] = heat_table_organised.loc['EGCG'].div(heat_table_organised.loc['DMSO'])
heat_table_organised.loc['Resveratrol'] = heat_table_organised.loc['Resveratrol'].div(heat_table_organised.loc['DMSO'])
heat_table_organised.loc['Genistein'] = heat_table_organised.loc['Genistein'].div(heat_table_organised.loc['DMSO'])
heat_table_organised.loc['Xanthohumol'] = heat_table_organised.loc['Xanthohumol'].div(heat_table_organised.loc['DMSO'])
heat_table_organised.loc['Bezafibrat'] = heat_table_organised.loc['Bezafibrat'].div(heat_table_organised.loc['DMSO'])

heat_table_organised.loc['Etomoxir'] = heat_table_organised.loc['Etomoxir'].div(heat_table_organised.loc['H2O'])
heat_table_organised.loc['DMSO'] = heat_table_organised.loc['DMSO'].div(heat_table_organised.loc['H2O'])
heat_table_organised.loc['H2O'] = heat_table_organised.loc['H2O'].div(heat_table_organised.loc['H2O'])

heatmap_mainbody_df = heat_table_organised.loc[order_mainbody]
heatmap_supplyment_df = heat_table_organised.loc[order_supplyment]

# heat_comparison_df = heat_comparison_df.apply(np.log)
# WHY no LOGARITHMS: it overtones the down-regulated abundance (mathematics)

# x = heat_comparison_df  # returns a numpy array
# min_max_scaler = preprocessing.MinMaxScaler()
# x_scaled = min_max_scaler.fit_transform(x)
# heat_comparison_df_normalised = pd.DataFrame(x_scaled, columns=heat_comparison_df.columns, index=heat_comparison_df.index).dropna(axis=1)
# WHY no MINMAX-NORMALISATION: not satisfied visualisation after doing this

# +++ plot for mainbody
sns.set_theme(style='ticks', font='Helvetica')
plt.figure(figsize=(12, 5))
ax = sns.heatmap(heatmap_mainbody_df, linewidth=.5, linecolor='w', square=True, xticklabels=True, vmin=-1, vmax=3,
                 cmap='twilight_shifted',
                 cbar_kws={'location': "right", 'shrink': .5, "fraction": .2, 'aspect': 20})
cbar = ax.collections[0].colorbar
cbar.set_ticks([-1, 1, 3])
cbar.set_ticklabels(['Up-regulated', 'Neutral', 'Down-regulated'])
plt.title('Mouse')
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/regulation_overview_heatmap_twilight_try9_mouse_mainbody.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# +++ plot for supplyment
plt.figure(figsize=(12, 5))
ax = sns.heatmap(heatmap_supplyment_df, linewidth=.5, linecolor='w', square=True, xticklabels=True, vmin=-1, vmax=3,
                 cmap='twilight_shifted',
                 cbar_kws={'location': "right", 'shrink': .5, "fraction": .2, 'aspect': 20})
cbar = ax.collections[0].colorbar
cbar.set_ticks([-1, 1, 3])
cbar.set_ticklabels(['Up-regulated', 'Neutral', 'Down-regulated'])
plt.title('Mouse')
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/regulation_overview_heatmap_twilight_try9_mouse_supplyment.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
# some colour choosing thoughts: cmap='RdYlBu_r'(更愉悦，中间黄色),'BrBG_r','RdBu_r（好看的红蓝分布）','twilight_shifted（低沉有魅力）'
