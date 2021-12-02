# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12

# 0.1 set up
from library_new_version import *

order = ['H2O', 'DMSO', 'Xanthohumol', 'Genistein', 'Resveratrol', 'EGCG']

# 1. DataFrame Processing:
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Datentabelle_humane Adipocyten 20110221_2014-01-27 corri.xls'
sheet_name = 'Behandlung'
human_df = pd.read_excel(file_path, sheet_name=sheet_name)  # open the correct sheet
human_df.columns = human_df.iloc[0]  # rename the columns
human_df = human_df.drop(index=0).reset_index(drop=True)  # delete the duplicate column names

human_df['Date'] = human_df['Sample ID'].map(lambda x: x.split('_')[0])
human_df['Cell'] = human_df['Sample ID'].map(lambda x: x.split('_')[1])
human_df['Unknown'] = human_df['Sample ID'].map(lambda x: x.split('_')[2])
human_df['Time'] = human_df['Sample ID'].map(lambda x: x.split('_')[3])
human_df['Repeat'] = human_df['Sample ID'].map(lambda x: x.split('_')[4])
human_df['Charge'] = human_df['Sample ID'].map(lambda x: x.split('_')[5])
human_df['Treatment'] = human_df['Sample ID'].map(lambda x: x.split('_')[6])

# 2. creat simplified and purified df for plotting
main_df = human_df
main_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
main_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
main_df.replace(',', '.', regex=True, inplace=True)
main_df.replace('Kontrolle', 'H2O', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '100-C0:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = main_df[['Treatment', 'Charge']]
for x in dict_LOD:
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN

# 3. processing data
# +++ normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # for preparation of normalisation by CSA

# +++ merge the int and str df
simplified_df = pd.concat([Extract_str_df, Extract_int_df_div], axis=1)

# +++ set new colums with treatments-charge from df
melt_df = pd.melt(simplified_df, id_vars=['Treatment', 'Charge'], var_name='C-Chains',
                  value_name='Normalised Value (mean)')
pair_df = pd.DataFrame(melt_df['Treatment'].str.cat(melt_df['Charge'], sep='&')).rename(
    columns={'Treatment': 'Pair'})  # pair 'Treatment' and 'Charge
final_df = pd.concat([melt_df, pair_df], axis=1).loc[:, ['C-Chains', 'Normalised Value (mean)', 'Pair']]
# all pairs / infos for plotting are included

# +++ groupby and calculate mean in set
grouped = final_df['Normalised Value (mean)'].groupby([final_df['Pair'], final_df['C-Chains']]).mean()
# some programing thoughts for groupby here:
# 'Normalised Value (mean)'是我传入的第一个数据，或者说我想要分析的数据。这里也可以不用写，因为就这一个剩下
# 如果想要一次传入多个数组，后面的就是多个column的引用作为key索引
# 如果直接groupby，他不会'智慧的'按照词条进行默认的索引分组

grouped_df = pd.DataFrame(grouped)
# some programing thoughts: 不需要使用这个把index变成列，下面的代码reset的时候，直接就分配到了列。
# alternative：grouped_df['Treatment'] = grouped_df.index
grouped_df.reset_index(inplace=True)

# +++ build abundance df
abundance_df = pd.merge(grouped_df, df_LOD, on='C-Chains').dropna()
# dropna removes C0:0, although later pivot function will remove C0:)'s NaN value anyway: So here dropna not a must
# some programing thoughts
# 这里其实可以不用on，默认相同序列，但是我这里写出来，方便复杂的情况下一次使用
abundance_df['Treatment'] = abundance_df['Pair'].map(lambda x: x.split('&')[0])
abundance_df['Charge'] = abundance_df['Pair'].map(lambda x: x.split('&')[1])
abundance_df.drop('Pair', axis=1, inplace=True)

# +++ save abundance info separately in excel
# abundance_high_df = abundance_df.loc[(abundance_df['Normalised Value (mean)'] >= 10 * abundance_df['LOD'])]
# abundance_medium_df = abundance_df.loc[(abundance_df['Normalised Value (mean)'] >= 4.5 * abundance_df['LOD']) & (
#         abundance_df['Normalised Value (mean)'] < 10 * abundance_df['LOD'])]
# abundance_low_df = abundance_df.loc[(abundance_df['Normalised Value (mean)'] >= 1 * abundance_df['LOD']) & (
#         abundance_df['Normalised Value (mean)'] < 4.5 * abundance_df['LOD'])]
# abundance_below_LOD_df = abundance_df.loc[(abundance_df['Normalised Value (mean)'] < 1 * abundance_df['LOD'])]
#
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210930_5h/abundance_all_treatments_human.xlsx'
# with pd.ExcelWriter(save_path) as writer:
#     abundance_below_LOD_df.to_excel(writer, 'below_LOD')
#     abundance_low_df.to_excel(writer, 'low')
#     abundance_medium_df.to_excel(writer, 'medium')
#     abundance_high_df.to_excel(writer, 'high')
# some programming thoughts here
# 多个dataframe需要写入同一个excel时，每次使用df.to_excel(文件名)的形式去写，系统都会重新创建一个新的文件。也就意味着前面的文件会被覆盖掉，你得到的只能是最后一个df写入的结果文件

abundance_df['Abundance'] = abundance_df['Normalised Value (mean)'] / abundance_df['LOD']
heat_oleic_yes_df = abundance_df.loc[abundance_df['Charge'] == 'Oleic+']
heat_oleic_no_df = abundance_df.loc[abundance_df['Charge'] == 'Oleic-']
heat_table = pd.pivot_table(heat_oleic_yes_df, values='Abundance', index=['Treatment'], columns=['C-Chains'])

# 4. visualisation
heat_table_recolumn = heat_table[list_custom]
heat_comparison_df = heat_table_recolumn.reindex(order)
# heat_comparison_df = heat_table_organised.apply(np.log)
# log function is needed to have a better visualisation
# !!! but log here is not preferred by journals, so I decide to change the colorbar on a log scale instead of data
# heat_comparison_percentage_df = heat_comparison_df.div(heat_comparison_df.stack().max())
# at the beginning, I prefer to do in %

# +++ calculation for colour bar
# min_value = heat_comparison_df.stack().min()
# max_value = heat_comparison_df.stack().max()
# zw_value = max_value - min_value
# if use the log above, I can set color bar in % (or to say, divide my colorbar properly, and define 'low', 'medium','high'...

# +++ plotting
sns.set_theme(style='ticks', font='Helvetica')
plt.figure(figsize=(12, 5))
ax = sns.heatmap(heat_comparison_df, linewidth=.5, linecolor='w', square=True, xticklabels=True,
                 cmap='YlGnBu', norm=LogNorm(vmin=100, vmax=100000),
                 # !!! to set a shared log colorbar, must set vmin and vmax inside LogNorm, NOT ouside!!!
                 # with this code here and same in the code of mice analysis, we can make two heatmap share a same colorbar
                 cbar_kws={'location': "right", 'shrink': 1, "fraction": .2, 'aspect': 40})
# function LogNorm here is essential, better than using the original colorbar arguments
# cbar_kws={'ticks':MaxNLocator(2), 'format':'%.e'}
cbar = ax.collections[0].colorbar
plt.tight_layout()
plt.title('Human')
cbar.set_ticks([100000, 10000, 1000, 100])
# cbar.set_ticks([min_value+.25 * zw_value, min_value+.5 * zw_value, min_value+.75 * zw_value, max_value])
cbar.set_ticklabels(['High \n$10^5$ * LOD ', 'Medium \n$10^4$ * LOD', 'Low \n$10^3$ * LOD', 'Very low \n$10^2$ * LOD', ])
# $ $ makes the exponential superscript available
save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211130_3h/abundance_human_v1130.pdf'
plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
