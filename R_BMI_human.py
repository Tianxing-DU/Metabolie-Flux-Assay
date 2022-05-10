# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12

# 0. set up
from library_R import *

order = ['H2O', 'DMSO', 'Xanthohumol', 'Genistein', 'Resveratrol', 'EGCG']

# 1. DataFrame Processing:
# 1. DataFrame Processing:
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Datentabelle_humane Adipocyten 20110221_2014-01-27 corri.xls'
sheet_name = 'Behandlung'
human_df = pd.read_excel(file_path, sheet_name=sheet_name)  # open the correct sheet
human_df.columns = human_df.iloc[0]  # rename the columns
human_df = human_df.drop(index=0).reset_index(drop=True)  # delete the duplicate column names
yaxislabel = 'Acylcarnitine products \n (normalized data on a log scale)'
xaxislabel = 'Oleic acid (µM)'
# +++ split data
human_df['Date'] = human_df['Sample ID'].map(lambda x: x.split('_')[0])
human_df['Cell'] = human_df['Sample ID'].map(lambda x: x.split('_')[1])
human_df['BMI'] = human_df['Sample ID'].map(lambda x: x.split('_')[2])
human_df['Time'] = human_df['Sample ID'].map(lambda x: x.split('_')[3])
human_df['Repeat'] = human_df['Sample ID'].map(lambda x: x.split('_')[4])
human_df['Charge'] = human_df['Sample ID'].map(lambda x: x.split('_')[5])
human_df['Treatment'] = human_df['Sample ID'].map(lambda x: x.split('_')[6])

# +++ creat simplified and purified df for plotting
main_df = human_df
main_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
main_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
main_df.replace(',', '.', regex=True, inplace=True)
main_df.replace('Kontrolle', 'H2O', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '100-C0:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = main_df[['Date', 'Treatment', 'Charge', 'BMI']]  # extract the real interested columns for plotting

for x in dict_LOD:
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN

# 2. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # for preparation of normalisation by CSA
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back

simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1)
simplified_df = simplified_df.loc[simplified_df['Charge'].isin(['Oleic+'])]
# 3. set new colums with treatments-charge from df
melt_df = pd.melt(simplified_df, id_vars=['Date', 'Treatment', 'Charge', 'BMI'], var_name='C-Chains',
                  value_name='Normalised Value (mean)')

date_order = ['2010-06-03', '2010-07-01', '2010-12-09', '2010-07-15', '2010-11-25', '2011-02-10']
# melt_df = melt_df[~(melt_df.Date == '2010-12-09')]
# df_BMI_a_27 = melt_df.loc[melt_df['BMI'].isin(['>27'])].drop(columns=['Date', 'Charge', 'BMI'])
# df_BMI_u_27 = melt_df.loc[melt_df['BMI'].isin(['<27'])].drop(columns=['Date', 'Charge', 'BMI'])


melt_df = melt_df.loc[(melt_df.Date != date_order[3]) | (melt_df.Treatment != 'EGCG')]
df_BMI_a_27 = melt_df.loc[melt_df['BMI'].isin(['>27'])].drop(columns=['Date', 'Charge', 'BMI'])
df_BMI_u_27 = melt_df.loc[melt_df['BMI'].isin(['<27'])].drop(columns=['Date', 'Charge', 'BMI'])

# 6. groupby
# +++ calculate mean in set
grouped_a_27 = df_BMI_a_27['Normalised Value (mean)'].groupby(
    [df_BMI_a_27['Treatment'], df_BMI_a_27['C-Chains']]).mean()
grouped_u_27 = df_BMI_u_27['Normalised Value (mean)'].groupby(
    [df_BMI_u_27['Treatment'], df_BMI_u_27['C-Chains']]).mean()
# some programing thoughts for groupby here:
# 'Normalised Value (mean)'是我传入的第一个数据，或者说我想要分析的数据。这里也可以不用写，因为就这一个剩下
# 如果想要一次传入多个数组，后面的就是多个column的引用作为key索引
# 如果直接groupby，他不会'智慧的'按照词条进行默认的索引分组

# +++ reset index
grouped_a_27_df = pd.DataFrame(grouped_a_27).reset_index()
grouped_u_27_df = pd.DataFrame(grouped_u_27).reset_index()
# +++ prepare pivot table for heatmap
heat_a_27 = pd.pivot_table(grouped_a_27_df, values='Normalised Value (mean)', index=['Treatment'],
                           columns=['C-Chains'])
heat_u_27 = pd.pivot_table(grouped_u_27_df, values='Normalised Value (mean)', index=['Treatment'],
                           columns=['C-Chains'])
heat_a_27_recolumn = heat_a_27[list_custom]
heat_u_27_recolumn = heat_u_27[list_custom]
heat_a_27_organised = heat_a_27_recolumn.reindex(order, columns=list_custom)
heat_u_27_organised = heat_u_27_recolumn.reindex(order, columns=list_custom)
# sort index and columns

# 4. visualisation for comparison
# +++ separated data processing to calculate relative values
# effect due to treatments is compared with their corresponding ctrl.
# The relative values are presented through doing division

# etomoxir is controlled by K/H2O, so an extra division process is needed
# attention to the dividing consequences
dfs = [heat_a_27_organised, heat_u_27_organised]
for heat_table_organised in dfs:
    heat_table_organised.drop(axis=0, index='H2O', inplace=True)
    heat_table_organised.loc['EGCG'] = heat_table_organised.loc['EGCG'].div(heat_table_organised.loc['DMSO'])
    heat_table_organised.loc['Resveratrol'] = heat_table_organised.loc['Resveratrol'].div(
        heat_table_organised.loc['DMSO'])
    heat_table_organised.loc['Genistein'] = heat_table_organised.loc['Genistein'].div(
        heat_table_organised.loc['DMSO'])
    heat_table_organised.loc['Xanthohumol'] = heat_table_organised.loc['Xanthohumol'].div(
        heat_table_organised.loc['DMSO'])
    heat_table_organised.loc['DMSO'] = heat_table_organised.loc['DMSO'].div(
        heat_table_organised.loc['DMSO'])  # 这里开始有所修改
    # heat_table_organised.loc['H2O'] = heat_table_organised.loc['H2O'].div(heat_table_organised.loc['H2O'])

names = [f'Huamn BMI >27 without {date_order[3]}', f'Human BMI <27 without {date_order[3]}']
for heat_table, name in zip(dfs, names):
    sns.set_theme(style='ticks', font='Helvetica')
#     plt.figure(figsize=(12, 5))
    ax = sns.heatmap(heat_table, linewidth=.5, linecolor='w', square=True, xticklabels=True, vmin=-1,
                     vmax=3,
                     cmap='twilight_shifted',
                     cbar_kws={'location': "right", 'shrink': 1, "fraction": .2, 'aspect': 40})
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([-1, 1, 3])
    cbar.set_ticklabels(['Down-regulated', 'Neutral', 'Up-regulated'])
    plt.title(name)
    plt.tight_layout()
    # save_path = f'/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20220127_3h/{name}_heatmap.pdf'
    # plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
    plt.show()
# cmap=''RdYlBu_r'(更愉悦，中间黄色),'BrBG_r','RdBu_r（好看的红蓝分布）','twilight_shifted（低沉有魅力）'
