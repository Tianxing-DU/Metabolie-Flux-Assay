# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12

# 0. import
from library_new_version import *

# 1. creat simplified and purified df for plotting
# +++ import data
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Datentabelle_humane Adipocyten 20110221_2014-01-27 corri.xls'
sheet_name = 'Konzentrationsabhängigkeit'
human_df = pd.read_excel(file_path, sheet_name=sheet_name)  # open the correct sheet
human_df.columns = human_df.iloc[0]  # rename the columns
human_df = human_df.drop(index=0).reset_index(drop=True)  # delete the duplicate column names
yaxislabel = 'Acylcarnitine products \n (normalized data on a log scale)'  # for preparing axis names in plots later
xaxislabel = 'Oleic acid (µM)'

# +++ split data
human_df['Date'] = human_df['Sample ID'].map(lambda x: x.split('_')[0])
human_df['Cell'] = human_df['Sample ID'].map(lambda x: x.split('_')[1])
human_df['Time'] = human_df['Sample ID'].map(lambda x: x.split('_')[2])
human_df['Repeat'] = human_df['Sample ID'].map(lambda x: x.split('_')[3])
human_df['Charge'] = human_df['Sample ID'].map(lambda x: x.split('_')[4])
human_df[xaxislabel] = human_df['Sample ID'].map(lambda x: x.split('_')[5])  # varied concentrations

# +++ purify data
main_df = human_df  # I do not change the name 'main.df' here to save time, as codes copied from mouse data processing
main_df.replace(' ', '', regex=True, inplace=True)
main_df.replace('µM', '', regex=True, inplace=True)
main_df.replace(',', '.', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '101-C2:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # simplify the column name
Extract_str_df = main_df[xaxislabel]  # extract the real interested columns for plotting
for x in dict_LOD:  # remove all values under LOD, and replace them with NaN
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN
# alternative way here for selection dict_sel = ({k: v} for k, v in dict_LOD.items() if v > 0.005) 

# 2. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back
Simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1)
df = pd.melt(Simplified_df, id_vars=[xaxislabel], var_name='C-Chains', value_name=yaxislabel)  # data is finally ready

# 3. visualisation
# +++ overview for having a feeling (not be used in publication)
dosis_ax = sns.catplot(x=xaxislabel, y=yaxislabel, data=df, kind="bar",
                       row='C-Chains', palette='summer_r', order=dose_order, row_order=list_custom,
                       errwidth=0.3, aspect=1.1, alpha=.8, dodge=True, edgecolor='black',
                       capsize=0.2, legend=True, sharey=False, sharex=False)
dosis_ax.set(yscale="log")
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210902_6h/dosis_all.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# +++ explore the df in short, medium, long chain range (bar plot)(not be used in publication)
a = df.replace(to_replace=short_chain, value='short_chain')
b = a.replace(to_replace=medium_chain, value='medium_chain')
c = b.replace(to_replace=long_chain, value='long_chain')
# some programming thoughts:
# 这里是通过replace的手段，替换metabolite 名字--> short、long chain
# 而下面的这种alternative的方式，是独立提取我想要的部分，通过函数 isin
# short_chain_df = df[df['C-Chains'].isin(short_chain)]
# middle_chain_df = df[df['C-Chains'].isin(middle_chain)]
# long_chain_df = df[df['C-Chains'].isin(long_chain)]
# 这里我尝试了很多次drop的 for loop都没有搞定，其实isin一下就搞定了

grouped = c[yaxislabel].groupby([c[xaxislabel], c['C-Chains']]).sum()  # sum by using groupby
grouped_df = pd.DataFrame(grouped)  # transfer series to df
grouped_df.reset_index(inplace=True)  # convert the merged index to column

dose_order = ['0', '0.3', '1', '3', '10', '30', '100', '300', '1000']  # x axis

sns.set_theme(font_scale=1.6, style='ticks', font='Helvetica')
sublist = ['short_chain', 'medium_chain', 'long_chain']
dosis_ax = sns.catplot(x=xaxislabel, y=yaxislabel, data=grouped_df, kind="bar",
                       col='C-Chains', saturation=1, palette='Blues', col_order=sublist, order=dose_order,
                       edgecolor='k', aspect=1.2, legend=False, sharey=False, sharex=False)
dosis_ax.set(yscale="log")
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211104/dose_dependence_chain_length_human.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# +++ plot some selected metabolite products (box plot)
sns.set_theme(font_scale=1.6,style='ticks', font='Helvetica')
sublist = ['C2:0','C10:0', 'C16:1']
dosis_ax = sns.catplot(x=xaxislabel, y=yaxislabel, data=df, kind="box",
                  col='C-Chains', saturation=1, palette='summer_r', col_order=sublist, order=dose_order,width=.6,
                  aspect=1.2,  legend=False, sharey=False, sharex=False,
                  flierprops={'markerfacecolor':'C7', 'markersize':7,'linestyle':'none', 'markeredgecolor':'g'})
dosis_ax.set(yscale="log")
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211022_23_24_5h/2_2_dose_dependence_human.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
