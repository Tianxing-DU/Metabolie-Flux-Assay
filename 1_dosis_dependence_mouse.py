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
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Datentabelle_3T3-L1 2011-12-20_2014-01-27_corri.xls'
sheet_name = '080609ens1 NGS Ens'
mouse_df = pd.read_excel(file_path, sheet_name=sheet_name)  # open the correct sheet
mouse_df.columns = mouse_df.iloc[0]  # rename the columns
mouse_df = mouse_df.drop(index=0).reset_index(drop=True)  # delete the duplicate column names
mouse_df = mouse_df.iloc[0:144]
yaxislabel = 'Acylcarnitine products \n (normalized data on a log scale)'
xaxislabel = 'Oleic acid (µM)'

# +++ split data
mouse_df['Date'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[0])
mouse_df['Cell'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[1])
mouse_df['Time'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[2])
mouse_df['Repeat'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[3])
mouse_df['Charge'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[4])
mouse_df['Day'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[5])
mouse_df[xaxislabel] = mouse_df['Sample ID'].map(lambda x: x.split('_')[6])

# +++ purify data
main_df = mouse_df
main_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
main_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
main_df.replace(',', '.', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '100-C0:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = main_df[xaxislabel]  # extract the real interested columns for plotting

for x in dict_LOD:  # remove all values under LOD, and replace them with NaN
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN

# 2. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(mouse_df['CSA'], axis=0)  # for preparation of normalisation by CSA
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back
Simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1)

df = pd.melt(Simplified_df, id_vars=[xaxislabel], var_name='C-Chains', value_name=yaxislabel)

# 3. visualisation
dose_order = ['0', '0.3', '1', '3', '10', '30', '100', '300', '1000']
sns.set_theme(font_scale=1.6,style='ticks', font='Helvetica')

# +++ overview for having a feeling (not be used in publication)
dosis_ax = sns.catplot(x=xaxislabel, y=yaxislabel, data=df, kind="bar",
                       row='C-Chains', palette='summer_r', order=dose_order, row_order=list_custom,
                       errwidth=0.3, aspect=1.1, alpha=.8, dodge=True, edgecolor='black',
                       capsize=0.2, legend=True, sharey=False, sharex=False)
dosis_ax.set(yscale="log")
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211011_3h/0_2_dose_dependence_all_mouse.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# +++ explore the df in short, medium, long chain range (bar plot)(not be used in publication)
a = df.replace(to_replace=short_chain, value='short_chain')
b = a.replace(to_replace=medium_chain, value='medium_chain')
c = b.replace(to_replace=long_chain, value='long_chain')
grouped = c[yaxislabel].groupby([c[xaxislabel], c['C-Chains']]).sum()  # sum by using groupby
grouped_df = pd.DataFrame(grouped)
grouped_df.reset_index(inplace=True)  # convert the merged index to column

sublist = ['short_chain','medium_chain', 'long_chain']
dosis_ax = sns.catplot(x=xaxislabel, y=yaxislabel, data=grouped_df, kind="bar",
                  col='C-Chains', saturation=1, palette='Blues', col_order=sublist, order=dose_order,
                  aspect=1.2,  legend=False, sharey=False, sharex=False, edgecolor='k'
                 )
dosis_ax.set(yscale="log")
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211104/dose_dependence_chain_length_mouse.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# +++ plot some selected metabolite products (box plot)
sublist = ['C2:0', 'C10:0', 'C16:1']
dosis_ax = sns.catplot(x=xaxislabel, y=yaxislabel, data=df, kind="box",
                       col='C-Chains', saturation=1, palette='summer_r', col_order=sublist, order=dose_order, width=.6,
                       aspect=1.2, legend=False, sharey=False, sharex=False,
                       flierprops={'markerfacecolor': 'C7', 'markersize': 7, 'linestyle': 'none',
                                   'markeredgecolor': 'g'})
dosis_ax.set(yscale="log")
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211022_23_24_5h/2_2_dose_dependence_mouse_box.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
