# 0. import
from library_new_version import *

# 1. data import and pre-processing:
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Datentabelle_3T3-L1 2011-12-20_2014-01-27_corri.xls'
sheet_name = 'Zeitabhängigkeit'
mouse_df = pd.read_excel(file_path, sheet_name=sheet_name)  # open the correct sheet
mouse_df.columns = mouse_df.iloc[0]  # rename the columns
mouse_df = mouse_df.drop(index=0).reset_index(drop=True)  # delete the duplicate column names
yaxislabel = 'Acylcarnitine products \n (normalized data on a log scale)'
xaxislabel = 'Incubation Time (min)'

# +++ split the column 'Sample ID' into multiple columns with '_'
mouse_df['Date'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[0])
mouse_df['Cell'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[1])
mouse_df[xaxislabel] = mouse_df['Sample ID'].map(lambda x: x.split('_')[2])
mouse_df['Repeat'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[3])
mouse_df['Charge'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[4])
mouse_df['Day'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[5])
mouse_df['Conc (charge)'] = mouse_df['Sample ID'].map(lambda x: x.split('_')[6])

# +++ creat simplified and purified df for plotting
mouse_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
mouse_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
mouse_df.replace(',', '.', regex=True, inplace=True)
Extract_int_df = mouse_df.loc[:, '100-C0:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = mouse_df[[xaxislabel, 'Charge']]

for x in dict_LOD:
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN

# 2. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(mouse_df['CSA'], axis=0)  # for preparation of normalisation by CSA
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back
Simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1)
df_melt = pd.melt(Simplified_df, id_vars=[xaxislabel, 'Charge'], var_name='C-Chains',
                  value_name=yaxislabel)
# 3. visualisation
df_melt = pd.melt(Simplified_df, id_vars=[xaxislabel, 'Charge'], var_name='C-Chains', value_name=yaxislabel)

time_order = ['30', '60', '90', '120']

sns.set_theme(font_scale=1.6, style='ticks', font='Helvetica')

# +++ plot the data with hue='Time', row='Charge' x=C0 to C18, to have an overview
ax = sns.catplot(x=xaxislabel, y=yaxislabel, data=df_melt, kind="bar", hue="Charge", hue_order=['Oleic+'],
                 row='C-Chains', palette='hot', order=time_order, row_order=list_custom,
                 errwidth=0.3, aspect=1.1, alpha=.8, dodge=True, edgecolor='black',
                 capsize=0.2, legend=True, sharey=False, sharex=False)
ax.set(yscale="log")
ax.set_xticklabels(rotation=90)
plt.legend(title='Time', loc=0)
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20210713/time_pointplot.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# +++ plot some selected metabolite products (box plot)
plt.figure(figsize=(10, 6))
sublist = ['C2:0', 'C8:0', 'C16:1']
ax1 = sns.catplot(x=xaxislabel, y=yaxislabel, data=df_melt.loc[df_melt['Charge'] == 'Oleic+'], kind="box",
                  col='C-Chains', saturation=1, palette='Blues', col_order=sublist, order=time_order, width=.6,
                  aspect=.7, legend=False, sharey=False, sharex=False,
                  flierprops={'markerfacecolor': 'C7', 'markersize': 7, 'linestyle': 'none', 'markeredgecolor': 'b'})
ax1.set(yscale="log")
plt.tight_layout()
# save_path = "/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211022_23_24_5h/1_3_time_dependence_mouse_box.pdf"
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
