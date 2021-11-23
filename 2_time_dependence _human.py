# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12

# 0. import
import matplotlib
from library_new_version import *

# 1. creat simplified and purified df for plotting
# +++ import data
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Datentabelle_humane Adipocyten 20110221_2014-01-27 corri.xls'
sheet_name = 'Zeitabhängigkeit'
human_df = pd.read_excel(file_path, sheet_name=sheet_name)  # open the correct sheet
human_df.columns = human_df.iloc[0]  # rename the columns
human_df = human_df.drop(index=0).reset_index(drop=True)  # delete the duplicate column names
yaxislabel = 'Acylcarnitine products \n (normalized data on a log scale)'
xaxislabel = 'Incubation Time (s)'

# +++ split data
human_df['Date'] = human_df['Sample ID'].map(lambda x: x.split('_')[0])
human_df['Cell'] = human_df['Sample ID'].map(lambda x: x.split('_')[1])
human_df[xaxislabel] = human_df['Sample ID'].map(lambda x: x.split('_')[2])  # varied incubation time
human_df['Repeat'] = human_df['Sample ID'].map(lambda x: x.split('_')[3])
human_df['Charge'] = human_df['Sample ID'].map(lambda x: x.split('_')[4])

# +++ purify data
main_df = human_df  # I do not change the name 'main.df' here to save time, as codes copied from mouse data processing
main_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
main_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
main_df.replace(',', '.', regex=True, inplace=True)
main_df.replace('Kontrolle', 'H2O', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '100-C0:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = main_df[[xaxislabel, 'Charge']]  # extract the real interested columns for plotting

for x in dict_LOD:
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN

# 2. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # for preparation of normalisation by CSA
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back
Simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1)
df_melt = pd.melt(Simplified_df, id_vars=[xaxislabel, 'Charge'], var_name='C-Chains', value_name=yaxislabel)

# 3. visualisation
# +++ plot the data with hue='Time', row='Charge' x=C0 to C18, overview for having a feeling (not be used in publication)
time_order = ['30', '60', '90', '120']

sns.set_theme(font_scale=1.6, style='ticks', font='Helvetica')
ax = sns.catplot(x=xaxislabel, y=yaxislabel, data=df_melt, kind="bar", hue="Charge", hue_order=['Oleic+'],
                 row='C-Chains', palette='hot', order=time_order, row_order=list_custom,
                 errwidth=0.3, aspect=1.1, alpha=.8, dodge=True, edgecolor='black',
                 capsize=0.2, legend=True, sharey=False, sharex=False)
ax.set(yscale="log")
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211011_3h/0_0_time_all_human.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# +++ bar plot (select representative metabolite products)(not be used in publication)
plt.figure(figsize=(12, 10))
sns.set_theme(font_scale=1.6, style='ticks', font='Helvetica')
sublist = ['C6:0', 'C12:0', 'C18:1']
ax1 = sns.catplot(x=xaxislabel, y=yaxislabel, data=df_melt, kind="bar", hue='Charge', hue_order=['Oleic+'],
                  col='C-Chains', palette='hot', col_order=sublist, order=time_order,
                  errwidth=0.8, aspect=.9, alpha=.7, edgecolor='black', legend=False,
                  capsize=0.15, sharey=False, sharex=False)
plt.legend(title='Charge', bbox_to_anchor=(1, 0.6), frameon=False, fontsize='large', title_fontsize='large')
ax1.set(yscale="log")
plt.tight_layout()
# save_path = "/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211011_3h/0_1_time_dependence_human.pdf"
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# +++ boxplot (same as above but in box plot)
plt.figure(figsize=(10, 6))
sns.set_theme(font_scale=1.6, style='ticks', font='Helvetica')
sublist = ['C2:0', 'C8:0', 'C16:1']
ax1 = sns.catplot(x=xaxislabel, y=yaxislabel, data=df_melt.loc[df_melt['Charge'] == 'Oleic+'], kind="box",
                  col='C-Chains', saturation=1, palette='Blues', col_order=sublist, order=time_order, width=.6,
                  aspect=.7, legend=False, sharey=False, sharex=False,
                  flierprops={'markerfacecolor': 'C7', 'markersize': 7, 'linestyle': 'none', 'markeredgecolor': 'b'})
ax1.set(yscale="log")
plt.tight_layout()
# save_path = "/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211022_23_24_5h/0_1_time_dependence_human_box.pdf"
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
