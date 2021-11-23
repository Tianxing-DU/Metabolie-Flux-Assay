# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12

# 0. set up infos

from library_new_version import *
order = ['H2O','DMSO', 'EGCG', 'Resveratrol', 'Genistein', 'Xanthohumol']

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
human_df['Unknown'] = human_df['Sample ID'].map(lambda x: x.split('_')[2])
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
Extract_str_df = main_df[['Treatment', 'Charge']]  # extract the real interested columns for plotting

for x in dict_LOD:
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN

# 2. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # for preparation of normalisation by CSA
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back

simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1)

# 3. set new colums with treatments-charge from df
melt_df = pd.melt(simplified_df, id_vars=['Treatment', 'Charge'], var_name='C-Chains', value_name='Normalised Value (mean)')
pair_df = pd.DataFrame(melt_df['Treatment'].str.cat(melt_df['Charge'], sep='-')).rename(columns={'Treatment': 'Pair'})
# pair 'Treatment' and 'Charge

final_df = pd.concat([melt_df, pair_df], axis=1).loc[:, ['C-Chains', 'Normalised Value (mean)', 'Pair']]
# all pairs / infos for plotting are included

# 4. barplot plotting

"""Xanthohumol_DMSO_with_oleic"""
sns.set_theme( style='ticks', font='Helvetica',font_scale=2)
plt.figure(figsize=(10, 21))
Xanthohumol_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair', orient='h',
                             hue_order=['Xanthohumol-Oleic+', 'DMSO-Oleic+'],
                             order=list_custom, palette='spring')
Xanthohumol_ax.set(xscale="log")
plt.legend(title='Treatment')
Xanthohumol_ax.get_legend().remove()
plt.tight_layout()
# plt.title('Xanthohumol_human', fontsize='large', fontweight='bold')
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/XN_DMSO_human.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""EGCG_DMSO_with_oleic"""

plt.figure(figsize=(10, 21))
EGCG_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair', orient='h',
                      hue_order=['EGCG-Oleic+', 'DMSO-Oleic+'],
                      order=list_custom, palette='spring')
EGCG_ax.set(xscale="log")
plt.legend(title='Treatment')
EGCG_ax.get_legend().remove()
plt.tight_layout()
# plt.title('EGCG', fontsize='large', fontweight='bold')
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/EGCG_DMSO_human.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""Genistein_DMSO_with_oleic"""

plt.figure(figsize=(10, 21))
Genistein_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair', orient='h',
                           hue_order=['Genistein-Oleic+', 'DMSO-Oleic+'],
                           order=list_custom, palette='spring')
Genistein_ax.set(xscale="log")
plt.legend(title='Treatment')
Genistein_ax.get_legend().remove()
plt.tight_layout()
# plt.title('Genistein', fontsize='large', fontweight='bold')
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/Genistein_DMSO_human.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""Resveratrol_DMSO_with_oleic"""
plt.figure(figsize=(10, 21))
Resveratrol_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair', orient='h',
                             hue_order=['Resveratrol-Oleic+', 'DMSO-Oleic+'],
                             order=list_custom, palette='spring')
Resveratrol_ax.set(xscale="log")
plt.legend(title='Treatment')
Resveratrol_ax.get_legend().remove()
plt.tight_layout()
# plt.title('Resveratrol', fontsize='large', fontweight='bold')
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/Resveratrol_DMSO_human.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""H2O_DMSO_with_oleic"""
plt.figure(figsize=(10, 21))
Kontrolle_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair', orient='h',
                           hue_order=['DMSO-Oleic+', 'H2O-Oleic+'],
                           order=list_custom, palette='spring')
Kontrolle_ax.set(xscale="log")
Kontrolle_ax.get_legend().remove()
plt.tight_layout()
# plt.title('DMSO vs H2O', fontsize='large', fontweight='bold')
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/DMSO_H2O_human.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()
