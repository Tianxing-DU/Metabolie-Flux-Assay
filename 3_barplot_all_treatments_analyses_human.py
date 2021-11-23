# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12


# 0. set up infos

from library_new_version import *
order = ['H2O','DMSO', 'Resveratrol', 'Genistein', 'Xanthohumol', 'EGCG','Bezafibrat', 'Etomoxir']

# 1. DataFrame Processing:
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Behandlung_3T3-L1 2011-12-20_2014-01-27_corr.xls'
main_df = pd.read_excel(file_path)

# +++ creat simplified and purified df for plotting
main_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
main_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
main_df.replace(',', '.', regex=True, inplace=True)
main_df.replace('K', 'H2O', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '100-C0:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = main_df[['Treatment', 'Charge']]  # extract the real interested columns for plotting

# 2. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(main_df['CSA'], axis=0)  # for preparation of normalisation by CSA
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back
simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1) # merge the int and str df

# 3. set new columns with treatments-charge from df
melt_df = pd.melt(simplified_df, id_vars=['Treatment', 'Charge'], var_name='C-Chains', value_name='Normalised Value (mean)')
pair_df = pd.DataFrame(melt_df['Treatment'].str.cat(melt_df['Charge'], sep='-')).rename(columns={'Treatment': 'Pair'})  # pair 'Treatment' and 'Charge
final_df = pd.concat([melt_df, pair_df], axis=1).loc[:, ['C-Chains', 'Normalised Value (mean)', 'Pair']]  # all pairs / infos for plotting are included


# 4. visualisation
"""Xanthohumol_DMSO_with_oleic"""
sns.set_theme( style='ticks', font='Helvetica',font_scale=2)
plt.figure(figsize=(10, 21))
Xanthohumol_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair',orient='h',
                             hue_order=['Xanthohumol-Oleic+', 'DMSO-Oleic+'],
                             order=list_custom, palette='spring')
Xanthohumol_ax.set(xscale="log")
plt.legend(title='Treatment')
Xanthohumol_ax.get_legend().remove()
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/XN_DMSO_mouse.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""EGCG_DMSO_with_oleic"""
plt.figure(figsize=(10, 21))
EGCG_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair',orient='h',
                      hue_order=['EGCG-Oleic+', 'DMSO-Oleic+'],
                      order=list_custom, palette='spring')
EGCG_ax.set(xscale="log")
plt.legend(title='Treatment')
EGCG_ax.get_legend().remove()
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/EGCG_DMSO_mouse.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""Genistein_DMSO_with_oleic"""
plt.figure(figsize=(10, 21))
Genistein_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair',orient='h',
                           hue_order=['Genistein-Oleic+', 'DMSO-Oleic+'],
                           order=list_custom, palette='spring')
Genistein_ax.set(xscale="log")
plt.legend(title='Treatment')
Genistein_ax.get_legend().remove()
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/Genistein_DMS_mouse.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""Resveratrol_DMSO_with_oleic"""
plt.figure(figsize=(10, 21))
Resveratrol_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair',orient='h',
                             hue_order=['Resveratrol-Oleic+', 'DMSO-Oleic+'],
                             order=list_custom, palette='spring')
Resveratrol_ax.set(xscale="log")
plt.legend(title='Treatment')
Resveratrol_ax.get_legend().remove()
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/Resveratrol_DMSO_mouse.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""H2O_DMSO_with_oleic"""
plt.figure(figsize=(10, 21))
Kontrolle_ax = sns.barplot(data=final_df, y='C-Chains', x='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair',orient='h',
                           hue_order=['DMSO-Oleic+', 'H2O-Oleic+'],  # here H2O works as control for DMSO treatment
                           order=list_custom, palette='spring')
Kontrolle_ax.set(xscale="log")
Kontrolle_ax.get_legend().remove()
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/DMSO_H2O_mouse.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""Bezafibrat_DMSO_with_oleic"""

plt.figure(figsize=(21, 10))
Resveratrol_ax = sns.barplot(data=final_df, x='C-Chains', y='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair',
                             hue_order=['Bezafibrat-Oleic+', 'DMSO-Oleic+'],
                             order=list_custom, palette='spring')
Resveratrol_ax.set(yscale="log")
plt.xticks(rotation=90)
plt.legend(title='Treatment')
Resveratrol_ax.get_legend().remove()
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/Bezafibrat_DMSO_mouse.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()

"""Etomoxir_H2O_with_oleic"""
plt.figure(figsize=(21, 10))
Resveratrol_ax = sns.barplot(data=final_df, x='C-Chains', y='Normalised Value (mean)', capsize=0.3, errwidth=0.8,
                             alpha=.6, edgecolor='black', dodge=False, hue='Pair',
                             hue_order=['Etomoxir-Oleic+', 'H2O-Oleic+'],
                             order=list_custom, palette='spring')
Resveratrol_ax.set(yscale="log")
plt.xticks(rotation=90)
plt.legend(title='Treatment')
Resveratrol_ax.get_legend().remove()
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211021_5h/Etomoxir_H2O_mouse.pdf'
# plt.savefig(save_path, dpi=1000, bbox_inches='tight', transparent=True)
plt.show()
