# python           : 3.9.0.final.0
# pandas           : 1.2.4
# numpy            : 1.20.3
# matplotlib       : 3.4.2
# scipy            : 1.6.3
# seaborn          : 0.11.1
# pingouin	       : 0.3.12

# Import
from library_new_version import *

# 1. creat simplified and purified df for plotting
# +++ import data
file_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/1_Skurk/Datentabelle_3T3-L1 2011-12-20_2014-01-27_corri.xls'
sheet_name = '080609ens1 NGS Ens'
mouse_kic_df = pd.read_excel(file_path, sheet_name=sheet_name)  # open the correct sheet
mouse_kic_df.columns = mouse_kic_df.iloc[0]  # rename the columns
mouse_kic_df = mouse_kic_df.drop(index=0).reset_index(drop=True)  # delete the duplicate column names
mouse_kic_df = mouse_kic_df.iloc[630:824]
yaxislabel = 'Acylcarnitine products \n (normalized data on a log scale)'
xaxislabel = 'EGCG (µM)'
# +++ reorganise data
# 这里没有对名字做修改（其实是mouse egcg)
mouse_kic_df['Date'] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[0])
mouse_kic_df['Cell'] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[1])
mouse_kic_df['Time'] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[2])
mouse_kic_df['Repeat'] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[3])
mouse_kic_df['Charge'] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[4])
mouse_kic_df['Day'] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[5])
mouse_kic_df['Loading'] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[6])
mouse_kic_df['Treatment'] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[7])
mouse_kic_df[xaxislabel] = mouse_kic_df['Sample ID'].map(lambda x: x.split('_')[8])

# 2. Purify data
main_df = mouse_kic_df
main_df.replace(' ', '', regex=True, inplace=True)  # remove ' ' in all values
main_df.replace('µM', '', regex=True, inplace=True)  # remove 'µM' for better plotting
main_df.replace(',', '.', regex=True, inplace=True)
main_df.replace('K', 'H2O', regex=True, inplace=True)
Extract_int_df = main_df.loc[:, '101-C2:0': '229-C18:0-DC']  # extract int C0 to C18
Extract_int_df.rename(columns=lambda x: x[x.find('C'):len(x)], inplace=True)  # purify the column name
Extract_str_df = main_df[[xaxislabel, 'Charge', 'Treatment']].reset_index(drop=True)
# extract the real interested columns for plotting
for x in dict_LOD:
    Extract_int_df.loc[Extract_int_df[x] < dict_LOD[x], x] = NaN

# 3. normalise the data with 'CSA'
Extract_int_df_div = Extract_int_df.div(mouse_kic_df['CSA'], axis=0)  # for preparation of normalisation by CSA
x = Extract_int_df_div.values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
Extract_int_df_normalised = pd.DataFrame(x_scaled, columns=Extract_int_df_div.columns)  # give my col.name back
Simplified_df = pd.concat([Extract_str_df, Extract_int_df_normalised], axis=1)

# 4. Post-Purification
needed_df = Simplified_df[Simplified_df['Treatment'].isin(['DMSO', 'EGCG'])]
# !!! isin 和 Boolean判断，== 等等都可以，但是有多个判断条件的时候，使用isin真的很方便，而且也不用反向判断，使用drop函数
needed_oleic_yes = needed_df[needed_df['Charge'].isin(['Oleic+'])]
needed_oleic_yes.drop(['Charge', 'Treatment'], axis=1, inplace=True)

df = pd.melt(needed_oleic_yes, id_vars=[xaxislabel], var_name='C-Chains',
             value_name=yaxislabel)

# 5. Visualisation for lineplot
# 不能使用order的功能，只能在之前排序一次
df['C-Chains'] = df['C-Chains'].astype('category')  # critical step
df['C-Chains'].cat.reorder_categories(list_custom, inplace=True)
# inplace=True在原数据上生效
df_sorted = df.sort_values('C-Chains')
# 对sort函数的部分这里理解的其实还不够


# 6. Visualisation for errorbar problems. I will use groupby to get mean value(Success!)
# +++ calculate mean in set
grouped = df_sorted[yaxislabel].groupby([df_sorted[xaxislabel], df_sorted['C-Chains']]).mean()
# 'Normalised Value (mean)'是我传入的第一个数据，或者说我想要分析的数据。这里也可以不用写，因为就这一个剩下
# 如果想要一次传入多个数组，后面的就是多个column的引用作为key索引
# 如果直接groupby，他不会'智慧的'按照词条进行默认的索引分组

# +++ reset index
grouped_df = pd.DataFrame(grouped)
# 不需要使用这个把index变成列，下面的代码reset的时候，直接就分配到了列。alternative：grouped_df['Treatment'] = grouped_df.index
grouped_df.reset_index(inplace=True)

# 7. Plotting
dose_order = ['0', '0.1', '1', '6.25', '10', '12.5', '25', '30']
sns.set_theme( style='ticks', font='Helvetica',font_scale=1.5)
# dosis_ax = sns.pointplot(x='C-Chains', y=yaxislabel, data=df, hue=xaxislabel,
#                         palette='summer_r', hue_order=dose_order, order=list_custom,
#                        errwidth=0.2, alpha=.8,dodge=True, scale=.5,
#                        capsize=0.2, legend=True)
dosis_ax = sns.relplot(x='C-Chains', y=yaxislabel, data=df_sorted, hue=xaxislabel, kind='line', col=xaxislabel, col_wrap=2,
                       palette="viridis", col_order=dose_order, linewidth=4, zorder=5, hue_order=dose_order,
                       height=6, aspect=1.8, legend=False)
# !!!这种思维其实特别值得学习，先画一个图，然后用for loop在上面继续修饰

# +++ Iterate over each subplot to customize further
for concentration, ax in dosis_ax.axes_dict.items():  # 8个items （ax）

    # Add the title as an annotation within the plot 
    # 可惜这里添加了不好看，但这留下这个选择
    # ax.text(.8, .85, concentration, transform=ax.transAxes, fontweight="bold")

    # Plot every concentration's effect as series in the background
    sns.lineplot(
        data=grouped_df, x='C-Chains', y=yaxislabel, units=xaxislabel,
        estimator=None, color="C7", linewidth=1, ax=ax)
dosis_ax.set_xticklabels(labels=list_custom, rotation=90)
# alternative简单化x轴： ax.set_xticks(ax.get_xticks()[::2])


dosis_ax.set(yscale="log")
plt.tight_layout()
# save_path = '/Users/tianxingdu/Documents/4_Software Data/6_PycharmProjects/Bioinformatics/0_1_Work/1_Metabolitenassay/0_Output/20211025_26_7h/dosis_EGCG__all_lineplot_new.pdf'
# plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
plt.show()
