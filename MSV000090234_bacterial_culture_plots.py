#Receives mass spectrometry data of N-acyl amides from bacterial culture data,
#and visualizes through various plotting methods

#Data originally obtained from Douglas Guzior's MSV000090234 Massive Dataset,
#and queried for N-acyl amides with MassQL 
#MZMine3 was used to pre-process data obtained from GNPS2

#Code was originally run in Google Colaboratory 
#Running on other platforms may require reformatting heatmap dimensions

#Contact Prajit Rajkumar (prajkumar@ucsd.edu) for any questions

#Note that you may have to install some of the imports using 
#pip or another installer
import pandas as pd
import numpy as np

import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
import matplotlib.colors as mcolors
from matplotlib import rcParams
import scipy.spatial.distance as pdist

import annotationdata
from correlationmetrics import pearson
from queryretriever import prot2_retrieve
from queryretriever import gnps_retrieve

#Data Input

#Read in queries from Proteomics2 and GNPS2
prot2_input_data = pd.read_csv("/content/MSV000090234_Proteomics2_Data.csv")
gnps_input_data = pd.read_csv("/content/MSV000090234_GNPS2_Data.csv")

#Read in quantification table for MZMine3 processed data (from GNPS2)
quant_table = pd.read_csv("/content/MSV000090234_quantification_table.csv")
quant_table.drop(columns=(["row ion mobility",
                           "row ion mobility unit",
                           "row CCS",
                           "correlation group ID",
                           "annotation network number",
                           "best ion",
                           "auto MS2 verify",
                           "identified by n=",
                           "partners",
                           "neutral M mass"]), inplace=True)

#Download unprocessed files from Proteomics2

successful_heads = []
prot2_scan_num = []
for item in prot2_input_data.iterrows():
    prot2_scan_num.append(prot2_retrieve(item[1][1], "_" + item[1][0]))
    successful_heads.append(item[1][0])

#Change Proteomics2 data into two column table, remove empty hits
prot2_hits = pd.DataFrame({"Heads":successful_heads,"Scans":prot2_scan_num})
prot2_hits = prot2_hits[(prot2_hits != 0).all(1)]

#Extract culture condition and intesity data
prot2_column_headers = (["Heads",
                         "Culture_Conditions",
                         "Strain",
                         "Intensity",
                         "Delta_Mass"
                         ])
prot2_plot_data = pd.DataFrame(columns = prot2_column_headers)

i = 0
for row in prot2_hits.iterrows():
    data = pd.read_csv("/content/prot2_queried_file_" + row[1][0] + ".tsv",
                       delimiter="\t")
    for another_row in data.iterrows():
        intensity = another_row[1][7]
        delta_mass = round(another_row[1][3]
                           - annotationdata.head_masses[row[1][0]])
        if another_row[1][12].find("Sterile") != -1:
            culture_condition = (another_row[1][12])[:-7]
            strain_condition = culture_condition
        elif another_row[1][12].find("Blank") != -1:
            culture_condition = (another_row[1][12])[:-8]
            strain_condition = culture_condition
        else:
            culture_condition = (another_row[1][12])[-10:-7]
            strain_condition = ((another_row[1][12])[:-11]
                                + " "
                                + culture_condition
                                )
        prot2_plot_data.loc[i] = ([row[1][0],
                                   culture_condition,
                                   strain_condition,
                                   intensity,
                                   delta_mass]
                                  )
        i = i + 1

#Plot the intensities of unprocessed head group data from Proteomics2
#against culture conditon

rcParams['figure.figsize'] = 39,8.27
sns.stripplot(data=prot2_plot_data,
              x="Culture_Conditions",
              y="Intensity",
              hue="Heads",
              size=5,
              dodge=True,
              palette="Spectral"
              )

#Plot the intensities of unprocessed head group data from Proteomics2
#against strain data

rcParams['figure.figsize'] = 65,8.27
prot2_strains_plot = sns.stripplot(data=prot2_plot_data,
                                   x="Strain",
                                   y="Intensity",
                                   hue="Heads",
                                   size=5,
                                   dodge=True,
                                   palette="Spectral"
                                   )
prot2_strains_plot.tick_params(axis="x", rotation=90)

#Create a pivot table to plot a heatmap of head group data against strains of
#bacteria using unprocessed data from Proteomics2

prot2_strains_heads_pivot = pd.pivot_table(prot2_plot_data,
                                           values='Intensity',
                                           index=['Heads'],
                                           columns=['Strain'],
                                           aggfunc="sum"
                                           )
prot2_strains_heads_pivot.fillna(0, inplace=True)

#Create heatmap using previous pivot table

#Create the colormap
colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)]  #R -> G -> B
cmap_name = "white_blue_red"
cmap_wbr = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors)

#Plot data
rcParams['figure.figsize'] = 90,20
prot2_strains_heads_heatmap = sns.clustermap(prot2_strains_heads_pivot,
                                             norm = LogNorm(),
                                             cmap = cmap_wbr,
                                             method = "single",
                                             metric = pearson,
                                             col_cluster=False,
                                             row_cluster=True,
                                             xticklabels = 1,
                                             yticklabels = 1
                                             )

#Format heatmap
hm = prot2_strains_heads_heatmap.ax_heatmap.get_position()
prot2_strains_heads_heatmap.ax_heatmap.set_position([hm.x0,
                                                     hm.y0,
                                                     hm.width*3,
                                                     hm.height
                                                     ])

for _, spine in prot2_strains_heads_heatmap.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in prot2_strains_heads_heatmap.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

#Create a pivot table to plot a heatmap of head group data against acyl groups
#of bacteria using unprocessed data from Proteomics2

prot2_strains_acyls_pivot = pd.pivot_table(prot2_plot_data,
                                           values='Intensity',
                                           index=['Strain'],
                                           columns=['Delta_Mass'],
                                           aggfunc="sum"
                                           )
prot2_strains_acyls_pivot.fillna(0, inplace=True)

#Create heatmap using previous pivot table

#Plot heatmap
delta_masses_strains_plot = sns.clustermap(prot2_strains_acyls_pivot,
                                           norm = LogNorm(),
                                           cmap = cmap_wbr,
                                           method = "single",
                                           metric = pearson,
                                           col_cluster=False,
                                           row_cluster=True,
                                           xticklabels = 1,
                                           yticklabels = 1,
                                           dendrogram_ratio=(0.8,0.8)
                                           )

#Annotate and format the plot
hm = delta_masses_strains_plot.ax_heatmap.get_position()
delta_masses_strains_plot.ax_heatmap.set_position([hm.x0,
                                                   hm.y0,
                                                   2*hm.width,
                                                   4*hm.height
                                                   ])

masses = delta_masses_strains_plot.ax_heatmap.get_xticklabels()
acyl_groups = annotationdata.acyl_group_annotations

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(acyl_groups[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
delta_masses_strains_plot.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in delta_masses_strains_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in delta_masses_strains_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

row_dendro = delta_masses_strains_plot.ax_row_dendrogram.get_position()
delta_masses_strains_plot.ax_row_dendrogram.set_position([row_dendro.x0,
                                                          row_dendro.y0,
                                                          row_dendro.width,
                                                          4*row_dendro.height
                                                          ])

#Download processed files from GNPS2

gnps_successful_heads = []
gnps_scan_num = []
for row in gnps_input_data.iterrows():
    gnps_scan_num.append(gnps_retrieve(row[1][1], row[1][0]))
    gnps_successful_heads.append(row[1][0])

#Change GNPS2 data into formatted column table
gnps_column_headers = (["Head Group",
                        "Culture Condition",
                        "Strain",
                        "Peak Area",
                        "Delta Mass"
                        ])
gnps_plot_df = pd.DataFrame(columns = gnps_column_headers)
gnps_empty_files = []
j = 0

#Extract culture condition, strain, peak area, and acyl group data
for row in gnps_input_data.iterrows():
    with open("gnps_queried_file" + row[1][0] + ".tsv",) as input_file:
        if (input_file.read(15) == "<!DOCTYPE html>"):
            gnps_empty_files.append(row[1][0])
            continue
    temp_data = pd.read_csv("gnps_queried_file" + row[1][0] + ".tsv",
                            delimiter = "\t"
                            )
    for another_row in temp_data.iterrows():
        a_row = quant_table.loc[(quant_table["row ID"]
                                           == another_row[1][11]
                                           )]
        k = 3
        try:
            delta = round(a_row.iloc[0,1]
                          - annotationdata.head_masses[row[1][0]]
                          )
        except:
            continue
        for value in a_row.iloc[0][3:]:
            strain_cond = quant_table.iloc[:,k].name[:-17]
            if strain_cond.find("Blank") != -1:
                strain_cond = "Blank"
                culture_condit = "Blank"
            elif strain_cond.find("Sterile") != -1:
                strain_cond = "Sterile " + strain_cond[-3:]
                culture_condit = "Sterile " + strain_cond[-3:]
            else:
                culture_condit = strain_cond[-3:]
            gnps_plot_df.loc[j] = ([row[1][0],
                                    culture_condit,
                                    strain_cond,
                                    value,
                                    delta
                                    ])
            k = k + 1
            j = j + 1

gnps_plot_df.dropna(inplace=True)

#Plot the intensities of MZmine3 processed head group data from GNPS2
#against culture conditon

rcParams['figure.figsize'] = 39,8.27
sns.stripplot(data=gnps_plot_df,
              x="Culture Condition",
              y="Peak Area",
              hue="Head Group",
              size=5,
              dodge=True,
              order=(["CAT",
                      "RCM",
                      "Sterile CAT",
                      "Sterile RCM",
                      "Blank"
                      ])
              )

#Uncomment to plot log of the y-axis
#plt.gca().set_yscale('log')

#Divide data by media

#CAT data
sterile_RCM = gnps_plot_df[gnps_plot_df['Culture Condition']
                           == 'Sterile RCM'].index
RCM = gnps_plot_df[gnps_plot_df['Culture Condition'] == 'RCM'].index
all_RCM = RCM.union(sterile_RCM)

gnps_plot_df_cat = gnps_plot_df.drop(all_RCM)

#RCM data
sterile_CAT = gnps_plot_df[gnps_plot_df['Culture Condition']
                           == 'Sterile CAT'].index
CAT = gnps_plot_df[gnps_plot_df['Culture Condition'] == 'CAT'].index
all_CAT = CAT.union(sterile_CAT)

gnps_plot_df_rcm = gnps_plot_df.drop(all_CAT)

#Create and normalize pivot tables to plot CAT and RCM data of heads associated
#with different strains of bacteria from the MZmine3 processed data from GNPS2

#CAT media
gnps2_strains_heads_pivot_CAT = pd.pivot_table(gnps_plot_df_cat,
                                               values='Peak Area',
                                               index=['Head Group'],
                                               columns=['Strain'],
                                               aggfunc="mean"
                                               )

gnps2_strains_heads_pivot_CAT.fillna(0, inplace=True)
for column in gnps2_strains_heads_pivot_CAT:
    if gnps2_strains_heads_pivot_CAT[column].name == "Blank":
        pass
    else:
        gnps2_strains_heads_pivot_CAT[column] = (gnps2_strains_heads_pivot_CAT[column]
                                                 - gnps2_strains_heads_pivot_CAT["Sterile CAT"]
                                                 )
gnps2_strains_heads_pivot_CAT.drop(columns=["Sterile CAT"], inplace=True)

#RCM media
gnps2_strains_heads_pivot_RCM = pd.pivot_table(gnps_plot_df_rcm,
                                           values='Peak Area',
                                           index=['Head Group'],
                                           columns=['Strain'],
                                           aggfunc="mean"
                                           )
gnps2_strains_heads_pivot_RCM.fillna(0, inplace=True)
for column in gnps2_strains_heads_pivot_RCM:
    if gnps2_strains_heads_pivot_RCM[column].name == "Blank":
        pass
    else:
        gnps2_strains_heads_pivot_RCM[column] = (gnps2_strains_heads_pivot_RCM[column]
                                                 - gnps2_strains_heads_pivot_RCM["Sterile RCM"]
                                                 )
gnps2_strains_heads_pivot_RCM.drop(columns=["Sterile RCM"], inplace=True)

#Strains vs Heads CAT

#Plot data
rcParams['figure.figsize'] = 90,20
gnps_strains_heads_cat = sns.clustermap(gnps2_strains_heads_pivot_CAT,
                                        norm = SymLogNorm(linthresh=10),
                                        cmap = "vlag",
                                        method = "single",
                                        metric = pearson,
                                        col_cluster=False,
                                        row_cluster=True,
                                        xticklabels = 1,
                                        yticklabels = 1
                                        )

#Format heatmap
hm = gnps_strains_heads_cat.ax_heatmap.get_position()
gnps_strains_heads_cat.ax_heatmap.set_position([hm.x0,
                                                hm.y0,
                                                hm.width*3,
                                                hm.height
                                                ])

for _, spine in gnps_strains_heads_cat.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in gnps_strains_heads_cat.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

#Strains vs Heads RCM

#Plot data
rcParams['figure.figsize'] = 90,20
gnps_strains_heads_rcm = sns.clustermap(gnps2_strains_heads_pivot_RCM,
                                        norm = SymLogNorm(linthresh=10,
                                                          vmin=-10e6,
                                                          vmax = 10e6
                                                          ),
                                        cmap = "vlag",
                                        method = "single",
                                        metric = pearson,
                                        col_cluster=False,
                                        row_cluster=True,
                                        xticklabels = 1,
                                        yticklabels = 1
                                        )

#Format heatmap
hm = gnps_strains_heads_rcm.ax_heatmap.get_position()
gnps_strains_heads_rcm.ax_heatmap.set_position([hm.x0,
                                                hm.y0,
                                                hm.width*3,
                                                hm.height
                                                ])

for _, spine in gnps_strains_heads_rcm.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in gnps_strains_heads_rcm.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

#Create and normalize pivot tables to plot CAT and RCM data of acyl groups
#associated with different strains of bacteria from the MZmine3 processed
#data from GNPS2

#CAT media
gnps2_strains_acyls_pivot_CAT = pd.pivot_table(gnps_plot_df_cat,
                                               values='Peak Area',
                                               index=['Delta Mass'],
                                               columns=['Strain'],
                                               aggfunc="mean"
                                               )
gnps2_strains_acyls_pivot_CAT.fillna(0, inplace=True)
for column in gnps2_strains_acyls_pivot_CAT:
    if gnps2_strains_acyls_pivot_CAT[column].name == "Blank":
        pass
    else:
        gnps2_strains_acyls_pivot_CAT[column] = (gnps2_strains_acyls_pivot_CAT[column]
                                                 -
                                                 gnps2_strains_acyls_pivot_CAT["Sterile CAT"]
                                                 )
gnps2_strains_acyls_pivot_CAT.drop(columns=["Sterile CAT"], inplace=True)
gnps2_strains_acyls_pivot_CAT = gnps2_strains_acyls_pivot_CAT.transpose()

#RCM media
gnps2_strains_acyls_pivot_RCM = pd.pivot_table(gnps_plot_df_rcm,
                                               values='Peak Area',
                                               index=['Delta Mass'],
                                               columns=['Strain'],
                                               aggfunc="mean"
                                               )
gnps2_strains_acyls_pivot_RCM.fillna(0, inplace=True)
for column in gnps2_strains_acyls_pivot_RCM:
    if gnps2_strains_acyls_pivot_RCM[column].name == "Blank":
        pass
    else:
        gnps2_strains_acyls_pivot_RCM[column] = (gnps2_strains_acyls_pivot_RCM[column]
                                                 -
                                                 gnps2_strains_acyls_pivot_RCM["Sterile RCM"]
                                                 )
gnps2_strains_acyls_pivot_RCM.drop(columns=["Sterile RCM"], inplace=True)
gnps2_strains_acyls_pivot_RCM = gnps2_strains_acyls_pivot_RCM.transpose()

#Strains vs Acyl Groups CAT

#Plot heatmap
delta_masses_strains_plot = sns.clustermap(gnps2_strains_acyls_pivot_CAT,
                                           norm = SymLogNorm(linthresh=10,
                                                             vmin=-10e6,
                                                             vmax=10e6
                                                             ),
                                           cmap = "vlag",
                                           method = "single",
                                           metric = pearson,
                                           col_cluster=False,
                                           row_cluster=True,
                                           xticklabels = 1,
                                           yticklabels = 1,
                                           dendrogram_ratio=(0.8,0.8)
                                           )

#Annotate and format the plot
hm = delta_masses_strains_plot.ax_heatmap.get_position()
delta_masses_strains_plot.ax_heatmap.set_position([hm.x0,
                                                   hm.y0,
                                                   2*hm.width,
                                                   4*hm.height
                                                   ])

masses = delta_masses_strains_plot.ax_heatmap.get_xticklabels()

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(annotationdata.acyl_group_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
delta_masses_strains_plot.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in delta_masses_strains_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in delta_masses_strains_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

row_dendro = delta_masses_strains_plot.ax_row_dendrogram.get_position()
delta_masses_strains_plot.ax_row_dendrogram.set_position([row_dendro.x0,
                                                          row_dendro.y0,
                                                          row_dendro.width,
                                                          4*row_dendro.height
                                                          ])

#Strains vs Acyl Groups RCM

#Plot heatmap
delta_masses_strains_plot = sns.clustermap(gnps2_strains_acyls_pivot_RCM,
                                           norm = SymLogNorm(linthresh=10,
                                                             vmin=-10e6,
                                                             vmax=10e6
                                                             ),
                                           cmap = "vlag",
                                           method = "single",
                                           metric = pearson,
                                           col_cluster=False,
                                           row_cluster=True,
                                           xticklabels = 1,
                                           yticklabels = 1,
                                           dendrogram_ratio=(0.8,0.8)
                                           )

#Annotate and format the plot
hm = delta_masses_strains_plot.ax_heatmap.get_position()
delta_masses_strains_plot.ax_heatmap.set_position([hm.x0,
                                                   hm.y0,
                                                   2*hm.width,
                                                   4*hm.height
                                                   ])

masses = delta_masses_strains_plot.ax_heatmap.get_xticklabels()

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(annotationdata.acyl_group_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
delta_masses_strains_plot.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in delta_masses_strains_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in delta_masses_strains_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

row_dendro = delta_masses_strains_plot.ax_row_dendrogram.get_position()
delta_masses_strains_plot.ax_row_dendrogram.set_position([row_dendro.x0,
                                                          row_dendro.y0,
                                                          row_dendro.width,
                                                          4*row_dendro.height
                                                          ])

#Create and normalize pivot tables to plot CAT and RCM data of acyl groups
#associated with different head groups from the MZmine3 processed
#data from GNPS2

cat_steriles = gnps_plot_df_cat[gnps_plot_df_cat['Culture Condition']
                                == 'Sterile CAT'].index
cat_blanks = gnps_plot_df_cat[gnps_plot_df_cat['Culture Condition']
                              == 'Blank'].index
cat_drops = cat_steriles.union(cat_blanks)
gnps_plot_df_cat_steriles = gnps_plot_df.loc[cat_steriles]
gnps_plot_df_cat_cultures = gnps_plot_df_cat.drop(cat_drops)


gnps2_heads_acyls_pivot_CAT_cultures = pd.pivot_table(gnps_plot_df_cat_cultures,
                                                      values='Peak Area',
                                                      index=['Head Group'],
                                                      columns=['Delta Mass'],
                                                      aggfunc="mean"
                                                      )

gnps2_heads_acyls_pivot_CAT_steriles = pd.pivot_table(gnps_plot_df_cat_steriles,
                                                      values='Peak Area',
                                                      index=['Head Group'],
                                                      columns=['Delta Mass'],
                                                      aggfunc="mean"
                                                      )

gnps2_heads_acyls_pivot_CAT_cultures.fillna(0, inplace=True)
gnps2_heads_acyls_pivot_CAT_steriles.fillna(0, inplace=True)
gnps2_heads_acyls_pivot_CAT = gnps2_heads_acyls_pivot_CAT_cultures.subtract(gnps2_heads_acyls_pivot_CAT_steriles)



rcm_steriles = gnps_plot_df_rcm[gnps_plot_df_rcm['Culture Condition']
                                == 'Sterile RCM'].index
rcm_blanks = gnps_plot_df_rcm[gnps_plot_df_rcm['Culture Condition']
                              == 'Blank'].index
rcm_drops = rcm_steriles.union(rcm_blanks)
gnps_plot_df_rcm_steriles = gnps_plot_df.loc[rcm_steriles]
gnps_plot_df_rcm_cultures = gnps_plot_df_rcm.drop(rcm_drops)


gnps2_heads_acyls_pivot_RCM_cultures = pd.pivot_table(gnps_plot_df_rcm_cultures,
                                                      values='Peak Area',
                                                      index=['Head Group'],
                                                      columns=['Delta Mass'],
                                                      aggfunc="mean"
                                                      )

gnps2_heads_acyls_pivot_RCM_steriles = pd.pivot_table(gnps_plot_df_rcm_steriles,
                                                      values='Peak Area',
                                                      index=['Head Group'],
                                                      columns=['Delta Mass'],
                                                      aggfunc="mean"
                                                      )

gnps2_heads_acyls_pivot_RCM_cultures.fillna(0, inplace=True)
gnps2_heads_acyls_pivot_RCM_steriles.fillna(0, inplace=True)
gnps2_heads_acyls_pivot_RCM = gnps2_heads_acyls_pivot_RCM_cultures.subtract(gnps2_heads_acyls_pivot_RCM_steriles)

#Acyl Groups vs Heads CAT

#Plot heatmap
delta_masses_strains_plot = sns.clustermap(gnps2_heads_acyls_pivot_CAT,
                                           norm = SymLogNorm(linthresh=10,
                                                             vmin=-10e6,
                                                             vmax=10e6
                                                             ),
                                           cmap = "vlag",
                                           method = "single",
                                           metric = pearson,
                                           col_cluster=False,
                                           row_cluster=True,
                                           xticklabels = 1,
                                           yticklabels = 1,
                                           dendrogram_ratio=(0.8,
                                                             0.8
                                                             )
                                           )

#Annotate and format data
hm = delta_masses_strains_plot.ax_heatmap.get_position()
delta_masses_strains_plot.ax_heatmap.set_position([hm.x0,
                                                   hm.y0,
                                                   2*hm.width,
                                                   4*hm.height
                                                   ])

#Annotate and format the plot
masses = delta_masses_strains_plot.ax_heatmap.get_xticklabels()

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(annotationdata.acyl_group_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
delta_masses_strains_plot.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in delta_masses_strains_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in delta_masses_strains_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

row_dendro = delta_masses_strains_plot.ax_row_dendrogram.get_position()
delta_masses_strains_plot.ax_row_dendrogram.set_position([row_dendro.x0,
                                                          row_dendro.y0,
                                                          row_dendro.width,
                                                          4*row_dendro.height
                                                          ])

#Acyl Groups vs Heads RCM

#Plot heatmap
delta_masses_strains_plot = sns.clustermap(gnps2_heads_acyls_pivot_RCM,
                                           norm = SymLogNorm(linthresh=10,
                                                             vmin=-10e6,
                                                             vmax = 10e6),
                                           cmap = "vlag",
                                           method = "single",
                                           metric = pearson,
                                           col_cluster=False,
                                           row_cluster=True,
                                           xticklabels = 1,
                                           yticklabels = 1,
                                           dendrogram_ratio=(0.8,
                                                             0.8
                                                             )
                                           )

hm = delta_masses_strains_plot.ax_heatmap.get_position()
delta_masses_strains_plot.ax_heatmap.set_position([hm.x0,
                                                   hm.y0,
                                                   2*hm.width,
                                                   4*hm.height
                                                   ])

#Annotate and format the plot
masses = delta_masses_strains_plot.ax_heatmap.get_xticklabels()

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(annotationdata.acyl_group_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
delta_masses_strains_plot.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in delta_masses_strains_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in delta_masses_strains_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

row_dendro = delta_masses_strains_plot.ax_row_dendrogram.get_position()
delta_masses_strains_plot.ax_row_dendrogram.set_position([row_dendro.x0,
                                                          row_dendro.y0,
                                                          row_dendro.width,
                                                          4*row_dendro.height
                                                          ])

#Create normalized dataframes to plot heads vs. deltas for B. fragilis data

#CAT
fragilis_CAT_unprocessed = gnps_plot_df_cat_cultures[gnps_plot_df_cat_cultures['Strain']
                                                     == 'B-fragilis-CAT']

fragilis_CAT_unprocessed_pivot = pd.pivot_table(fragilis_CAT_unprocessed,
                                                values='Peak Area',
                                                index=['Head Group'],
                                                columns=['Delta Mass'],
                                                aggfunc="mean"
                                                )
fragilis_CAT_unprocessed_pivot.fillna(0, inplace=True)

fragilis_CAT_normalized = fragilis_CAT_unprocessed_pivot.subtract(gnps2_heads_acyls_pivot_CAT_steriles)

#RCM
fragilis_RCM_unprocessed = gnps_plot_df_rcm_cultures[gnps_plot_df_rcm_cultures['Strain']
                                                     == 'B-fragilis-RCM']

fragilis_RCM_unprocessed_pivot = pd.pivot_table(fragilis_RCM_unprocessed,
                                                values='Peak Area',
                                                index=['Head Group'],
                                                columns=['Delta Mass'],
                                                aggfunc="mean"
                                                )
fragilis_RCM_unprocessed_pivot.fillna(0, inplace=True)

fragilis_RCM_normalized = fragilis_RCM_unprocessed_pivot.subtract(gnps2_heads_acyls_pivot_RCM_steriles)

#Acyl Groups vs Heads CAT

#Plot heatmap
#Clustering is not applied as most rows have nearly identical values
fragilis_CAT_deltas_strains_plt = sns.clustermap(fragilis_CAT_normalized,
                                                 norm = SymLogNorm(linthresh=10,
                                                                   vmin=-10e6,
                                                                   vmax=10e6
                                                                   ),
                                                 cmap = "vlag",
                                                 method = "single",
                                                 metric = pearson,
                                                 col_cluster=False,
                                                 row_cluster=False,
                                                 xticklabels = 1,
                                                 yticklabels = 1,
                                                 dendrogram_ratio=(0.8,
                                                                   0.8
                                                                   )
                                                 )

#Annotate and format data
hm = fragilis_CAT_deltas_strains_plt.ax_heatmap.get_position()
fragilis_CAT_deltas_strains_plt.ax_heatmap.set_position([hm.x0,
                                                         hm.y0,
                                                         2*hm.width,
                                                         4*hm.height
                                                         ])

#Annotate and format the plot
masses = fragilis_CAT_deltas_strains_plt.ax_heatmap.get_xticklabels()

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(annotationdata.acyl_group_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
fragilis_CAT_deltas_strains_plt.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in fragilis_CAT_deltas_strains_plt.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in fragilis_CAT_deltas_strains_plt.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

row_dendro = fragilis_CAT_deltas_strains_plt.ax_row_dendrogram.get_position()
fragilis_CAT_deltas_strains_plt.ax_row_dendrogram.set_position([row_dendro.x0,
                                                                row_dendro.y0,
                                                                row_dendro.width,
                                                                4*row_dendro.height
                                                                ])

#Acyl Groups vs Heads RCM

#Plot heatmap
#Clustering is not applied to make comparison easier
fragilis_RCM_deltas_strains_plt = sns.clustermap(fragilis_RCM_normalized,
                                                 norm = SymLogNorm(linthresh=10,
                                                                   vmin=-10e6,
                                                                   vmax=10e6
                                                                   ),
                                                 cmap = "vlag",
                                                 method = "single",
                                                 metric = pearson,
                                                 col_cluster=False,
                                                 row_cluster=False,
                                                 xticklabels = 1,
                                                 yticklabels = 1,
                                                 dendrogram_ratio=(0.8,
                                                                   0.8
                                                                   )
                                                 )

#Annotate and format data
hm = fragilis_RCM_deltas_strains_plt.ax_heatmap.get_position()
fragilis_RCM_deltas_strains_plt.ax_heatmap.set_position([hm.x0,
                                                         hm.y0,
                                                         2*hm.width,
                                                         4*hm.height
                                                         ])

#Annotate and format the plot
masses = fragilis_RCM_deltas_strains_plt.ax_heatmap.get_xticklabels()

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(annotationdata.acyl_group_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
fragilis_RCM_deltas_strains_plt.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in fragilis_RCM_deltas_strains_plt.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in fragilis_RCM_deltas_strains_plt.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

row_dendro = fragilis_RCM_deltas_strains_plt.ax_row_dendrogram.get_position()
fragilis_RCM_deltas_strains_plt.ax_row_dendrogram.set_position([row_dendro.x0,
                                                                row_dendro.y0,
                                                                row_dendro.width,
                                                                4*row_dendro.height
                                                                ])