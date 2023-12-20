#Receives mass spectrometry data of N-acyl amides from bacterial culture data,  
#and visualizes through various plotting methods

#Data originally obtained from Emily Gentry's MSV000084475 Massive Dataset, 
#RCM Media, and queried for N-acyl amides with MassQL 
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
from matplotlib.colors import LogNorm, Normalize
import matplotlib.colors as mcolors
from matplotlib import rcParams
import scipy.spatial.distance as pdist

import annotationdata
import queryretriever
from correlationmetrics import pearson

#Data Input

#Read in queries from Proteomics2 and GNPS2
prot2_input_data = pd.read_csv("/content/MSV000084475_RCM_Proteomics2_Data.csv")

gnps_input_data = pd.read_csv("/content/MSV000084475_RCM_GNPS2_Data.csv")

#Read in quantification table for MZMine3 processed data (from GNPS2)
quant_table = pd.read_csv("/content/MSV000084475_RCM_quantification_table.csv")

#Read in metadata
culture_metadata = ({row[0] : row[1] for _,
                     row in
                     pd.read_csv("/content/MSV000084475"
                                 + "_RCM_culture_metadata.csv"
                                 ).iterrows()
                     })

species = ({row[0] : row[11] for _,
            row in pd.read_csv("/content/MSV000084475_strain_metadata.tsv",
                               delimiter = "\t").iterrows()})
genus = ({row[0] : row[10] for _,
         row in pd.read_csv("/content/MSV000084475_strain_metadata.tsv",
                            delimiter = "\t").iterrows()})
culture_cond = ({row[0] : row[16] for _,
                 row in pd.read_csv("/content/MSV000084475_strain_metadata.tsv",
                                    delimiter = "\t").iterrows()})

#Process quantification table to format for data display

#Contains MZMine3 processed files without a label
gnps_unlabeled_files = []

#Create processed quantification table
processed_quant_table = pd.DataFrame()
processed_quant_table["row ID"] = quant_table.iloc[:,0]
processed_quant_table["row m/z"] = quant_table.iloc[:,1]
processed_quant_table["row retention time"] = quant_table.iloc[:,2]

#Loop through original quantification table and correctly format columns
i = 1000
for column in quant_table.iloc[:,3:]:
    try:
        if (species[quant_table[column].name[:-10]] != "not applicable"):
            taxon = species[quant_table[column].name[:-10]].replace(' ', '\\ ')
            strain = ("$\it{"
                      + taxon
                      + "}$ "
                      + culture_metadata[quant_table[column].name[:-10]]
                      )
        elif (genus[quant_table[column].name[:-10]] != "not applicable"):
            taxon = genus[quant_table[column].name[:-10]].replace(' ', '\\ ')
            strain = ("$\it{"
                      + taxon
                      + "}$ "
                      + culture_metadata[quant_table[column].name[:-10]]
                      )
        elif (culture_cond[quant_table[column].name[:-10]] == "blank_media"):
            strain = "media blank"
        elif (culture_cond[quant_table[column].name[:-10]]
              == "blank_extraction"):
            strain = "extraction blank"
        else:
            strain = culture_metadata[quant_table[column].name[:-10]]
        strain = strain.replace("culture ", "")
        quant_table[column].name = strain + str(i)
        new_column = pd.DataFrame(quant_table[column])
        processed_quant_table = pd.concat([processed_quant_table,
                                           new_column],
                                          axis = 1
                                          )
        i = i + 1
    except KeyError:
        gnps_unlabeled_files.append(quant_table[column].name[:-10])
        quant_table.drop(quant_table[column].name, axis = 1)

#Download data from Proteomics2
prot2_successful_heads = []
prot2_scan_num = []
for item in prot2_input_data.iterrows():
    prot2_scan_num.append(queryretriever.prot2_retrieve(item[1][1],
                                                        "_" + item[1][0]))
    prot2_successful_heads.append(item[1][0])

#Change Proteomics2 data into two column table, remove empty hits
prot2_hits = pd.DataFrame({"Heads":prot2_successful_heads,
                           "Scans":prot2_scan_num})
prot2_hits = prot2_hits[(prot2_hits != 0).all(1)]

#Extract culture condition, strain, intensity, and acyl group data
prot2_column_headers = (["Head Group",
                         "Culture Condition",
                         "Strain",
                         "Intensity",
                         "Delta Mass"
                         ])
prot2_plot_data = pd.DataFrame(columns = prot2_column_headers)

i = 0
for row in prot2_hits.iterrows():
    data = pd.read_csv("/content/prot2_queried_file_" + row[1][0] + ".tsv",
                       delimiter="\t")
    for another_row in data.iterrows():
        intensity = another_row[1][7]
        delta_mass = round(another_row[1][3]
                           - annotationdata.head_masses[row[1][0]]
                           )
        culture_condition = culture_cond[(another_row[1][12])[:-2] + "XML"]
        if culture_condition == "blank_media":
            culture_condition = "media blank"
            strain_cond = culture_condition
        elif culture_condition == "blank_extraction":
            culture_condition = "extraction blank"
            strain_cond = culture_condition
        else:
            culture_condition = culture_metadata[(another_row[1][12])[:-2]
                                                 + "XML"]
            if (species[(another_row[1][12])[:-2] + "XML"] != "not applicable"):
                strain_cond = species[(another_row[1][12])[:-2] + "XML"]
                strain_cond = "$\it{" + strain_cond.replace(' ', '\\ ') + "}$"
            elif (genus[(another_row[1][12])[:-2] + "XML"] != "not applicable"):
                strain_cond = genus[(another_row[1][12])[:-2] + "XML"]
                strain_cond = "$\it{" + strain_cond.replace(' ', '\\ ') + "}$"
            else:
                strain_cond = "Unavailable"
        prot2_plot_data.loc[i] = ([row[1][0],
                                  culture_condition,
                                  strain_cond,
                                  intensity,
                                  delta_mass
                                  ])
        i = i + 1

#Create a scatterplot of the head groups by culture condition for
#Proteomics2 Data
rcParams['figure.figsize'] = 39,8.27
plot_palette = sns.color_palette("Spectral", 30)
sns.stripplot(data=prot2_plot_data,
              x="Culture Condition",
              y="Intensity",
              hue="Head Group",
              size=5,
              dodge=True,
              order = (['72h',
                        'media blank',
                        'extraction blank'
                        ])
              )

#Reformat the Proteomics2 head group and strain data into a pivot table
prot2_strains_heads_pivot = pd.pivot_table(prot2_plot_data,
                                           values='Intensity',
                                           index=['Head Group'],
                                           columns=['Strain'],
                                           aggfunc="sum"
                                           )
prot2_strains_heads_pivot.fillna(0, inplace=True)

#Plot a heatmap of the head group and strain data for Proteomics2 Data

#Create the colormap
colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)] #R -> G -> B
cmap_name = "white_blue_red"
cmap_wbr = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors)

#Plot the heatmap
rcParams['figure.figsize'] = 90,20
strains_heads_plot = sns.clustermap(prot2_strains_heads_pivot,
                                    norm = LogNorm(),
                                    cmap = cmap_wbr,
                                    method = "single",
                                    metric = pearson,
                                    col_cluster=False,
                                    row_cluster=True,
                                    xticklabels = 1,
                                    yticklabels = 1
                                    )

#Format the plot
hm = strains_heads_plot.ax_heatmap.get_position()
strains_heads_plot.ax_heatmap.set_position([hm.x0,
                                            hm.y0,
                                            hm.width*3,
                                            hm.height
                                            ])

for _, spine in strains_heads_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in strains_heads_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

#Reformat the Proteomics2 delta masses and strain data into a pivot table
prot2_strains_acyls_pivot = pd.pivot_table(prot2_plot_data,
                                           values='Intensity',
                                           index=['Strain'],
                                           columns=['Delta Mass'],
                                           aggfunc="sum"
                                           )
prot2_strains_acyls_pivot.fillna(0, inplace=True)

#Plot a heatmap of the delta masses and strain data for the unprocessed
#Proteomics2 Data

#Create the colormap
colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)]  #R -> G -> B
cmap_name = "white_blue_red"
cmap_wbr = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors)

deltas_strains_plot = sns.clustermap(prot2_strains_acyls_pivot,
                                     norm = LogNorm(),
                                     cmap = cmap_wbr,
                                     method = "single",
                                     metric = pearson,
                                     col_cluster=False,
                                     row_cluster=True,
                                     xticklabels = 1,
                                     yticklabels = 1
                                     )

hm = deltas_strains_plot.ax_heatmap.get_position()
deltas_strains_plot.ax_heatmap.set_position([hm.x0,
                                             hm.y0,
                                             hm.width,
                                             hm.height
                                             ])

#Annotate and format the plot
masses = deltas_strains_plot.ax_heatmap.get_xticklabels()
acyl_annotations = annotationdata.acyl_group_annotations

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(acyl_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
deltas_strains_plot.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in deltas_strains_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in deltas_strains_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

#Download processed files from GNPS2

gnps_successful_heads = []
gnps_scan_num = []
for row in gnps_input_data.iterrows():
    gnps_scan_num.append(queryretriever.gnps_retrieve(row[1][1], row[1][0]))
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
        a_row = processed_quant_table.loc[(processed_quant_table["row ID"]
                                           == another_row[1][11]
                                           )]
        k = 3
        delta = round(a_row.iloc[0,1] - annotationdata.head_masses[row[1][0]])
        for value in a_row.iloc[0][3:]:
            strain_cond = processed_quant_table.iloc[:,k].name[:-4]
            culture_condit = strain_cond
            if (culture_condit != 'extraction blank'
                and culture_condit != 'media blank'
                and culture_condit != '72h'
                and culture_condit != '0h'):
                culture_condit = culture_condit[culture_condit.find("}$") + 2:]
            gnps_plot_df.loc[j] = ([row[1][0],
                                    culture_condit,
                                    strain_cond,
                                    value,
                                    delta
                                    ])
            k = k + 1
            j = j + 1

gnps_plot_df['Culture Condition'] = (gnps_plot_df['Culture Condition'].
                                     replace(['72h'], ' 72h')
                                     )

#Create a scatterplot of the head groups by culture condition for
#GNPS2 Data
rcParams['figure.figsize'] = 39,8.27
sns.stripplot(data = gnps_plot_df,
              x="Culture Condition",
              y="Peak Area",
              hue="Head Group",
              size=5,
              dodge=True,
              order=([' 72h',
                      ' 0h',
                      'media blank',
                      'extraction blank'
                      ])
              )
plt.gca().set_yscale('log')

#Prepare the GNPS2 data to allow for normalization relative to time 0 hours
#to focus on production of metabolites

#Define function to format new columns
def aggregation_formatter(strain_condition):
    if (strain_condition != 'extraction blank'
        and strain_condition != 'media blank'
        and strain_condition.find("}$") != -1):
        return strain_condition[:strain_condition.find("}$") + 2]
    else:
        return strain_condition

#Create dataframe of cultures taken at 72h
agg_gnps_df_72h = gnps_plot_df.copy(deep=True)
start_values = agg_gnps_df_72h[agg_gnps_df_72h['Culture Condition']
                               == ' 0h'].index
agg_gnps_df_72h.drop(start_values , inplace=True)
agg_gnps_df_72h['Strain'] = agg_gnps_df_72h['Strain'].apply(aggregation_formatter)

#Create dataframe of cultures taken at 0h
agg_gnps_df_0h = gnps_plot_df.copy(deep=True)
end_values = agg_gnps_df_0h[agg_gnps_df_0h['Culture Condition'] == ' 72h'].index
agg_gnps_df_0h.drop(end_values , inplace=True)
agg_gnps_df_0h.loc[agg_gnps_df_0h['Culture Condition'] == 'extraction blank',
                   'Peak Area'
                   ] = 0
agg_gnps_df_0h.loc[agg_gnps_df_0h['Culture Condition'] == 'media blank',
                   'Peak Area'
                   ] = 0
agg_gnps_df_0h['Strain'] = agg_gnps_df_0h['Strain'].apply(aggregation_formatter)

#Reformat the GNPS2 head group and strain data into a pivot table for 72h and
#0h, and subtract the two to normalize and prepare for plotting

#Pivot table for 72h
gnps_strains_heads_pivot_72h = pd.pivot_table(agg_gnps_df_72h,
                                           values='Peak Area',
                                           index=['Head Group'],
                                           columns=['Strain'],
                                           aggfunc="sum"
                                           )
gnps_strains_heads_pivot_72h.fillna(0, inplace=True)
gnps_strains_heads_pivot_72h.drop(columns=['72h'], inplace=True)

#Pivot table for 0h
gnps_strains_heads_pivot_0h = pd.pivot_table(agg_gnps_df_0h,
                                           values='Peak Area',
                                           index=['Head Group'],
                                           columns=['Strain'],
                                           aggfunc="sum"
                                           )
gnps_strains_heads_pivot_0h.fillna(0, inplace=True)

gnps_strains_heads_pivot_plot = (gnps_strains_heads_pivot_72h
                                 .sub(gnps_strains_heads_pivot_0h)
                                 )

#Plot a heatmap of the head group and strain data for MZmine3 processed
#GNPS2 Data

#Create the colormap
colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)] #R -> G -> B
cmap_name = "white_blue_red"
cmap_wbr = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors)

#Plot the heatmap
rcParams['figure.figsize'] = 90,20
strains_heads_plot = sns.clustermap(gnps_strains_heads_pivot_plot,
                                    norm = LogNorm(),
                                    cmap = cmap_wbr,
                                    method = "single",
                                    metric = pearson,
                                    col_cluster=False,
                                    row_cluster=True,
                                    xticklabels = 1,
                                    yticklabels = 1
                                    )

#Format the plot
hm = strains_heads_plot.ax_heatmap.get_position()
strains_heads_plot.ax_heatmap.set_position([hm.x0,
                                            hm.y0,
                                            hm.width*3,
                                            hm.height
                                            ])

for _, spine in strains_heads_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in strains_heads_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

#Reformat the GNPS2 delta mass and strain data into a pivot table for 72h and
#0h, and subtract the two to normalize and prepare for plotting

#Pivot table for 72h
gnps_strains_deltas_pivot_72h = pd.pivot_table(agg_gnps_df_72h,
                                           values='Peak Area',
                                           index=['Strain'],
                                           columns=['Delta Mass'],
                                           aggfunc="sum"
                                           )
gnps_strains_deltas_pivot_72h.fillna(0, inplace=True)
gnps_strains_deltas_pivot_72h.drop(['72h'], inplace=True)

#Pivot table for 0h
gnps_strains_deltas_pivot_0h = pd.pivot_table(agg_gnps_df_0h,
                                           values='Peak Area',
                                           index=['Strain'],
                                           columns=['Delta Mass'],
                                           aggfunc="sum"
                                           )
gnps_strains_deltas_pivot_0h.fillna(0, inplace=True)

gnps_strains_deltas_pivot_plot = (gnps_strains_deltas_pivot_72h
                                 .sub(gnps_strains_deltas_pivot_0h)
                                 )

#Plot a heatmap of the delta masses and strain data for the MZmine3 processed
#GNPS2 Data

#Create the colormap
colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)]  #R -> G -> B
cmap_name = "white_blue_red"
cmap_wbr = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors)

#Plot the data
gn_deltas_strains_plot = sns.clustermap(gnps_strains_deltas_pivot_plot,
                                        norm = LogNorm(),
                                        cmap = cmap_wbr,
                                        method = "single",
                                        metric = pearson,
                                        col_cluster=False,
                                        row_cluster=True,
                                        xticklabels = 1,
                                        yticklabels = 1
                                       )

hm = gn_deltas_strains_plot.ax_heatmap.get_position()
gn_deltas_strains_plot.ax_heatmap.set_position([hm.x0,
                                                hm.y0,
                                                hm.width,
                                                1.3*hm.height
                                                ])

#Annotate and format the plot
masses = gn_deltas_strains_plot.ax_heatmap.get_xticklabels()
acyl_annotations = annotationdata.acyl_group_annotations

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(acyl_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
gn_deltas_strains_plot.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in gn_deltas_strains_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in gn_deltas_strains_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

#Reformat the GNPS2 heads and deltas data into a pivot table for 72h and
#0h, and subtract the two to normalize and prepare for plotting

#Pivot table for 72h
gnps_head_delta_pivot_72h = pd.pivot_table(agg_gnps_df_72h,
                                           values='Peak Area',
                                           index=['Head Group'],
                                           columns=['Delta Mass'],
                                           aggfunc="sum"
                                           )
gnps_head_delta_pivot_72h.fillna(0, inplace=True)

#Pivot table for 0h
gnps_head_delta_pivot_0h = pd.pivot_table(agg_gnps_df_0h,
                                          values='Peak Area',
                                          index=['Head Group'],
                                          columns=['Delta Mass'],
                                          aggfunc="sum"
                                          )
gnps_head_delta_pivot_0h.fillna(0, inplace=True)

gnps_head_delta_pivot_plot = (gnps_head_delta_pivot_72h
                              .sub(gnps_head_delta_pivot_0h)
                              )

#Plot a heatmap of the delta masses and strain data for the MZmine3 processed
#GNPS2 Data

#Create the colormap
colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)]  #R -> G -> B
cmap_name = "white_blue_red"
cmap_wbr = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors)

#Plot the data
gn_deltas_strains_plot = sns.clustermap(gnps_head_delta_pivot_plot,
                                        norm = LogNorm(),
                                        cmap = cmap_wbr,
                                        method = "single",
                                        metric = pearson,
                                        col_cluster=False,
                                        row_cluster=True,
                                        xticklabels = 1,
                                        yticklabels = 1
                                       )

hm = gn_deltas_strains_plot.ax_heatmap.get_position()
gn_deltas_strains_plot.ax_heatmap.set_position([hm.x0,
                                                hm.y0,
                                                hm.width,
                                                1.3*hm.height
                                                ])

#Annotate and format the plot
masses = gn_deltas_strains_plot.ax_heatmap.get_xticklabels()
acyl_annotations = annotationdata.acyl_group_annotations

for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(acyl_annotations[float(mass_number)]
                             + "  "
                             + mass_number
                             )
    except KeyError:
        mass = mass.set_text(mass_number)
gn_deltas_strains_plot.ax_heatmap.set_xticklabels(masses,rotation=90)

for _, spine in gn_deltas_strains_plot.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in gn_deltas_strains_plot.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)