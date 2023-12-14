#Receives mass spectrometry data of N-acyl amides from bacterial culture data,  
#and visualizes through various plotting methods

#Data originally obtained from Emily Gentry's MSV000084475 Massive Dataset, 
#FCM Media, and queried for N-acyl amides with MassQL 
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
from correlationmetrics import pearson
from queryretriever import prot2_retrieve
from queryretriever import gnps_retrieve

#Data Input

#Read in queries from Proteomics2 and GNPS2
prot2_input_data = pd.read_csv("/content/MSV000084475_FCM_Proteomics2_Data.csv")

gnps_input_data = pd.read_csv("/content/MSV000084475_FCM_GNPS2_Data.csv")

#Read in quantification table for MZMine3 processed data (from GNPS2)
quant_table = pd.read_csv("/content/MSV000084475_FCM_quantification_table.csv")

#Read in metadata file
culture_metadata = ({row[0] : row[1] for _,
                    row in
                    pd.read_csv("/content/MSV000084475_FCM_culture_metadata.tsv"
                    , delimiter = "\t").iterrows()
                    })

species = ({row[0] : row[11] for _,
            row in pd.read_csv("/content/MSV000084475_strain_metadata.tsv",
                               delimiter = "\t").iterrows()})
genus = ({row[0] : row[10] for _,
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
        if(species[quant_table[column].name[:-10]] != "not applicable"):
            taxon = species[quant_table[column].name[:-10]].replace(' ', '\\ ')
            strain = ("$\it{"
                      + taxon
                      + "}$ "
                      + culture_metadata[quant_table[column].name[:-10]]
                      )
        elif(genus[quant_table[column].name[:-10]] != "not applicable"):
            taxon = genus[quant_table[column].name[:-10]].replace(' ', '\\ ')
            strain = ("$\it{"
                      + taxon
                      + "}$ "
                      + culture_metadata[quant_table[column].name[:-10]]
                      )
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

#Download unprocessed files from Proteomics2

prot2_successful_heads = []
prot2_scan_num = []
for item in prot2_input_data.iterrows():
    prot2_scan_num.append(prot2_retrieve(item[1][1], "_" + item[1][0]))
    prot2_successful_heads.append(item[1][0])

#Change Proteomics2 data into two column table, remove empty hits
hits = pd.DataFrame({"Heads":prot2_successful_heads,"Scans":prot2_scan_num})
hits = hits[(hits != 0).all(1)]

#Extract culture condition and intesity data
prot2_column_headers = ["Head Group", "Culture Condition", "Intensity"]
prot2_plot_df = pd.DataFrame(columns = prot2_column_headers)
prot2_empty_files = []

i = 0
for row in hits.iterrows():
    data = pd.read_csv("/content/prot2_queried_file_" + row[1][0] + ".tsv",
                       delimiter="\t"
                       )
    for another_row in data.iterrows():
        intensity = another_row[1][7]
        try:
            culture_cond = culture_metadata[((another_row[1][12])[:-2]
                                             + "XML"
                                             )]
            if culture_cond.find("t=") != -1:
                culture_cond = culture_cond[culture_cond.find("t=") + 2:]
        except KeyError:
            prot2_empty_files.append((another_row[1][12])[:-2] + "XML")
        prot2_plot_df.loc[i] = [row[1][0], culture_cond, intensity]
        i = i + 1

#Plot the intensities of Proteomics2 head group data against culture conditon

rcParams['figure.figsize'] = 12.5,8
sns.stripplot(data=prot2_plot_df,
              x="Culture Condition",
              y="Intensity",
              hue="Head Group",
              size = 5,
              dodge=True,
              order = ['72h','0h','media blank', 'extraction blank']
              )

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
        a_row = processed_quant_table.loc[(processed_quant_table["row ID"]
                                           == another_row[1][11]
                                           )]
            #(formatted_quant_table["row ID"] == nested_row[1][11]) &
            #(formatted_quant_table["row m/z"] == nested_row[1][8])
        k = 3
        delta = a_row.iloc[0,1] - annotationdata.head_masses[row[1][0]]
        for value in a_row.iloc[0][3:]:
            strain_cond = processed_quant_table.iloc[:,k].name[:-4]
            culture_cond = strain_cond
            if culture_cond.find("t=") != -1:
                culture_cond = culture_cond[culture_cond.find("t=") + 2:]
            gnps_plot_df.loc[j] = ([row[1][0],
                                    culture_cond,
                                    strain_cond,
                                    value,
                                    delta
                                    ])
            k = k + 1
            j = j + 1

#Subtract head group masses from total mass to obtain acyl group masses

#Plot the intensities of Proteomics2 head group data against culture conditon
sns.stripplot(data=gnps_plot_df,
              x="Culture Condition",
              y="Peak Area",
              hue="Head Group",
              size = 5,
              dodge=True,
              order = ['72h','0h','media blank', 'extraction blank']
              )

#Order the strains vs. head group intensity plot labels

culture_order = list(processed_quant_table.columns)
for i in range(len(culture_order)):
    culture_order[i] = culture_order[i][:-4]

culture_order.remove("ro")
culture_order.remove("row")
culture_order.remove("row retention ")
culture_order.remove("t=0h")
culture_order.remove("t=72h")

culture_order.remove("media blank")
culture_order.remove("extraction blank")

culture_order = list(set(culture_order))
culture_order.sort()
culture_order.append("media blank")
culture_order.append("extraction blank")

#Plot strains vs. head group intensity

rcParams['figure.figsize'] = 65,8.27
a = sns.stripplot(data=gnps_plot_df,
                  x="Strain",
                  y="Peak Area",
                  hue="Head Group",
                  size=5,
                  dodge=True,
                  order=culture_order
                  )
a.tick_params(axis="x", rotation=90)

#Aggregate unique delta mass and head group data

#Create a list of unique delta masses and their respective head groups
unique_entries = []
for row in gnps_plot_df.iterrows():
    if row[1][1] != 'extraction blank' and row[1][1] != 'media blank':
        unique_entries.append((row[1][0], row[1][1], row[1][4]))
unique_entries = list(set(unique_entries))
heads_deltas_cols = ["Head Group","Delta Mass","Relative Peak Area"]
heads_deltas_df = pd.DataFrame(columns = heads_deltas_cols)

#Iterate through the main dataframe and aggregate and normalize data
#belonging to each unique N-acyl-lipid
i = 0
for item in unique_entries:
    relative_peak_area = 0
    relative_peak_area_end = 0
    relative_peak_area_zero = 0
    endpoint_number = 0
    zero_number = 0
    for row in gnps_plot_df.iterrows():
        if row[1][0] == item[0] and row[1][4] == item[2]:
            if row[1][1] == '72h':
                relative_peak_area_end = relative_peak_area_end + row[1][3]
                endpoint_number = endpoint_number + 1
            elif row[1][1] == '0h':
                relative_peak_area_zero = relative_peak_area_zero + row[1][3]
                zero_number = zero_number + 1
    relative_peak_area = ((relative_peak_area_end / (endpoint_number + 0.01))
                           - (relative_peak_area_zero / (zero_number + 0.01)))
    heads_deltas_df.loc[i] = [item[0], round(item[2]), relative_peak_area]
    i = i + 1

heads_deltas_df = heads_deltas_df.drop_duplicates()

#Convert data to a pivot table to allow for heatmap plotting
heads_deltas_plot_df = pd.pivot_table(heads_deltas_df,
                                         values='Relative Peak Area',
                                         index=['Head Group'],
                                         columns=['Delta Mass'],
                                         aggfunc="sum"
                                         )
heads_deltas_plot_df.fillna(0)

#Plot normalized N-acyl-lipid aggregation data

#Define color normalization parameters
vcenter = 0
vmin, vmax = -2000, 2000
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

#Plot data
rcParams['figure.figsize'] = 90,8.27
d = sns.clustermap(heads_deltas_plot_df,
                   cmap = "vlag",
                   norm=normalize,
                   col_cluster=False,
                   row_cluster=False,
                   xticklabels = 1,
                   yticklabels = 1
                   )

#Annotate and format the plot
masses = d.ax_heatmap.get_xticklabels()
for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(
            annotationdata.acyl_group_annotations[round(float(mass_number))]
            + "  "
            + str(round(float(mass_number)))
            )
    except KeyError:
        mass = mass.set_text(str(round(float(mass_number))))
d.ax_heatmap.set_xticklabels(masses)

for _, spine in d.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in d.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

#Create a plot of strain against acyl group for the particularly abundant
#head group Phenylalanine

phenylala_data = gnps_plot_df[gnps_plot_df["Head Group"] == "Phenylalanine"]
phenylala_data = phenylala_data[(phenylala_data["Strain"]
                                 == "extraction blank"
                                 )
                                 | (phenylala_data["Strain"]
                                    == "media blank"
                                    )
                                 | (phenylala_data["Peak Area"]
                                    != 0
                                    )]

phenylala_plot_df = pd.pivot_table(phenylala_data,
                                   values='Peak Area',
                                   index=['Strain'],
                                   columns=['Delta Mass'],
                                   aggfunc="sum"
                                   )

#Plot acyl group against strain data for Phenylalanine

#Plot data
rcParams['figure.figsize'] = 90,8.27
c = sns.clustermap(phenylala_plot_df,
                   cmap = "vlag",
                   norm=normalize,
                   col_cluster=False,
                   row_cluster=False,
                   xticklabels = 1,
                   yticklabels = 1
                   )

#Annotate and format the plot
masses = c.ax_heatmap.get_xticklabels()
for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(
            annotationdata.acyl_group_annotations[round(float(mass_number))]
            + "  "
            + str(round(float(mass_number)))
            )
    except KeyError:
        mass = mass.set_text(str(round(float(mass_number))))
c.ax_heatmap.set_xticklabels(masses)

for _, spine in c.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in c.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)