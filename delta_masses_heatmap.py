#Receives mass spectrometry data of n-acyl lipids,  
#visualizes it as a heatmap, and outputs information on synthesis targets
#Code was originally run in Google Colaboratory 
#Running on other platforms may require reformatting heatmap dimensions
#Contact Prajit Rajkumar (prajkumar@ucsd.edu) for any questions

#Note that you may have to install some of the imports using 
#pip or another installer
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
from matplotlib.colors import LogNorm, Normalize
import matplotlib.colors as pltclr
import scipy.spatial.distance as pdist

import annotationdata
from correlationmetrics import pearson
from structureretriever import retriever

#Input

#Read in head group names and masses
head_group_masses = pd.read_csv("/content/headgroupmasses.csv")

#Read in mass data
heatmap_data = pd.DataFrame()
for i in range(57):
    j = i + 1
    df = pd.read_csv("/content/deltamassfile" + str(j) + ".tsv", sep='\t')
    try:
        df = df["_dyn_#precmz"]
        df = df.rename("precmz")
    except:
        df = df["precmz"]
    df = df.map(
        lambda precmz : precmz - head_group_masses.loc[i]["SubtractMass"]
        )
    df = df.rename(head_group_masses.loc[i]["Head"])
    heatmap_data = pd.concat([heatmap_data, df], axis = 1)

heatmap_data = heatmap_data.transpose()

#Reformat data to two column table
heads = []
deltas = []
a = 0
for row in heatmap_data.itertuples(True):
    for i in range(len(row) - 1):
        heads.append(row[0])
        deltas.append(row[i + 1])

two_col_data = pd.DataFrame({"Heads": heads, "Deltas": deltas})
two_col_data["Deltas"] = two_col_data["Deltas"].round()
two_col_data.dropna(axis=0)

pre_plotted_data = two_col_data.pivot_table(index='Heads',
    columns='Deltas', aggfunc=lambda x: len(x))
pre_plotted_data.fillna(0, inplace=True)

#Plot heatmap

#Create the colormap
colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)]  #R -> G -> B
#Discretizes the interpolation into bins
n_bins = 8
cmap_name = "white_blue_red"
cmap_wbr = pltclr.LinearSegmentedColormap.from_list(cmap_name, 
    colors, 
    N=n_bins
    )

#Create the plot
a = sns.clustermap(pre_plotted_data, 
    cmap = cmap_wbr, 
    method = "single", 
    metric = pearson, 
    norm=LogNorm(), 
    col_cluster=False, 
    cbar_pos = (3.88, 0.2, 0.05, 0.5), 
    xticklabels = 1,
    yticklabels = 1, 
    cbar_kws = {'location': 'right', 'aspect': 30, 'anchor': (1, 10)}, 
    dendrogram_ratio = (0.99, 0.1))
hm = a.ax_heatmap.get_position()
a.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*400, hm.height])
col = a.ax_col_dendrogram.get_position()

col = a.ax_col_dendrogram.get_position()
masses = a.ax_heatmap.get_xticklabels()

#Annotate and format the plot
for mass in masses:
    mass_number = mass.get_text()
    try:
        mass = mass.set_text(
            annotationdata.acyl_group_annotations[float(mass_number)] 
            + "  " 
            + mass_number
            )
    except KeyError:
        mass = mass.set_text(mass_number)
a.ax_heatmap.set_xticklabels(masses)
a.ax_col_dendrogram.set_position([
    col.x0, 
    col.y0, 
    col.width*0.5, 
    col.height*0.25
    ])

for _, spine in a.ax_heatmap.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

for _, spine in a.ax_cbar.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)

a.ax_col_dendrogram.set_position([hm.x0, hm.y0, hm.width*8, hm.height])

plt.tight_layout()
plt.show()

#Retrieve specific data for synthesis information

#List out all unmodified and modified decarxboxylated amino acids in dataset
unmodified_decarboxylated_AA = ([
    "Tryptamine", 
    "Tyramine", 
    "Serotonin", 
    "Putrescine", 
    "Histamine", 
    "GABA", 
    "Dopamine", 
    "Cadaverine", 
    "Agmatine", 
    "2-phenethylamine"
    ])
modified_decarboxylated_AA = ([
    "N-carbamoylputrescine", 
    "N-acetylputrescine", 
    "N-acetylcadaverine", 
    "5-Methoxytryptamine",
    "3-methylhistamine",
    "3-Methoxytyramine"
    ])

#For synthesis purposes, chains of length 10 or less were considered
truncated_pre_plotted_data = pre_plotted_data.loc[:,:"154.0"]
largest_deltas_data = pd.DataFrame(columns=([
    "Head", 
    "1st", 
    "1st_SMILES", 
    "2nd", 
    "2nd_SMILES", 
    "3rd", 
    "3rd_SMILES"
    ]))

#Output a table that contains names of 3 most common N-acylamides of 
#respective head groups and their SMILES
i = 0
names = annotationdata.acyl_group_names
annotations = annotationdata.acyl_group_annotations
for item in unmodified_decarboxylated_AA:
    largest_deltas = truncated_pre_plotted_data.loc[item].nlargest(3)
    if item != "GABA":
        formatted_item = item
    else:
        formatted_item = "Gamma-aminobutyric acid"
    first = (
        "N-" 
        + names[annotations[largest_deltas.index[0]]] 
        + formatted_item
        )
    second = (
        "N-" 
        + names[annotations[largest_deltas.index[1]]] 
        + formatted_item
        )
    third = (
        "N-" 
        + names[annotations[largest_deltas.index[2]]] 
        + formatted_item
        )
    smiles = retriever([first, second, third])
    largest_deltas_data.loc[i] = ([
        item, 
        first, 
        smiles.iloc[0,1], 
        second, 
        smiles.iloc[1,1], 
        third, 
        smiles.iloc[2,1]
        ])
    i = i + 1