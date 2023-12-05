#Receives a list of molecules and outputs a datframe with SMILES and links
#to structure drawings of the molecules
#Note that you may have to expand some names in the input to their long form 
#(eg GABA -> Gamma aminobutyric acid)
#Contact Prajit Rajkumar (prajkumar@ucsd.edu) for any questions

#Note that you may have to install some of the imports using 
#pip or another installer
import numpy as np
import pandas as pd
from urllib.request import urlopen
from urllib.parse import quote

def retriever(list_of_names):
    """
    Retrieve the SMILES ID and images of a specified list of molecules.

    :param list_of_names: List of n molecules, preferably using IUPAC 
                          naming conventions
    :type list_of_names: list[str] or None
    :return: An n by 3 dataframe containing the SMILES IDs and images 
             of each molecule
    :rtype: pandas.core.frame.DataFrame
    """
    structures_df = pd.DataFrame(columns=["Molecule","SMILES", "Image"])
    i = 0
    for item in list_of_names:
        formatted_item = item.replace(" ","%20")
        try:
            url = ('http://cactus.nci.nih.gov/chemical/structure/' 
                   + formatted_item + '/smiles')
            smiles_id = urlopen(url).read().decode('utf8')
            image_url = ('https://cactus.nci.nih.gov/chemical/structure/' 
                         + smiles_id + '/image')
        except:
            smiles_id = 'Unavailable'
            image_url = 'Unavailable'
        structures_df.loc[i] = [item, smiles_id, image_url]
        i = i + 1
    return structures_df
