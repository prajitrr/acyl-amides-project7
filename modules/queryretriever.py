#Receives a link to either a Proteomics2 or GNPS2 job query and downloads the
#the results file of the job as a tsv
#Specify a file label to uniquely label each file
#Contact Prajit Rajkumar (prajkumar@ucsd.edu) for any questions

#Note that you may have to install some of the imports using 
#pip or another installer
import time
import ast
import numpy as np
import pandas as pd
import requests
from urllib.request import urlopen
from bs4 import BeautifulSoup
from urllib.request import urlopen

def prot2_retrieve(url, file_label):
    """
    Retrieve the result file from a Proteomics2 job query in tsv format.

    :param str url: URL to a Proteomics2 job query
    :param str file_label: A label used to uniquely identify the file
    :return: The number of scans in the result file; returns 0 for empty files
    :rtype: int
    """
    dataset_loader = url + "&view=extract_results"
    load_url = requests.get(dataset_loader, allow_redirects=True)
    not_loaded = True
    while not_loaded:
        opened_dataset = urlopen(dataset_loader)
        site_text = opened_dataset.read().decode("UTF-8")
        if "cellpadding=\"" in site_text:
            not_loaded = False
        elif "Job failed due to" in site_text:
            return 0
    task_ID = url[url.find("task=") + 5:]
    processed_url = ("https://proteomics2.ucsd.edu/ProteoSAFe/"
                     + "QueryResult?task=" + task_ID + "&file=extract_results"
                     + "-main_dynamic_extracted.db&pageSize="
                     + "-1&offset=0&query=&_=")
    r = requests.get(processed_url, allow_redirects=True)
    page = r.content
    soup = BeautifulSoup(page, 'html.parser')
    text = soup.find_all(string=True)
    text = "".join(text)
    if (text[0:35] == "Apache Tomcat/6.0.24 - Error report"):
        return 0
    text = text[text.find("row_data") + 10:-1]
    text = text.replace("\n", "")
    text = ast.literal_eval(text)
    text = pd.DataFrame(text)
    text.to_csv("prot2_" + "queried_file" + file_label + ".tsv", sep="\t")
    return len(text)

def gnps_retrieve(url, file_label):
    """
    Download the result file from a GNPS2 job query in tsv format.

    :param str url: URL to a GNPS2 job query
    :param str file_label: A label used to uniquely identify the file
    :return: The number of scans in the result file; returns 0 for empty files
    :rtype: int
    """
    task_ID = url[url.find("task=") + 5:]
    processed_url = ("https://gnps2.org/resultfile?task=" + task_ID 
                     + "&file=nf_output/msql/merged_query_results.tsv")
    r = requests.get(processed_url, allow_redirects=True)
    file_name = "gnps_" + "queried_file" + file_label + ".tsv"
    open(file_name, "wb").write(r.content)
    return len(pd.DataFrame(r.content))
