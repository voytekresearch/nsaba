"""
geneinfo.py: methods for querying, saving
and loading gene information for NIH
database.

Author: Torben Noto & Simon Haxby
"""

import pandas as pd
import os
import re
import random
import urllib2
from time import sleep
from collections import namedtuple

from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC


def gene_info(eid):
    """
    Pulls gene data based on Entrez ID from the NIH and returns summary.

    Parameters
    ----------
    eid : int
        Entrez ID of interest

    Returns
    -------
    (gene_name, gene_description) : (string, string)
        gene_description contains a string of appropriately 30-80
        characters describing function, relevance and attribution
        of the gene specified by eid.
    """
    if isinstance(eid, str):
        try:
            page_name = "http://www.ncbi.nlm.nih.gov/gene/?term=" + eid
            page = urllib2.urlopen(page_name)
            sleep(1+random.random())
            soup = BeautifulSoup(page, 'lxml')
            contents = []
            for ana in soup.findAll('dd'):
                if ana.parent.name == 'dl':
                    contents.append(ana.contents)
            gene_name = contents[1][0]
            gene_description = contents[9]
            if not len(gene_description[0]) > 1:
                gene_description = 'No description found'
            return gene_name, gene_description
        except IndexError:
            print "%s isn't registered with the NIH" % eid
            return 'No Gene identification found', 'No description found'
    else:
        raise TypeError("gene no must be a string")

def load_gene_file(path='.'):
    """
    Loads file containing gene descriptions of genes specified by
    their Entrez IDs from: http://www.ncbi.nlm.nih.gov/gene/ .

    Parameters
    ----------
    path : str, optional
        Specifies path to gene_info.csv.

    Returns
    -------
    pandas.DataFrame
        Returns a DataFrame where each row contains three fields:
        'Entrez', 'Gene Name' and 'Gene Description'. Where 'Entrez'
        specifies the gene's Entrez ID and the last two fields are
        of the same form as gene_info()'s returns.

        NOTE: This assumes that correct CSV has been loaded.
    """
    if isinstance(path, str):
        gene_file = os.path.join(path, 'gene_info.csv')
        df = pd.read_csv(gene_file)
        return df
    else:
        raise TypeError("Gene-file path must be a string")

def get_gene_info(path, gene_ids):
    """
    Extracts gene information from DataFrame created by
    load_gene_file() for specific genes based on list
    of Entrez IDs.

    Parameters
    ---------
    path : str
        Specifies path to gene_info.csv.

    gene_ids : list [ int ]
        List of Entrez IDs of gene descriptions to be fetched

    Returns
    -------
    output : list [ gi_tuple (long, str, u-str) ]
        Returns a list of gene information for specified
        Entrez IDs in form: ('Entrez', 'Gene Name' 'Gene Description').
    """
    gi_tuple = namedtuple("gi_tuple", "entrez name description")
    df = load_gene_file(path)
    output = []
    for gene_id in gene_ids:
        if gene_id in df['Entrez']:
            gi = df[df['Entrez'] == gene_id].as_matrix()[0]
            output.append(gi_tuple(gi[0], gi[1], gi[2]))
        else:
            print 'Gene %s not found in NIH database' % gene_id
    return output

def fetch_entrez_ids(term, id_num):
    """
    Returns Entrez IDs of genes most strongly associated a term from
    www.genecard.org.

    Parameters
    ----------
    term : str
        Term for

    id_num : str/int
        The number of highest scored Entrez IDs for the specified term.

    Returns
    -------
    entrez_ids : list
        A list of 'id_num' top Entrez IDs associated with 'term'.
    """
    entrez_ids = []
    strs = ["http://www.genecards.org/Search/Keyword?startPage=0&queryString=",
            term, "&pageSize=", str(id_num)]
    search_url = ''.join(strs)
    driver = webdriver.PhantomJS()
    driver.get(search_url)
    print (search_url)
    try:
        # Waits for DOM to render the search results table before accessing elements.
        check = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.ID, "searchResults")))
    finally:
        search_table = driver.find_element_by_id("searchResults")

    top_genes = search_table.find_elements_by_class_name("gc-gene-symbol")
    gene_urls = []
    for gene in top_genes:
        el = gene.find_element_by_tag_name("a")
        gene_urls.append(el.get_attribute("href"))

    for gene_url in gene_urls:
        driver.get(gene_url)
        try:
            # Waits for DOM to render the page table before accessing elements.
            check = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.CLASS_NAME, "gc-subsection")))
        finally:
            subsections = driver.find_elements_by_class_name("gc-subsection")

        text = subsections[1].text
        result = re.search("Entrez\sGene:\s([0-9]*)", text)
        entrez_ids.append(result.group(1))
        print (entrez_ids[-1])

    driver.close()

    return entrez_ids