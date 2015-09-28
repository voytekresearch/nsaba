import urllib2
from bs4 import BeautifulSoup
from time import sleep
import random


def gene_info(eid):
    """Pulls gene data based on Entrez ID from the NIH and returns summary"""
    if isinstance(eid, str):
        try:
            page_name = "http://www.ncbi.nlm.nih.gov/gene/?term=" + eid
            page = urllib2.urlopen(page_name)
            sleep(1+random.random())
            soup = BeautifulSoup(page)
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
