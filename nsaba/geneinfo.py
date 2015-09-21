import urllib2
from bs4 import BeautifulSoup

def gene_info(eid):
    """Pulls gene data based on Entrez ID from the NIH and returns summary"""
    if isinstance(eid, str):
        try:
            page_name = "http://www.ncbi.nlm.nih.gov/gene/?term=" + eid
            page = urllib2.urlopen(page_name)
            soup = BeautifulSoup(page)
            contents = []
            for ana in soup.findAll('dd'):
                if ana.parent.name == 'dl':
                    contents.append(ana.contents)
            return contents[9]
        except IndexError:
            print "%s isn't registered with the NIH" % eid
    else:
        raise TypeError("gene no must be a string")


def gene_id(eid):
    """Pulls gene data based on Entrez ID from the NIH and returns summary"""
    if isinstance(eid, str):
        try:
            page_name = "http://www.ncbi.nlm.nih.gov/gene/?term=" + eid
            page = urllib2.urlopen(page_name)
            soup = BeautifulSoup(page)
            contents = []
            for ana in soup.findAll('dd'):
                if ana.parent.name == 'dl':
                    contents.append(ana.contents)
            return contents[1][0]
        except IndexError:
            print "%s isn't registered with the NIH" % eid
    else:
        raise TypeError("gene no must be a string")
