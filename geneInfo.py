import urllib2
from bs4 import BeautifulSoup


def gene_info(gene_no):
    """Pulls gene data from the NIH and returns a blurb about a gene. BeautifulSoup4 and urllib2 dependencies"""
    if isinstance(gene_no, str):
        try:
            page_name = "http://www.ncbi.nlm.nih.gov/gene/?term=" + gene_no
            page = urllib2.urlopen(page_name)
            soup = BeautifulSoup(page)
            contents = []
            for ana in soup.findAll('dd'):
                if ana.parent.name == 'dl':
                    contents.append(ana.contents)
            return contents[9]
        except IndexError:
            print "This gene isn't registered with the NIH"
    else:
        raise TypeError("gene no must be a string")