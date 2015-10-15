"""
geneinfo.py: methods for querying, saving
and loading gene information for NIH
database.

Author: Torben Noto
"""

import pandas as pd

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

def load_gene_file(path='.'):
    if isinstance(path, str):
        gene_file = path+'gene_info.csv'
        df = pd.read_csv(gene_file)
        return df
    else:
        raise TypeError("gene no must be a string")

def get_gene_info(path, gene_ids):
    df = load_gene_file(path)
    output = []
    for gene_id in gene_ids:
        if isinstance(gene_id, str):
            if gene_id in df['Entrez']:
                output.append((df[df['Entrez'] == int(gene_id)].as_matrix()[0]))
            else:
                print 'Gene '+gene_id+' not found in NIH database'
        else:
            print str(gene_id)+' must be a str'
    return output
