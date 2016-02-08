# Nsaba (NeuroSynth, Allen Brain Atlas)

Methods to analyze genome-scale gene expression data from the Allen Human Brain Atlas in conjunction with fMRI activation maps from Neurosynth


## Local Setup Instructions

- Download the H0351.2001 file from http://human.brain-map.org/static/download and save the file in a place that you will remember (/Users/User/Documents/Nsaba_analysis/).
This file is genome-scale gene expression at ~900 sites on one human brain.

- Download current_data.tar.gz from https://github.com/neurosynth/neurosynth-data This is a bank of thousands of fMRI studies where a significant region of activation was found. These activated regions are matched with how often certain terms were used in these papers (i.e 'aging', 'schizophrenia', 'fear')
Save the decompressed output of this file to somewhere you will remember (/Users/User/Documents/Nsaba_analysis/).


## Example Usage:

    from nsaba.nsaba import Nsaba

Specifying paths to Allen Brain Atlas and Neurosynth data files.

    aba_path = '/Users/User/Documents/Nsaba_analysis/normalized_microarray_donor9861/'
    ns_path = '/Users/User/Documents/Nsaba_analysis/current_data/'

    Nsaba.aba_load(aba_path)
    Nsaba.ns_load(ns_path)

Declaring Entrez IDs for genes of interest: in this case [1813](http://www.ncbi.nlm.nih.gov/gene/?term=1813) and 
[1816](http://www.ncbi.nlm.nih.gov/gene/?term=1816).

    entrez_ids = [1813, 1816]
	my_term = 'attention'
    att = Nsaba()


Generating gene expression and term tf-idf statistics.

    att.estimate_aba_ge(entrez_ids)
    att.estimate_ns_act(my_term, knn_args={'n_neighbors':10})
    analysis_matrix = att.matrix_builder(my_term, entrez_ids)
    
Nsaba comes with an analysis module to aid with large-scale significance testing and information
retrieval.

    import nsaba.analyis as na
     
    ....
    
    att_stats  = na.NsabaAnalysis(att)
    tttest_metrics = att_stats.term_ge_ttest_multi('reward',quant=90)
    att_stats.fetch_gene_descriptions(tttest_metrics)
    
See for this [notebook](https://github.com/voytekresearch/nsaba/blob/master/notebooks/demos/Nsaba_Demonstration.ipynb) for a full demonstration.    

## Futures

 * A [pip](https://pypi.python.org/pypi/pip) package is currently under development as well as a 
 hosted [binder](http://mybinder.org/).
 
 * Multi-brain analysis.
 
 * Full suite of demonstration and example notebooks.