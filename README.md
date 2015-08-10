# nsaba (NeuroSynth, Allen Brain Atlas)

Methods to analyze genome-scale gene expression data from the Allen Human Brain Atlas in conjunction with fMRI activation maps from Neurosynth


##Setup Instructions

- Download the H0351.2001 file from http://human.brain-map.org/static/download and save the file in a place that you will remember (/Users/User/Documents/Nsaba_analysis/).
This file is genome-scale gene expression at ~900 sites on one human brain. Multi-brain analyses coming soon.

- Download current_data.tar.gz from https://github.com/neurosynth/neurosynth-data This is a bank of thousands of fMRI studies where a significant region of activation was found. These activated regions are matched with how often certain terms were used in these papers (i.e 'aging', 'schizophrenia', 'fear')
Save the decompressed output of this file to somewhere you will remember (/Users/User/Documents/Nsaba_analysis/).


## Example Usage:

    from nsaba import Nsaba

    aba_path = '/Users/User/Documents/Nsaba_analysis/normalized_microarray_donor9861/' #path to ABA data
    ns_path = '/Users/User/Documents/Nsaba_analysis/current_data/' #path to NS data

    Nsaba.aba_load(aba_path)
    Nsaba.ns_load(ns_path)

    entrez_ids = [740, 88]
	my_term = 'alzheimer'
    alz = Nsaba()
    
    alz.get_ge(entrez_id)
    alz.get_ns_act(myterm, thresh=0, method='knn', k=15)
    analysis_matrix = alz.make_ge_ns_mat(myterm, entrez_ids)
    
    # Run Science
    
## Science Sloth!    
![Science Sloth!](http://gifts.worldwildlife.org/gift-center/Images/large-species-photo/large-Three-toed-Sloth-photo.jpg 
    "Science Sloth preparing to do Science!")


