# Nsaba (NeuroSynth, Allen Brain Atlas)

Methods to analyze genome-scale gene expression data from the Allen Human Brain Atlas in conjunction with fMRI activation maps from Neurosynth


## Local Setup Instructions

- Download the H0351.2001 file from http://human.brain-map.org/static/download and save the file in a place that you will remember (/Users/User/Documents/Nsaba_analysis/).
This file is genome-scale gene expression at ~900 sites on one human brain.

- Download current_data.tar.gz from https://github.com/neurosynth/neurosynth-data This is a bank of thousands of fMRI studies where a significant region of activation was found. These activated regions are matched with how often certain terms were used in these papers (i.e 'aging', 'schizophrenia', 'fear')
Save the decompressed output of this file to somewhere you will remember (/Users/User/Documents/Nsaba_analysis/).


## Example Usage:

See for this [notebook](https://github.com/voytekresearch/nsaba/blob/master/notebooks/demos/Nsaba_Long_Demo.ipynb) for a full demonstration.    



## Futures

 * A [pip](https://pypi.python.org/pypi/pip) package is currently under development as well as a 
 hosted [binder](http://mybinder.org/).
 
 * Multi-brain analysis.
 
 * Full suite of demonstration and example notebooks.
