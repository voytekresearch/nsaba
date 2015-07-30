# nsaba (NeuroSynth, Allen Brain Atlas)

Methods to analyze human gene expression data from Neurosynth in conjunction with fMRI activation maps from Neurosynth
(Allen and Yarkoni sitting in a tree, C-O-D-I-N-G.)

### Example Usage:

    from nsaba import Nsaba

    ns_path = 'Big Project/NS'
    aba_path = 'Big Project/ABA/Files'

    Nsaba.aba_load(aba_path)
    Nsaba.ns_load(ns_path)

    entrez_ids = ['7889', '3', '24']
    alz = Nsaba()
    alz.get_ge(entrez_ids)
    alz.get_ns_act('alzheimer', .005)
    mat = alz.make_ge_ns_mat()

    # Run Science
