# Multi-band-ICA-for-EEG-and-MEG
Multi-band ICA is a combination of a filter bank, PCA and ICA for artifact removal and component extraction in EEG (and MEG).
The conventional wide-band ICA is only able to find the strongest sources and data from weaker sources are not always recovered. The multi-band ICA is superior to the conventional ICA in terms of separation performance, number of the extracted components, and the quality of the reconstructed sources. 
The background data provided in this repository is from a resting state of one person. Three sources are simulated at 5 Hz, 12 Hz, and 30 Hz with SNRs of 0.1, 0.04, and 0.04, respectively.

How to run:
1- Download the fieldtrip toolbox and add its path to your MATLAB (folders and subfolders)
ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/
2- Download the codes from this repository and run the code MultiBandICA_Demonstration.m.

It runs the multi-band ICA and the wide-band ICA on sample dataset provided on here. 

Question:
y.jonmo@auckland.ac.nz
