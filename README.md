# TRUHiC
## Introduction
TRUHiC is a Hi-C data resolution enhancement method that integrates a customized and lightweight transformer architecture embedded into a U-2 Net architecture to augment low-resolution Hi-C data for the characterization of 3D chromatin structure. 

This repository contains codes and processed files for the manuscript entitled *"TRUHiC: A TRansformer-embedded U-2 Net to enhance Hi-C data for 3D chromatin structure characterization."*. (link to be added)


## Getting started
Codes for the main experimental analysis are provided in the <code>Models.zip</code> file with instructions included in a readme file inside. All required input files for a demo can be found in the <code><b>Data</b></code> folder zipped and can be extracted using the 7zip tool.   

**Please ignore the following instructions for running the experiments at this point**

To get started, users can download the scripts and run them on their local machines. To run this on the HPC, after connecting to the user's HPC account, install the libraries as suggested in the Getting Started section. The user can run the same code on their HPC server. XXX

### Installation
TRUHiC can be downloaded by
```
git clone https://github.com/shilab/TRUHiC
```

### Prerequisites:
Python >= 3.7.3  
Jupyterlab >= 4.2.3

Install required dependencies 
```
pip3 install pandas==1.2.4 numpy==1.20.2 scipy==1.7.3 matplotlib==3.5.3 statsmodels==0.13.5 seaborn==0.11.1 scikit_posthocs==0.8.1 jupyterlab
```

Ensure that the virtual environment meets the following dependencies:  
Pandas 1.2.x, Numpy 1.20.x, SciPy 1.7.x, Matplotlib 3.5.x, statsmodels 0.13.x, seaborn 0.11.x, scikit_posthocs 0.8.x. 

Users can download the project repository and start the jupyter lab to experiment with the analysis
```
git clone https://github.com/shilab/TRUHiC.git
cd XXX
cd XXX
```

The <code><b>Data</code></b> folder contains the necessary datasets that are needed for running the main analyses included in our study. A *README* file for the detailed description of each file can be found under the data folder.

Please note that the scripts are specifically designed and organized for this study publication. All the input files and formats are specified in the scripts. Users are welcome to download and run the provided scripts on their own machines to replicate our results. It is possible that the programs may not run on the user's device due to environmental differences or bugs. Therefore, to use the scripts with the user's own data, please consider this repository as an experimental notebook and update the respective directory paths and input files accordingly. 

## Contact
We welcome your questions, suggestions, requests for additional information, or collaboration interests. Please feel free to reach out to us via the following email addresses and we will respond as soon as possible:  
:email: Chong Li:   tun53987@temple.edu or lichong0710@gmail.edu (personal email)  
:email: Mohammad Erfan Mowlaei:   mohammad.erfan.mowlaei@temple.edu  
:email: Dr. Mindy Shi:   mindyshi@temple.edu

## References
#### If you find our results useful in your research, please cite our work as:
Chong Li, Mohammad Erfan Mowlaei, Human Genome Structural Variation Consortium (HGSVC), HGSVC Functional Analysis Working Group, Vincenzo Carnevale, Sudhir Kumar, Xinghua Shi. “TRUHiC: A Transformer-embedded U-2 Net to enhance Hi-C data for chromatin structure characterization.”

