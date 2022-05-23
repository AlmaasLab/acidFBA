# acidFBA

The acidFBA approach expands on the enzyme-constrained framework of [GECKO](https://github.com/SysBioChalmers/GECKO) (GEM with Enzymatic Constraints using Kinetic and Omics data) by also accounting for the usage of proteinogenic amino acids. This repository contains all data and code necessary to construct an acidFBA-GEM of the concensus enzyme-constrained yeast GEM Yeast v. 8.3.4 and reproduce the results of the paper: "Phenotypic response of yeast metabolic network to availability of proteinogenic amino acids", Simensen et al. 2022.

## Requirements
- Matlab (>= 2020a)
- RAVEN Toolbox 2
- COBRA Toolbox v3.0
- python 3 (>= 3.5)
- pip package manager

## Installation
No installation is required. You can simply clone or download the repository, install the necessary python packages and run the scripts.

To clone the repository, run the following in your command line interface at a chosen directory

`git clone https://github.com/AlmaasLab/acidFBA.git`

Use pip to install the necessary Python packages listed in *requirements.txt*

`pip install -r requirements.txt`

## Usage
- Running the `initAcidFBA.m` script, the acidFBA-GEM is constructed and used to reproduce all the results presented in the article.
- The figures found under the `figures` folder can also be created by running the Jupyter notebook `fig_notebook.ipynb`.

# Contributors
- [Vetle Simensen](https://www.ntnu.no/ansatte/vetle.simensen), Norwegian University of Science and Technology, Norway
