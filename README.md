# acidFBA

The acidFBA approach expands on the enzyme-constrained framework of [GECKO](https://github.com/SysBioChalmers/GECKO) (GEM with Enzymatic Constraints using Kinetic and Omics data) by also accounting for the usage of proteinogenic amino acids. This repository contains the necessary functions to build an acidFBA-GEM from any enzyme-constrained GEM in the GECKO format. Additionally, the repository contains all data and code necessary to construct an acidFBA-GEM of the concensus enzyme-constrained yeast GEM Yeast v. 8.3.4 and reproduce the results of the paper: "Phenotypic response of yeast metabolic network to availability of proteinogenic amino acids", Simensen et al. 2022.

## Requirements
- Matlab (>= 2020a)
- RAVEN Toolbox 2
- COBRA Toolbox v3.0
- Gurobi (>= 8.1.1)
- python 3 (>= 3.5)
- pip package manager

## Installation
No installation is required. You can simply clone or download the repository, install the necessary python packages and run the scripts.

To clone the repository, run the following in your command line interface at a chosen directory

`git clone https://github.com/AlmaasLab/acidFBA.git`

Use pip to install the necessary Python packages listed in *requirements.txt*

`pip install -r requirements.txt`

## Usage
- Running the `initAcidFBA.m` script, the yeast acidFBA-GEM is constructed and used to reproduce all results presented in the article.
- The figures found under the `/figures` folder can be recreated via the Jupyter notebook `fig_notebook.ipynb`.
- To create an acidFBA-GEM from any other enzyme-constrained GEM (GECKO framework), the following is required:
    - In the `/proteins` folder, exchange the yeast-specific `protSeqs.xlsx` to one of your selected organism. This file should contain two columns: (1) identifiers used as part of the metabolite IDs of the protein pseudo-metabolites in the model (`Identifier`, e.g., UniProt identifiers), (2) amino acid sequences in single letter format (`Sequence`).
    - In `getParameters.m`, set the following fields to your organism-specific values: average enzyme saturation factor `sigma`, total cellular protein content `Ptot` [g/gDW], mass fraction of proteins accounted for in the model `f` [g/g], and the protein pool exchange reaction identifier `protPool`. Allows for the use of `constrainProtPool.m` to appropriately constrain the protein pool of the model.
    - Run the function `createAcidFBAmodel(model)` where the variable `model` is a model struct of the enzyme-constrained GEM in GECKO format.

# Contributors
- [Vetle Simensen](https://www.ntnu.no/ansatte/vetle.simensen), Norwegian University of Science and Technology, Norway
