%% Initializes the construction and analysis of the yeast acidFBA model
% Script for constructing and analyzing an acidFBA-implementation of the
% yeast GEM, Yeast 8.3.
%
% .. Authors:
%       - Vetle Simensen 25/01/21


%% Build the acidFBA model
% Define solver, as well as feasibility and optimality tolerance 
changeCobraSolver('gurobi','all');
changeCobraSolverParams('LP','feasTol',1.0000e-06);
changeCobraSolverParams('LP','optTol',1.0000e-06);

% Read model
currDir = cd;
cd models
fprintf('Reading SBML file ...');
geckoPool = importModel('ecYeastGEM_batch.xml');
cd(currDir);
fprintf('DONE \n');

% Construct the acidFBA model
fprintf('Constructing acidFBA model ...');
[acidFBAmodel,acidFBAuniprot] = createAcidFBAmodel(geckoPool);
fprintf('DONE \n');


%% Correlation with experimental data (Fig. 1 and S1 Fig.)
aminoAcidCorrelation(acidFBAmodel,acidFBAuniprot)


%% Growth rate-dependency of amino acid usage (Fig. 2)
aminoAcidGrowthRates(acidFBAmodel)


%% Sample nutrient environments (Fig. 3, 4)
nSamples = 5000;
[aaMatDflt,aaMinDflt,aaMaxDflt] = sampleNutrientEnvironments(acidFBAmodel,nSamples,'Default');


%% Species-specific amino acid distributions (Fig. 5) 
minAminoAcidDiff(acidFBAmodel);


%% Sample nutrient environments, universal amino acid distribution (S2 Fig)
fprintf('Constructing acidFBA model (universal and homogeneous amino acid distribution) ...');
univBool = true;
[acidFBAmodelUniv,acidFBAuniprotUniv] = createAcidFBAmodel(geckoPool,univBool);
fprintf('DONE \n');
[aaMatUniv,aaMinUniv,aaMaxUniv] = sampleNutrientEnvironments(acidFBAmodelUniv,nSamples,'Universal');

