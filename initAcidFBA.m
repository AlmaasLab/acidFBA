%% Initializes the construction and analysis of the yeast acidFBA model
% Script for constructing and analyzing an acidFBA-implementation of the
% yeast GEM, Yeast 8.3. All generated output files are stored in the 
% /data folder.
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

% Add metabolic subsystems
fprintf('Adding metabolic subsystems ...');
geckoPool = addMetabolicSubsystems(geckoPool);
fprintf('DONE \n');

% Construct the acidFBA model
fprintf('Constructing acidFBA model ...');
[acidFBAmodel,protArray] = createAcidFBAmodel(geckoPool);
fprintf('DONE \n');


%% Correlation with experimental data (Fig. 1 and S1 Fig.)
aminoAcidCorrelation(acidFBAmodel,protArray);
% Output file aaCorrDat.csv, with the following headers:
%   qAAmassFrac:   experimental amino acid mass fractions (g/gDW)
%   simAAmassFrac: simulated amino acid mass fractions (g/gDW)
%   aaMax:         max simulated amino acid mass fraction (g/gDW) at 99% of optimal
%                  growth (FVA)
%   aaMin:         min simulated amino acid mass fraction (g/gDW) at 99% of optimal
%                  growth (FVA)


%% Growth rate-dependency of amino acid usage (Fig. 2)
aminoAcidGrowthRates(acidFBAmodel);
% Output file aaUsage.csv, containing the simulated amino acid mass fractions (g/gDW)
% for every amino acid (column) at each growth rate (row)

% Output file growthRates.csv, containing the corresponding growth rates
% (1/h).

% Output file absRelDiff.csv, with the following headers:
%   Amino acids:    three letter amino acid codes
%   AbsRelDiff:     absolute, relative mass fraction difference between a
%                   fully fermentative versus a fully respiratory metabolism


%% Sample nutrient environments (Fig. 3, 4)
nSamples = 5000;
[aaMatDflt,aaMinDflt,aaMaxDflt] = sampleNutrientEnvironments(acidFBAmodel,nSamples,'Default');
% Output file aaMatDefault.csv, containing the simulated amino acid (column) mass
% fractions (g/gDW) for the <nSamples> sampled conditions (row)

% Output file aaMinMatDefault.csv, containing the min simulated amino acid
% mass fraction (g/gDW) at 99% of optimal growth (FVA). Same format as
% aaMatDefault.csv.

% Output file aaMaxMatDefault.csv, containing the max simulated amino acid
% mass fraction (g/gDW) at 99% of optimal growth (FVA). Same format as
% aaMatDefault.csv.


%% Species-specific amino acid distributions (Fig. 5) 
minAminoAcidDiff(acidFBAmodel);
% Output file orgRelGrowthRates.csv, containing the relative growth rates.

% Output file aaSce.csv, aaEco.csv, aaBsu.csv, and aaPse.csv, containing
% the absolute relative difference of each amino acid between the
% the simulated amino acid profile and the species-specific amino acid 
% distributions. 

%% Sample nutrient environments, universal amino acid distribution (S2 Fig)
fprintf('Constructing acidFBA model (universal and homogeneous amino acid distribution) ...');
univBool = true;
[acidFBAmodelUniv,protArrayUniv] = createAcidFBAmodel(geckoPool,univBool);
fprintf('DONE \n');

[aaMatUniv,aaMinUniv,aaMaxUniv] = sampleNutrientEnvironments(acidFBAmodelUniv,nSamples,'Universal');
% Output files aaMatUniversal.csv, aaMinMatUniversal.csv, and 
% aaMaxMatUniversal.csv with the same formats as described above. 

