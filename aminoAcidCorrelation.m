function aminoAcidCorrelation(model,protArray)
%aminoAcidCorrelation Plots the acidFBA model-predicted amino acid and protein
%distribution against experimental measurements.
%   This functions calculates and plots the model-predicted amino acid
%   and protein distribution against experimentally measurements from 
%   Bartolomeo et al. 2020 (exponential, glucose-limited growth).
%
% USAGE:
%
%    aminoAcidCorrelation(model,protArray)
%
% INPUTS:
%    model:                    AcidFBA model structure
%
%    protArray:                Cell array containing the identifiers
%                              and sequences of the GECKO-implemented proteins.
%                              Has the following columns
%                                   1. Identifiers
%                                   2. Amino acid sequences (single letter codes)
%
% .. Authors:
%       - Vetle Simensen 27/01/21

%% Initialize
params = getParameters();
sigma = 0.439;     % numerical issues with sigma = 0.44

% UniProt IDs and amino acid sequences
uniprotSeqHash = containers.Map(protArray(:,1),protArray(:,2));

%% Predicted amino acid mass fractions 
fprintf('Simulating amino acid usage ...');
simModel = constrainProtPool(model,params.Ptot,params.f,sigma);
sol = optimizeCbModel(simModel,'max','one');  % optimize growth, minimize sum of fluxes
simAAmassFrac = sol.v(findRxnIDs(simModel,params.aaDrains))';     % amino acid mass fraction usage
protIdxs = ~cellfun(@isempty,regexp(model.rxns,'prot_\w*_drain','match'));
geckoProts = cellfun(@(x) regexprep(x,{'prot_','_drain'},''),model.rxns(protIdxs),'UniformOutput',false);
simProtMoleFrac = sol.v(protIdxs)';     % protein mole fraction usage
simProtMoleFracHash = containers.Map(geckoProts,simProtMoleFrac);
[aaMin,aaMax] = fluxVariability(simModel,99,'max',params.aaDrains);
fprintf('DONE \n');

%% Process experimental proteomics measurements
fprintf('Reading proteomics data ...');
currDir = cd;
cd proteins
expDat = readcell('Bartolomeo2020_proteomics.xlsx');
cd(currDir);

% Calculate means and SDs from three biological replicates
protReplicates = cell2mat(expDat(3:end,5:7));
protMeans = mean(protReplicates,2);
protSDs = std(protReplicates,0,2);
expDat(3:end,5) = num2cell(protMeans);
expDat(3:end,6) = num2cell(protSDs);
expDat(2,5:6) = {'Mean (g/gDW)','SD (g/gDW)'};
expDat(1,:) = [];
expDat(:,[2 3 7:18]) = [];
fprintf('DONE \n');

%% Measured amino acid mass fractions
fprintf('Calculating experimental amino acid usage ...');
idxs = regexp(model.rxns,'prot_\w*_drain');     % find GECKO proteins (drain reactions)
idxs(cellfun(@isempty, idxs)) = {0};
protDrains = model.rxns(logical(cell2mat(idxs)));
geckoProt = cellfun(@(x) regexprep(x,{'prot_','_drain'},''),protDrains,'UniformOutput',false);
nProt = length(geckoProt);
qAAmassFrac = zeros(1,20);    % g/gDW of each amino acid
unmeasuredProt = 0;     % flux-carrying, but unmeasured protein mass

% Calculate mass fractions (g/gDW) of amino acids from experimental
% measurements
for i = 1:nProt
    uniprotID = geckoProt{i};
    aaSeq = uniprotSeqHash(uniprotID);
    [aaMassFrac,protMW] = findAaMassFrac(aaSeq);
    idx = find(ismember(expDat(:,1),uniprotID),1);
    
    % Add mass contribution of each amino acid (g/gDW)
    if ~isempty(idx)
        protLevel = expDat{idx,3};
        qAAmassFrac = qAAmassFrac + (aaMassFrac * protLevel);     % amino acid mass fraction (g/gDW)
    else
        unmeasuredProt = unmeasuredProt + (simProtMoleFracHash(uniprotID) * protMW * 1/1000);
    end
    
end
fprintf('DONE \n');

% Print overall protein usage and mass of active and unmeasured GECKO proteins
qAAmassFrac = qAAmassFrac * sigma;
fprintf('Experimental metabolic protein mass (g/gDW): %5.4f\n',sum(qAAmassFrac));
fprintf('Simulated metabolic protein mass (g/gDW): %5.4f\n',sum(simAAmassFrac));
fprintf('Active and unmeasured GECKO protein mass (g/gDW): %5.4f\n',unmeasuredProt);

% Correlation and mean relative error
[rho,pval] = corrcoef(qAAmassFrac,simAAmassFrac);
relError = mean(abs((simAAmassFrac - qAAmassFrac) ./ qAAmassFrac));
fprintf('Correlation of %4.3f (p-val = %d) and a mean relative error of %3.2f\n',round(rho(1,2),3),pval(1,2),relError);

% Save data
t = table(qAAmassFrac', simAAmassFrac', aaMax, aaMin,'VariableNames',{'qAAmassFrac','simAAmassFrac','aaMax','aaMin'});
writetable(t,'data/aaCorrDat.csv');

end

function [aaMassFrac,protMW] = findAaMassFrac(seq)
%findAaMassFrac Calculates the mass fraction of all amino acids in seq as
%gram amino acid per gram protein (g/g), as well as the molecular mass
%(g/mol) of the protein
%
% INPUTS:
%    seq:                      Protein sequence (char array)
%
% OUTPUTS:                     
%    aaMassFrac:               Amino acid mass fractions (g/g)
%
%    protMW:                   Protein molecular weight (g/mol)

params = getParameters();
nTerm = 2 * 1.007825;   % mass of 2 H
cTerm = 15.999405;  % mass of 1 O
aaStart = seq(1);   % N-terminal amino acid
aaEnd = seq(end);   % C-terminal amino acid
startIdx = find(ismember(params.aaCodes,aaStart));
endIdx = find(ismember(params.aaCodes,aaEnd));

% Calculate protein molecular weight (g/mol)
aaMoles = cellfun(@(x) length(strfind(seq,x)),params.aaCodes);
aaMasses = aaMoles .* params.aaMWprot;
aaMasses(startIdx) = aaMasses(startIdx) + nTerm;
aaMasses(endIdx) = aaMasses(endIdx) + cTerm;
protMW = sum(aaMasses);
aaMassFrac = aaMasses ./ protMW;    % amino acid mass fractions (g/g)

end
