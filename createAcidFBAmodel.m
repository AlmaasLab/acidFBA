function [acidFBAmodel,protArray] = createAcidFBAmodel(model,univBool)
%createAcidFBAmodel Converts a GECKO model to an acidFBA model.
%   Creates an acidFBA-GEM from the GECKO model structure. It begins by
%   constructing the Xi matrix in which every element Xi_(i,k) denotes the 
%   molecular weight of amino acid k with respect to protein i. The
%   individual protein drain reactions are then removed and replaced by the
%   corresponding amino acid drain reactions. Finally, amino acid source
%   reactions are added for all 20 amino acids which individually draws from
%   the protein pool reaction.
%
% USAGE:
%
%    [acidFBAmodel,protArray] = createAcidFBAmodel(model,univBool)
%
% INPUTS:
%    model:             GECKO model structure (protein pool)
%
%    univBool:          True if a universal and homogeneous amino acid
%                       distribution is to be used, false otherwise 
%                       (default).
%
% OUTPUTS:
%    acidFBAmodel:      acidFBA model struct
%
%    protArray:         Filtered cell array containing the identifiers
%                       and sequences of the GECKO-implemented proteins.
%                       Has the following columns
%                           1. Identifiers
%                           2. Amino acid sequences (single letter codes)
%
% .. Authors:
%       - Vetle Simensen 18/08/20

params = getParameters;
acidFBAmodel = model;

if ~exist('universalBool','var')
    univBool = false;
end

%% Construct amino acid composition matrix Xi
% Amino acid IUPAC codes and molecular weights (protein-bound and free)
aaCodes = params.aaCodes;
aaMWprot = params.aaMWprot;

% Read protein sequence data
currDir = cd;
cd proteins
uniprotArray = readcell('protSeqs.xlsx');
cd(currDir);

% Construct Xi matrix - molecular weight of amino acid l with respect to protein k
nProteins = length(uniprotArray);
nAminoAcids = length(aaCodes);
xi = zeros(nProteins,nAminoAcids);
nTerm = 2 * 1.007825;   % mass of 2 H
cTerm = 15.999405;  % mass of 1 O

% Universal amino acid distribution
if univBool
    aaFracs = [0.0495, 0.0103, 0.0606, 0.0678, 0.0575, 0.0362, 0.0272, 0.0694, 0.0847, 0.0917, 0.0245, 0.0463, 0.0385, 0.0361, 0.0592, 0.0540, 0.0521, 0.0675, 0.0181, 0.0486];
    aaFracs = aaFracs ./ sum(aaFracs);
end

for k = 2:nProteins
    sequence = uniprotArray{k,2};
    aaStart = sequence(1);   % N-terminal amino acid
    aaEnd = sequence(end);   % C-terminal amino acid
    
    for l = 1:nAminoAcids
        addTerms = [0 0];
        
        if strcmp(aaCodes{l},aaStart)
            addTerms(1) = 1;
        end
        
        if strcmp(aaCodes{l},aaEnd)
            addTerms(2) = 1;
        end
        
        % Universal and homogeneous amino acid distribution
        if univBool
            aaMole = length(sequence) * aaFracs(l);
        else
            aaMole = length(strfind(sequence,aaCodes{l}));
        end
        
        % Add to Xi matrix
        xi(k - 1,l) = ((aaMole * aaMWprot(l)) + (addTerms(1) * nTerm) + (addTerms(2) * cTerm)) * 1/1000;    % convert to [g/mmol] using scaling factor
    end
end

%% Extract information and sequences for the GECKO-implemented proteins
% Identify the set of proteins that are incorporated into the GECKO formulation
rxnPrefix = 'draw_prot_';
protIdxs = startsWith(acidFBAmodel.rxns,rxnPrefix);
protExcRxns = acidFBAmodel.rxns(protIdxs);
geckoProteins = cellfun(@(x) regexprep(x,['^' rxnPrefix],''),protExcRxns,'UniformOutput',false);
[~,geneIdx] = intersect(uniprotArray(:,1),char(geckoProteins));
protArray = uniprotArray(geneIdx,:);    % exract protein sequences for GECKO proteins only
xi = xi(geneIdx - 1,:);     % do the same for the Xi matrix


%% Modify GECKO model
% Remove the individual protein drain reactions
rxns = acidFBAmodel.rxns;
aaMets = cellfun(@(x) [x '[c]'],aaCodes,'UniformOutput',false);
acidFBAmodel = addMultipleMetabolites(acidFBAmodel,aaMets);
protIdxs = find(protIdxs);
for k = 1:length(protIdxs)
    enzymeMet = findMetsFromRxns(acidFBAmodel,rxns{protIdxs(k)});  % enzyme k
    uniprotID = regexp(enzymeMet{1},'_\w{6}','match');
    uniprotID = strrep(uniprotID{1},'_','');
    geneIdx = ismember(protArray(:,1),uniprotID);
    acidFBAmodel = removeRxns(acidFBAmodel,rxns{protIdxs(k)}); % remove old reaction
    
    % Add amino acid pool reactions
    mets = [aaMets enzymeMet{1}];
    stoich = [-xi(geneIdx,:) 1];
    acidFBAmodel = addReaction(acidFBAmodel,['prot_' uniprotID '_drain'],'metaboliteList',mets,...
        'stoichCoeffList',stoich,'reactionName',...
        ['prot_' uniprotID '_drain'],'reversible',false,'geneRule','');
end

% Add amino acid source reactions
rxnIDs = cellfun(@(x) ['AAR_pool_' x],aaCodes,'UniformOutput',false);
metList = ['prot_pool[c]' aaMets];
stoich = [(-1 * ones(20,1)) eye(20)]';
acidFBAmodel = addMultipleReactions(acidFBAmodel,rxnIDs,metList,stoich,'lb',zeros(length(rxnIDs),1),'ub',ones(length(rxnIDs),1) * 1000);

end

