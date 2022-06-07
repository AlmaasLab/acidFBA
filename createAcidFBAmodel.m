function [acidFBAmodel,uniProtTable] = createAcidFBAmodel(model,univBool)
%createAcidFBAmodel Converts a GECKO model to an acidFBA model.
%   Creates an acidFBA GEM from the GECKO model struct model. It begins by
%   constructing the Xi matrix in which every element Xi_(i,k) denotes the 
%   molecular weight of amino acid k with respect to protein i. The
%   individual protein drain reactions are then removed and replaced by the
%   corresponding amino acid drain reactions. Finally, amino acid source
%   reactions are added for all 20 amino acids which individually draws from
%   the protein pool reaction.
%
% USAGE:
%
%    [acidFBAmodel,acidFBAuniprot] = createAcidFBAmodel(model,univBool)
%
% INPUTS:
%    model:             GECKO model struct (protein pool)
%
%    univBool:          True if a universal and homogeneous amino acid
%                       distribution is to be used, false otherwise 
%                       (default).
%
% OUTPUTS:
%    acidFBAmodel:      acidFBA model struct
%
%    uniProtTable:      Filtered uniprotArray containing relevant 
%                       information on GECKO-implemented proteins only
%                       with the following columns:
%                           1. UniProt ID
%                           2. Sequence length
%                           3. Yeast gene name
%                           4. Amino acid sequence
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

% Read UniProt data
currDir = cd;
cd proteins
uniprotArray = readcell('uniprotYeast.xlsx');
uniprotArray(cellfun(@(x) isa(x,'missing'), uniprotArray)) = {''};
cd(currDir);

% Construct Xi matrix - molecular weight of amino acid l with respect to protein k
nProteins = length(uniprotArray);
nAminoAcids = length(aaCodes);
xi = zeros(nProteins,nAminoAcids);
nTerm = 2 * 1.007825;   % mass of 2 H
cTerm = 15.999405;  % mass of 1 O

% Universal and homogeneous amino acid distribution
if univBool
    aaFracs = [0.0495, 0.0103, 0.0606, 0.0678, 0.0575, 0.0362, 0.0272, 0.0694, 0.0847, 0.0917, 0.0245, 0.0463, 0.0385, 0.0361, 0.0592, 0.0540, 0.0521, 0.0675, 0.0181, 0.0486];
    aaFracs = aaFracs ./ sum(aaFracs);
end

for k = 2:nProteins
    sequence = uniprotArray{k,8};
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


%% Modify GECKO model
mets = acidFBAmodel.mets;
metNames = acidFBAmodel.metNames;
comps = acidFBAmodel.comps;
compNames = acidFBAmodel.compNames;
rxnPrefix = 'draw_prot_';
colIdx = 1;

% Add compartment suffixes to metNames
compsReg = cellfun(@(x) ['\[' x '\]$'],comps,'UniformOutput',false);

for i = 1:length(mets)
    compIdx = ~cellfun(@isempty,regexp(mets{i},compsReg,'match'));
    metNames{i} = [metNames{i} ' [' compNames{compIdx} ']'];
end

% Swap met names and met IDs
acidFBAmodel.mets = metNames;
acidFBAmodel.metNames = mets;

%% Extract information and sequences for the GECKO-implemented proteins
% Identify the set of genes that are incorporated into the GECKO formulation
protIdxs = startsWith(acidFBAmodel.rxns,rxnPrefix);
protExcRxns = acidFBAmodel.rxns(protIdxs);
geckoProteins = cellfun(@(x) regexprep(x,['^' rxnPrefix],''),protExcRxns,'UniformOutput',false);
[~,geneIdx] = intersect(uniprotArray(:,colIdx),char(geckoProteins));
uniProtTable = uniprotArray(geneIdx,:);    % exract the UniProt table for GECKO proteins only
uniProtTable(:,[2 3 4 7]) = [];    % delete irrelevant data
xi = xi(geneIdx - 1,:);     % do the same for the Xi matrix

%% Add subsystems
% Query KEGG REST API for S. cerevisiae pathways
keggCharArray = split(webread('http://rest.kegg.jp/list/pathway/sce'));
pathwayIdxs = regexp(keggCharArray,'path:');
pathwayIdxs(cellfun(@isempty, pathwayIdxs)) = {0};
scePathways = keggCharArray(logical(cell2mat(pathwayIdxs)));
scePathways = strrep(scePathways,'path:','');

% Read yeast-GEM .xlsx-file for subsystem annotations
currDir = cd;
cd models
yeastAnnotations = readcell('yeast-GEM.xlsx','Sheet',1);
yeastAnnotations(cellfun(@(x) isa(x,'missing'), yeastAnnotations)) = {''};
cd(currDir);

% Add subsystem information to model struct
for i = 2:length(yeastAnnotations)
    rxn = yeastAnnotations{i,2};
    annotation = yeastAnnotations{i,11};
    modelIdxs = find(contains(acidFBAmodel.rxns,rxn));
    
    if ~isempty(modelIdxs)
        pathIdxs = regexp(annotation,scePathways);
        pathIdxs(cellfun(@isempty, pathIdxs)) = {0};
        subsystems = scePathways(logical(cell2mat(pathIdxs)));
        acidFBAmodel = addSubSystemsToReactions(acidFBAmodel,acidFBAmodel.rxns(modelIdxs),subsystems);
    end
end

%% Modify GECKO model
% Remove the individual protein drain reactions
rxns = acidFBAmodel.rxns;
aaMets = cellfun(@(x) [x '[cytoplasm]'],aaCodes,'UniformOutput',false);
acidFBAmodel = addMultipleMetabolites(acidFBAmodel,aaMets);
protIdxs = find(protIdxs);
for k = 1:length(protIdxs)
    enzymeMet = findMetsFromRxns(acidFBAmodel,rxns{protIdxs(k)});  % enzyme k
    uniprotID = regexp(enzymeMet{1},'_\w{6}','match');
    uniprotID = strrep(uniprotID{1},'_','');
    geneIdx = ismember(uniProtTable(:,1),uniprotID);
    acidFBAmodel = removeRxns(acidFBAmodel,rxns{protIdxs(k)}); % remove old reaction
    
    % Add the amino acid pool reactions
    mets = [aaMets enzymeMet{1}];
    stoich = [-xi(geneIdx,:) 1];
    acidFBAmodel = addReaction(acidFBAmodel,['prot_' uniprotID '_drain'],'metaboliteList',mets,...
        'stoichCoeffList',stoich,'reactionName',...
        ['prot_' uniprotID '_drain'],'reversible',false,'geneRule','');
end

% Add amino acid source reactions
rxnIDs = cellfun(@(x) ['AAR_pool_' x],aaCodes,'UniformOutput',false);
metList = ['prot_pool [cytoplasm]' aaMets];
stoich = [(-1 * ones(20,1)) eye(20)]';
acidFBAmodel = addMultipleReactions(acidFBAmodel,rxnIDs,metList,stoich,'lb',zeros(length(rxnIDs),1),'ub',ones(length(rxnIDs),1) * 1000);

end

