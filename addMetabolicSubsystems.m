function newModel = addMetabolicSubsystems(model)
%v Queries the KEGG REST API for metabolic pathway
%informatino on S. cerevisiae.
%
% USAGE:
%
%    newModel = addMetabolicSubsystems(model)
%
% INPUTS:
%    model:             Model structure
%
% OUTPUTS:
%    newModel:          Updated model structure with metabolic subsystem
%                       information.
%
% .. Authors:
%       - Vetle Simensen 20/07/22

newModel = model;

% Query KEGG REST API for S. cerevisiae pathways
options = weboptions;
options.CertificateFilename = ('');
keggCharArray = split(webread('http://rest.kegg.jp/list/pathway/sce', options));
pathwayIdxs = regexp(keggCharArray,'path:');
pathwayIdxs(cellfun(@isempty, pathwayIdxs)) = {0};
scePathways = keggCharArray(logical(cell2mat(pathwayIdxs)));
scePathways = strrep(scePathways,'path:','');

% Read yeast-GEM.xlsx-file for subsystem annotations
currDir = cd;
cd models
yeastAnnotations = readcell('yeast-GEM.xlsx','Sheet',1);
yeastAnnotations(cellfun(@(x) isa(x,'missing'), yeastAnnotations)) = {''};
cd(currDir);

% Add subsystem information to model structure
for i = 2:length(yeastAnnotations)
    rxn = yeastAnnotations{i,2};
    annotation = yeastAnnotations{i,11};
    modelIdxs = find(contains(newModel.rxns,rxn));
    
    if ~isempty(modelIdxs)
        pathIdxs = regexp(annotation,scePathways);
        pathIdxs(cellfun(@isempty, pathIdxs)) = {0};
        subsystems = scePathways(logical(cell2mat(pathIdxs)));
        newModel = addSubSystemsToReactions(newModel,newModel.rxns(modelIdxs),subsystems);
    end
end

end

