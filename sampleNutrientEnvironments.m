function [aaMat,aaMinMat,aaMaxMat] = sampleNutrientEnvironments(model,nSamples,fileNameSuffix)
%sampleNutrientEnvironments Simulates growth on multiple environmental
%conditions and predicts the amino acid composition of the metabolic
%proteome.
%
% USAGE:
%
%    environmentalGrowth(model)
%
% INPUTS:
%    model:             AcidFBA model structure
%
%    nSamples:          Number of nutrient combinations
%
%    fileNameSuffix:    Suffix of the filename to store the output data
%
% OUTPUTS:
%    aaMat:             nSamples x 20 matrix of amino acid distributions
%
%    aaMinMat:          nSamples x 20 matrix of min amino acid flux
%                       (FVA)
%
%    aaMaxMat:          nSamples x 20 matrix of max amino acid flux
%                       (FVA)
%
% .. Authors:
%       - Vetle Simensen 08/09/21

%% Initialize
rng('default');     % set random seed for reproducibility
params = getParameters;
aaDrains = params.aaDrains;
aaRxnIds = findRxnIDs(model,params.aaDrains);

% Constrain protein availability
protAmount = 0.5;
f = 0.4461;
sigma = 0.44;
model = constrainProtPool(model,protAmount,f,sigma);

% Wild type FBA
wtFBA = optimizeCbModel(model);
wtGrowth = wtFBA.f;

%% Identify viable carbon, nitrogen, phosphorus, and sulfur sources
% By default, the model uses glucose, ammonia, phosphate, and
% sulphate. Test the remaining exchange reactions for viable 
% candidates for each elemental nutrient source.
fprintf('Identifying viable nutrient sources... ');
dfltExc = {'r_1714_REV','r_1654_REV','r_2005_REV','r_2060_REV'};
elements = {'C','N','P','S'};
excRxns = model.rxns(findExcRxns(model));
excRxns = excRxns(contains(excRxns,'_REV'));   % only uptake rxns
excRxns(ismember(excRxns,'r_2111_REV')) = [];   % remove biomass exchange
nExcRxns = length(excRxns);
eleSrcs = zeros(nExcRxns,4);

environment = getEnvironment();
parfor i = 1:nExcRxns
    restoreEnvironment(environment);
    excRxn = excRxns{i,1};
    excMet = findMetsFromRxns(model,excRxn);
    excMetIdx = findMetIDs(model,excMet);
    metFormula = model.metFormulas{excMetIdx};
    for j = 1:4
        if contains(metFormula,elements{j})
            growthBool = testNutrientSource(model,excRxn,dfltExc{j},j);
            eleSrcs(i,j) = growthBool;
        end
    end
end

% All nutrient combinations
nutrComb = {excRxns(logical(eleSrcs(:,1))),excRxns(logical(eleSrcs(:,2))),excRxns(logical(eleSrcs(:,3))),excRxns(logical(eleSrcs(:,4)))};
n = length(nutrComb);
[nutrComb{:}] = ndgrid(nutrComb{:});
nutrComb = reshape(cat(n + 1,nutrComb{:}),[],n);
nutrComb = nutrComb(randi(length(nutrComb),round(nSamples * 1.5),1),:);    % 50% as many combinations in case of infeasible solutions
nNutrComb = length(nutrComb);
fprintf('DONE\n');

%% Calculate amino acid usage when changing only carbon source
fprintf('Simulating optimal growth on various nutrient sources...');
aaMinMat = zeros(nNutrComb,20);
aaMaxMat = zeros(nNutrComb,20);
aaMat = zeros(nNutrComb,20);

environment = getEnvironment();
parfor i = 1:nNutrComb
    restoreEnvironment(environment);
    simModel = model;
    simModel = changeRxnBounds(simModel,dfltExc,0,'u');   % remove default nutrients
    simModel = changeRxnBounds(simModel,nutrComb(i,:),1000,'u');  % add new nutrient sources
    sol = optimizeCbModel(simModel,'max');

    if sol.f > wtGrowth * 0.01 && sol.stat == 1
        try
            [aaMin,aaMax] = fluxVariability(simModel,99,'max',aaDrains);
        catch
            aaMin = zeros(1,length(aaDrains));
            aaMax = zeros(1,length(aaDrains));
        end

        aaMat(i,:) = sol.x(aaRxnIds);
        aaMinMat(i,:) = aaMin;
        aaMaxMat(i,:) = aaMax;
    end
end
fprintf('DONE\n');

% Remove superfluous/infeasible samples
aaMat(all(aaMat == 0,2),:) = [];
aaMinMat(all(aaMinMat == 0,2),:) = [];
aaMaxMat(all(aaMaxMat == 0,2),:) = [];

if length(aaMat) > nSamples
    aaMat(nSamples + 1:end,:) = [];
    aaMinMat(nSamples + 1:end,:) = [];
    aaMaxMat(nSamples + 1:end,:) = [];
end

% Save files
writematrix(aaMat,['data/aaMat' fileNameSuffix '.csv']);
writematrix(aaMinMat,['data/aaMinMat' fileNameSuffix '.csv']);
writematrix(aaMaxMat,['data/aaMaxMat' fileNameSuffix '.csv']);

end

function growthBool = testNutrientSource(model,newExc,oldExc,j)
    growthBool = 0;
    model = changeRxnBounds(model,oldExc,0,'u');
    model = changeRxnBounds(model,newExc,1000,'u');
    try
        sol = optimizeCbModel(model);
        if sol.f > 0 && sol.stat == 1
            growthBool = 1;
        end
    catch
        fprintf('Failed (%s)\n',str(j));
    end
end
