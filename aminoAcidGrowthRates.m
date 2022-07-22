function aminoAcidGrowthRates(model)
%aminoAcidGrowthRates Calculates the amino acid utilization against varying
%growth rates. Simulates the data presented in Fig. 2 of the manuscript.
%
% USAGE:
%
%    aminoAcidGrowthRates(model)
%
% INPUTS:
%    model:             AcidFBA model structure
%
% .. Authors:
%       - Vetle Simensen 22/11/21

%% Initialize
params = getParameters;
aaDrains = params.aaDrains;
aaDrainsRxnIdxs = findRxnIDs(model, aaDrains);
growthRates = round(linspace(0.0333, 0.4, 12),4);     % growth rates between 0.0333 and 0.4

%% Non-linear regression
% Fit a non-linear regression model to growth-rate dependent protein fractions of
% biomass from D. Pejin and R. Razmovski (1993).
dilRates = [0.05 0.07 0.10 0.20 0.42];
protFracs = [48.8 51.6 52.7 53.9 54.2] / 100;
modelFun = @(b, x) b(1) .* (1 - exp(-b(2) .* x));   % scaled exponential curve
initCoeff = [54; 50];
mdl = fitnlm(dilRates, protFracs, modelFun, initCoeff);
growthProteins = modelFun(mdl.Coefficients.Estimate, growthRates);  % corresponding protein fractions

%% Simulate metabolic phenotype
aaUsage = zeros(length(growthRates), length(aaDrains));     % data for scatter plots

for i = 1:length(growthRates)
    % Constrain growth and protein availability
    tempModel = changeRxnBounds(model, model.rxns(logical(model.c)), growthRates(i), 'b');
    tempModel = constrainProtPool(tempModel, growthProteins(i), params.f, 0.46);
    
    % Block L-serine transport from mitochondria to cytoplasm
    tempModel = changeRxnBounds(tempModel, 'r_2045_REV', 0, 'u');
    
    % Constrain secretion of pyruvate, (R,R)-2,3-butanediol, acetaldehyde,
    % and glycine to physiological levels
    tempModel = changeRxnBounds(tempModel, {'r_2033', 'r_1549', 'r_1631', 'r_1810'}, 1e-5, 'u');
    
    % Minimize substrate uptake (glucose) and constrain to minimum (+/-
    % 0.1% flexibility)
    tempModel = changeObjective(tempModel, 'r_1714_REV');
    solSub = optimizeCbModel(tempModel, 'min');
    tempModel = changeRxnBounds(tempModel, 'r_1714_REV', solSub.x(findRxnIDs(tempModel, 'r_1714_REV')) * 0.999, 'l');
    tempModel = changeRxnBounds(tempModel, 'r_1714_REV', solSub.x(findRxnIDs(tempModel, 'r_1714_REV')) * 1.001, 'u');
    
    % Minimize overall enzyme usage
    tempModel = changeObjective(tempModel, 'prot_pool_exchange');
    sol = optimizeCbModel(tempModel, 'min');
    
    aaUsage(i,:) = sol.x(aaDrainsRxnIdxs);  % amino acid drain fluxes
end

% Save data
writematrix(aaUsage,'data/aaUsage.csv');
writematrix(growthRates,'data/growthRates.csv');

%% Simulate respiratory vs. fermentative amino acid distribution
% Aerobic model
aerobModel = model;
aerobModel = changeRxnBounds(aerobModel,aerobModel.rxns(logical(aerobModel.c)),0.2,'u');    % limit maximal growth rate to 0.2 (prevent respiro-fermentative metabolism)
constrainProtPool(aerobModel,params.Ptot,params.f,0.44);

% Block L-serine transport from mitochondria to cytoplasm
aerobModel = changeRxnBounds(aerobModel, 'r_2045_REV', 0, 'u');
    
% Constrain secretion of pyruvate, (R,R)-2,3-butanediol, acetaldehyde,
% and glycine to physiological levels
aerobModel = changeRxnBounds(aerobModel, {'r_2033', 'r_1549', 'r_1631', 'r_1810'}, 1e-5, 'u');

% Anaerobic model
anaerobModel = convertToAnaerob(aerobModel);

% Minimize total sum of fluxes
solAerob = optimizeCbModel(aerobModel,'max','one');
solAnaerob = optimizeCbModel(anaerobModel,'max','one');

% Relative amino acid distribution
aerobAA = solAerob.x(aaDrainsRxnIdxs);
anaerobAA = solAnaerob.x(aaDrainsRxnIdxs);
aerobAA = aerobAA ./ sum(aerobAA);
anaerobAA = anaerobAA ./ sum(anaerobAA);
aaBar = abs(1 - anaerobAA ./ aerobAA);

% Sort, and reorder amino acid codes
[Y,I] = sort(aaBar,'descend');
X = categorical(params.aaThreeCodes(I));
X = reordercats(X,params.aaThreeCodes(I));

% Save data
respFermTable = table(X',Y,'VariableNames',{'Amino acids','AbsRelDiff'});
writetable(respFermTable,'data/absRelDiff.csv');

end
