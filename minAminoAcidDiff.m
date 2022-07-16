function minAminoAcidDiff(model)
%minAminoAcidDiff Minimizes the Euclidean distance to the species-specific
%amino acid distributions across a range of growth rates. Calculates the
%absolute relative difference for each amino acid.
%
% USAGE:
%
%    minAminoAcidDiff(model)
%
% INPUTS:
%    model:             AcidFBA model structure
%
% .. Authors:
%       - Vetle Simensen 29/06/22

% Initialize
params = getParameters();
changeCobraSolver('gurobi','all');
changeCobraSolverParams('QP','feasTol',1.0000e-09);
changeCobraSolverParams('QP','optTol',1.0000e-09);

% Get amino acid mass fractions (S. cerevisiae, E. coli, B. subtilis)
currDir = cd;
cd data
aaFrac = readmatrix('aaMassFracsOrganisms.xlsx');
aaFrac = aaFrac(:,2:5);     % fix for more organisms
[nAA,nOrg] = size(aaFrac);
cd(currDir);

% Wild type yeast reference conditions (experimental amino acid profile)
protIdx = findRxnIDs(model,'prot_pool_exchange');
aaYeastProfile = aaFrac(:,1);
wtModel = changeRxnBounds(model,params.aaDrains,aaYeastProfile * model.ub(protIdx),'b');
wtSol = optimizeCbModel(wtModel);

% Define range of relative growth rates
nVals = 45;
relGrowthRates = linspace(0.1,2.5,nVals);
yvals = zeros(nAA,nVals,nOrg);

% Each organism
for i = 1:nOrg
    % Select experimental organism-specific amino acid distribution
    aaProfile = model.ub(protIdx) * aaFrac(:,i);
    
    % Each fraction of optimal yeast growth
    for j = 1:length(relGrowthRates)
        simModel = model;
        
        % Set fraction of optimal growth
        simModel = changeRxnBounds(simModel,simModel.rxns(logical(simModel.c)),relGrowthRates(j) * wtSol.f,'b');
        
        % Construct QP problem
        QPproblem = buildLPproblemFromModel(simModel);
        [~, n] = size(simModel.S);
        objIndxs = findRxnIDs(model,params.aaDrains);
        
        % Formulate quadratic objective term
        quadMat = sparse(n,n);
        objIndxsVec = zeros(n,1);
        objIndxsVec(objIndxs) = 1;
        quadMat = quadMat + diag(objIndxsVec);
        QPproblem.F = quadMat;

        % Formulate linear objective term
        linVec = zeros(n,1);
        linVec(objIndxs) = -aaProfile;
        QPproblem.c = linVec;
        
        % Add linear constraints
        QPproblem.A = simModel.S;  % LHS matrix
        QPproblem.b = simModel.b;  % RHS vector
        QPproblem.lb = simModel.lb;    % lower flux bounds
        QPproblem.ub = simModel.ub;    % upper flux bounds
        QPproblem.osense = 1;   % objective sense linear part (minimize)
        QPproblem.csense = buildLPproblemFromModel(simModel).csense;   % constraint senses (all 'E')

        % Solve QP problem
        momaSol = solveCobraQP(QPproblem);
        
        % Store solution
        yvals(:,j,i) = abs(((momaSol.full(objIndxs) ./ momaSol.full(protIdx)) - aaFrac(:,i)) ./ (aaFrac(:,i)));
    end
end

% Save data
writematrix(relGrowthRates,'data/orgRelGrowthRates.csv');   % relative growth rates
fileSuffixes = {'Sce' 'Eco' 'Bsu' 'Pse'};

for i = 1:nOrg
    writematrix(yvals(:,:,i),['data/aa' fileSuffixes{i} '.csv']);
end

end

