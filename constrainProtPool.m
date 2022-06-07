function model = constrainProtPool(model,protAmount,f,sigma)
%constrainProtPool Constrains the protein pool reaction according to the
%given amount.
%   Detailed explanation goes here
%
% USAGE:
%
%    modelOut = constrainProtPool(model,protAmount,f,sigma)
%
% INPUTS:
%    model:             GECKO/acidFBA model struct (protein pool)
%    protAmount:        Total protein dry weight fraction (g/gDW) (double)
%
% OPTIONAL IPUTS:
%    f:                 Mass fraction of proteins accounted for in the model
%    sigma:             Average enzyme saturation factor
%   
% OUTPUTS:
%    model:             GECKO/acidFBA model struct with constrained 
%                       protein pool reaction
%
%
% .. Authors:
%       - Vetle Simensen 25/08/20

params = getParameters();

if ~exist('f','var')
    f = params.f;
end

if ~exist('sigma','var')
    sigma = params.sigma;
end

if strcmp(model.id,'emodel_Saccharomyces_cerevisiae')
    protPoolRxn = params.protPool2;
    fluxBound = -0.2300 * protAmount / 0.5;     % default lower bound of -0.2300
    model = changeRxnBounds(model,protPoolRxn,fluxBound,'l');
elseif strcmp(model.id,'ecYeastGEM_batch_v8.3.4')
    protPoolRxn = params.protPool;
    fluxBound = protAmount * f  * sigma;
    model = changeRxnBounds(model,protPoolRxn,fluxBound,'u');
end

end

