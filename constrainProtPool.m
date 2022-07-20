function model = constrainProtPool(model,protAmount,f,sigma)
%constrainProtPool Constrains the protein pool reaction according to the
%given amount.
%
% USAGE:
%
%    modelOut = constrainProtPool(model,protAmount,f,sigma)
%
% INPUTS:
%    model:             GECKO/acidFBA model structure (protein pool)
%
%    protAmount:        Total protein dry weight fraction (g/gDW)
%
% OPTIONAL INPUTS:
%    f:                 Mass fraction of proteins accounted for in the model
%
%    sigma:             Average enzyme saturation factor
%   
% OUTPUTS:
%    model:             GECKO/acidFBA model structure with constrained 
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

% Change flux bound
protPoolRxn = params.protPool;
fluxBound = protAmount * f  * sigma;
model = changeRxnBounds(model,protPoolRxn,fluxBound,'u');

end
