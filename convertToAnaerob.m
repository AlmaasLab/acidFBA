function anaerobModel = convertToAnaerob(model)
%convertToAnaerob Converts the acidFBAmodel of the Yeast GEM to an anaerob
%model with the ability to grow without oxygen. Based on the approach defined
%by the yeast-GEM community: 
%https://github.com/SysBioChalmers/yeast-GEM/blob/main/code/otherChanges/anaerobicModel.m
%
% USAGE:
%
%    convertToAnaerob(model)
%
% INPUTS:
%    model:             AcidFBA model structure
%
% .. Authors:
%       - Vetle Simensen 07/12/21

anaerobModel = model;

% Constrain oxygen uptake
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_1992_REV')) = 0;

% Change biomass composition (cofactors)
met_index = ismember(anaerobModel.mets,{'s_3714[c]','s_1198[c]','s_1203[c]','s_1207[c]','s_1212[c]','s_0529[c]'});
anaerobModel.S(met_index,strcmp(anaerobModel.rxns,'r_4598')) = 0;

% Allow for fatty acid and sterol exchange
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_1757_REV')) = 1000;    %ergosterol
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_1915_REV')) = 1000;    %lanosterol
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_1994_REV')) = 1000;    %palmitoleate
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_2106_REV')) = 1000;    %zymosterol
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_2134_REV')) = 1000;    %14-demethyllanosterol
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_2137_REV')) = 1000;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_2189_REV')) = 1000;    %oleate

% Update GAM and NGAM
GAM = 30.49;  % data from Nissen et al. 1997
NGAM = 0;  % data from Nissen et al. 1997
met_index = ismember(anaerobModel.mets,{'s_0394[c]','s_0434[c]','s_0794[c]',...
    's_0803[c]','s_1322[c]'});
anaerobModel.S(met_index,strcmp(anaerobModel.rxns,'r_4041')) = [GAM -GAM GAM -GAM GAM];
anaerobModel = changeRxnBounds(anaerobModel,'r_4046',NGAM,'b');

% Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_0713_REVNo1')) = 0; %Mithocondria
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_0714_REVNo1')) = 0; %Cytoplasm

% Block glycerol dehydroginase (only acts in microaerobic conditions)
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_0487No1')) = 0;

% Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
anaerobModel.ub(strcmp(anaerobModel.rxns,'r_0472No1')) = 0;

end
