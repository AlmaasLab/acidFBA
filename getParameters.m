function parameters = getParameters()
%getParameters Get parameters that are used by the acidFBAmodel generation
%and analysis pipeline.
%
% USAGE:
%
%    parameters = getParameters()
%
% OUTPUT:
%    parameters:             Structure of parameters. See descriptions below
%
% .. Authors:
%       - Vetle Simensen 27/08/20

% Names of all 20 proteinogenic amino acids
parameters.aminoAcids = {'alanine','cysteine','aspartate','glutamate','phenylalanine','glycine','histidine','isoleucine','lysine',...
    'leucine','methionine','asparagine','proline','glutamine','arginine','serine','threonine','valine','tryptophan','tyrosine'};

% Single letter amino acid codes
parameters.aaCodes = {'A','C','D','E','F','G','H','I','K','L','M','N', ...
            'P','Q','R','S','T','V','W','Y'}; 
        
% Three letter amino acid codes
parameters.aaThreeCodes = {'Ala','Cys','Asp','Glu','Phe','Gly','His','Ile','Lys','Leu','Met','Asn', ...
            'Pro','Gln','Arg','Ser','Thr','Val','Trp','Tyr'}; 
        
% Corresponding model metabolites
parameters.aaMets = {'A[c]','C[c]','D[c]','E[c]','F[c]','G[c]','H[c]','I[c]','K[c]','L[c]','M[c]',...
    'N[c]','P[c]','Q[c]','R[c]','S[c]','T[c]','V[c]','W[c]','Y[c]'};

% Amino acid source reactions
parameters.aaDrains = {'AAR_pool_A','AAR_pool_C','AAR_pool_D','AAR_pool_E','AAR_pool_F','AAR_pool_G',...
    'AAR_pool_H','AAR_pool_I','AAR_pool_K','AAR_pool_L','AAR_pool_M','AAR_pool_N','AAR_pool_P',...
    'AAR_pool_Q','AAR_pool_R','AAR_pool_S','AAR_pool_T','AAR_pool_V','AAR_pool_W','AAR_pool_Y'}; 

% Molar mass (g/mol) of free amino acids assuming pH ~7.3
parameters.aaMW = [89.09 121.16 132.09 146.12 165.19 75.07 155.15 131.17 147.19 ...
    131.17 149.21 132.12 115.13 146.14 175.21 105.09 119.12 117.15 204.22 181.19];

% Molar mass (g/mol) of protein-bound amino acids assuming pH ~7.3
parameters.aaMWprot = [71.08 103.14 114.08 128.11 147.17 57.05 137.14 113.16 129.18 ...
    113.16 131.20 114.10 97.11 128.13 157.19 87.08 101.10 99.13 186.21 163.17];

% Average enzyme saturation factor
parameters.sigma = 0.5;

% Total protein content in the cell [g protein/gDw]
parameters.Ptot = 0.5;      % assumed constant

% Mass fraction of proteins accounted for in the model
parameters.f = 0.4461;

% Protein pool exchange
parameters.protPool = 'prot_pool_exchange';

% Relevant exchange reactions
parameters.glcUpt = 'r_1714_REV';
parameters.nh4Upt = 'r_1654_REV';
parameters.etohSec = 'r_1761';
parameters.acSec = 'r_1634';
parameters.pyrSec = 'r_2033';
parameters.succSec = 'r_2056';
parameters.glySec = 'r_1808';
parameters.co2Sec = 'r_1672';
parameters.oxyUpt = 'r_1992_REV';

end
