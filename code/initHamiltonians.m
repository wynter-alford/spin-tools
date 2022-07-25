%% initHamiltonians.m
% Wynter Alford
% 
% Adapted from pre-existing code into its own script for ease of use

if N==400&&couplingsCount<=4
    load('Uniform02_Dipole_Hamiltonians(4).mat','Hdips');
elseif ~exist('Hdips','var')
     Hdips = cell(couplingsCount,1);
     % Generate [couplingsCount] dipole Hamiltonians, each with a different coupling matrix
     for j=1:couplingsCount
         dip = random('Uniform',0.1,2,N);
         dip = triu(dip,1) + triu(dip,1)';
         Hdips{j} = getHdip(N, dim, x, y, z, dip);
     end
end