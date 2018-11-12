function [mapO,mapD,mapN,vmapO,vmapD,vmapN] = BuildBCMaps2D(BCType,Nfp,vmapM)

% Notes:
% - bnodes is required only in BuildBCMaps2D for constructing specific node maps for different boundary conditions (mapD, mapN, mapW etc)

Out = 2; Dirichlet = 6; Neumann = 7; %Associating boundary condition type with an integer

bct    = BCType';
bnodes = ones(Nfp, 1)*bct(:)';
bnodes = bnodes(:);

% find location of boundary nodes in face and volume node lists
mapO = find(bnodes==Out);          vmapO = vmapM(mapO);
mapD = find(bnodes==Dirichlet);    vmapD = vmapM(mapD);
mapN = find(bnodes==Neumann);      vmapN = vmapM(mapN);
