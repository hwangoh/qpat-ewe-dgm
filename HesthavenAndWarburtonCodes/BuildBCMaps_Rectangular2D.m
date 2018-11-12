function [mapO,mapD,mapN,vmapO,vmapD,vmapN,mapTop,mapLeft,mapRight,mapBttm,vmapTop,vmapLeft,vmapRight,vmapBttm] = BuildBCMaps_Rectangular2D(BCType,Nfp,vmapM,TopBndInfo,LeftBndInfo,RightBndInfo,BttmBndInfo)

% Notes:
% - bnodes is required only in BuildBCMaps2D for constructing specific node maps for different boundary conditions (mapD, mapN, mapW etc)
% - this code has the extra function (to the original Hesthaven/Warburton codes) that generates the list of boundary nodes for each of the edges.

Out = 2; Dirichlet = 6; Neumann = 7; %Associating boundary condition type with an integer

bct    = BCType';
bnodes = ones(Nfp, 1)*bct(:)';
bnodes = bnodes(:);

% find location of boundary nodes in face and volume node lists
mapO = find(bnodes==Out);          vmapO = vmapM(mapO);
mapD = find(bnodes==Dirichlet);    vmapD = vmapM(mapD);
mapN = find(bnodes==Neumann);      vmapN = vmapM(mapN);

% find indices of top, left, right and bottom edges
top    = TopBndInfo';
TopNodes = ones(Nfp, 1)*top(:)';
TopNodes = TopNodes(:);

left   = LeftBndInfo';
LeftNodes = ones(Nfp, 1)*left(:)';
LeftNodes = LeftNodes(:);

right   = RightBndInfo';
RightNodes = ones(Nfp, 1)*right(:)';
RightNodes = RightNodes(:);

bttm   = BttmBndInfo';
BttmNodes = ones(Nfp, 1)*bttm(:)';
BttmNodes = BttmNodes(:);

% find location of boundary nodes in face and volume node lists
mapTop = find(TopNodes==1);           vmapTop = vmapM(mapTop);
mapLeft = find(LeftNodes==2);         vmapLeft = vmapM(mapLeft);
mapRight = find(RightNodes==3);       vmapRight = vmapM(mapRight);
mapBttm = find(BttmNodes==4);         vmapBttm = vmapM(mapBttm);

