function sensors = EWE_DGM2D_ConstructSensorArray(RunOptions,SensorCoords,VX,VY,EToV,Bnd_ElmInd,Bnd_NodeInd,pinfo,Norder)

% EWE_DGM2D_ConstructSensorArray forms the sensor array by finding elements that 
% contains and receivers.
%
% Inputs: 
%    SensorCoords: Number of sensors by 2 matrix where each row represents the coordinates of the sensor
%    VX: x-coordinates of element vertices [m]
%    VY: y-coordinates of element vertices [m]
%    EToV: N_Elm by 3 matrix storing the indices of the nodes corresponding to the vertices of the triangular elements of Mesh.
%    Bnd_ElmInd: Indices of boundary elements, to ensure that the correct element is found
%    Bnd_NodeInd: Indices of boundary nodes, to ensure that the correct element is found
%    pinfo:
%        - Np: Number of nodes per p-th order element on the DGM mesh
%        - K: Number of p-th order elements on the DGM mesh
%        - x: x-coordinate of the nodes of the DGM inversion mesh [m] corresponding to nodes in p-th order elements
%        - y: y-coordinate of the nodes of the DGM inversion mesh [m] corresponding to nodes in p-th order elements
%        - ids: Np(p)*K(p) indices for nodes corresponding to pth-order elements 
%        - Fmask, r, s , Dr, Ds, LIFT, V, rx, sx, ry, sy, nx, ny, sJ, J, nxny, nxnx, nyny, fmapM, interpP: objects required for computation of PDE discretized by the DGM
%    Norder: Max order of polynomial interpolation
%
% Output: 
%     sensors: 1 by number of sensors array containing 4 cells 
%        - id: Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%        - xy: coordinates of the sensor
%        - l_iatsensor: 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%
% Created by Timo Lahivaara, his function is called LocateSourceAndReceivers
% Modified: Timo L. Mar. 16, 2016
%
% Adapted by Hwan Goh, University of Auckland, New Zealand 19/11/2017
% Last Modified: 19/11/2017
%
% Timo's Notes:
% does tsearchn work here?
% there should a check if tsearchn returns nan
%
% this code can be optimized so that it first finds elements that are found
% for each basis order - not so major effect if the number of receivers
% remains "small"
% this will can change the ordering on receiver coordinates and so it must
% be taken into account if modifying the code
%
% toimisko joku find(ismember(pinf.Kids, all)) / global element style might
% be faster - should be modified

%% =======================================================================%
%                      Finding Elements for Receivers
%=========================================================================%
for ii = 1:length(SensorCoords(:,1))
    elm_rec(ii) = tsearchn([VX(:) VY(:)],EToV,SensorCoords(ii,:)); %list where the ith entry represents the element the ith sensor lies in
    if isnan(elm_rec(ii)) %if the sensor does not lie in an element
        warning(' something wrong with the receiver setup since element that holds the receiver cannot be found --> using the closest available node');
        [~, id] = min((VX(:)-SensorCoords(ii,1)).^2 + (VY(:)-SensorCoords(ii,2)).^2);
        fprintf(' old receiver: [%0.5f %0.5f], new receiver: [%0.5f %0.5f]\n', SensorCoords(ii,1), SensorCoords(ii,2), VX(id), VY(id));
        SensorCoords(ii,:) = [VX(id) VY(id)];
        elm_rec(ii) = tsearchn([VX(:) VY(:)],EToV,SensorCoords(ii,:));
    end
    if RunOptions.UseFullDomainData ~=1 && ~any(elm_rec(ii)==Bnd_ElmInd) %If the element found by tsearchn is not a boundary element
        for j=1:3
            if any(EToV(elm_rec(ii),j)==Bnd_NodeInd)
                NodeInd = EToV(elm_rec(ii),j);
            end
        end
        [BndElm,~] = find(NodeInd==EToV(Bnd_ElmInd,:));
        elm_rec(ii) = Bnd_ElmInd(BndElm(1));
    end
end

% %=== Plotting Elements Where Sensors Are Located ===%
% figure
% triplot(EToV,VX,VY)
% hold on
% %=== Boundary Elements ===%
% for ii = 1:size(elm_rec,2) %Change this to highlight specific elements
%     coord1 = VX(EToV(elm_rec(ii),:));
%     coord2 = VY(EToV(elm_rec(ii),:));
%     patch(coord1(:),coord2(:),[1 0 0])
% end
% axis 'image'

%% =======================================================================%
%                            Constructing Sensors
%=========================================================================%
if ~isempty(elm_rec) %this ensures that elements that are already found are not used in the loop for higher order basis
    elmstofind_rec = 1:length(elm_rec);
else
    elmstofind_rec = [];
end
x = []; y = [];
for N1=1:max(Norder)
    pinf = pinfo(N1);
    if(pinf.K>0)        
        % collect x and y values from each order
        x(pinf.ids) = pinf.x;
        y(pinf.ids) = pinf.y;
        % interpolation for receivers
        found_rec_id = [];
        % for finding normal vectors. We reshape vmapP from Nfp by 3*K to 3*Nfp by K. This way each column represents an element.
        for ii = elmstofind_rec;%1:length(elm_rec)
            elm = find(ismember(pinf.ks, elm_rec(ii))); %pinf.ks is the list of element numbers that have N1-th order polynomial basis
            if elm
                sensors{ii}.id = pinf.ids(:,elm);
                sensors{ii}.xy = [SensorCoords(ii,1) SensorCoords(ii,2)];
                %=== compute interpolation ===%
                [r0, s0] = EWE_DGM2D_Getrs([VX(EToV(elm_rec(ii),:)) VY(EToV(elm_rec(ii),:))], SensorCoords(ii,1), SensorCoords(ii,2));
                Psi = Vandermonde2D(Norder(elm_rec(ii)), r0, s0); %1 by Np Vandermonde matrix using the single [r0,s0] point on the reference element. This is row vector of modal basis functions [Psi_1(r_0,s_0),Psi_2(r_0,s_0),...,Psi_Np(r_0,s_0)]
                %=== compute nodal set ===%
                [xs,ys] = Nodes2D(Norder(elm_rec(ii))); %warping nodes
                [r,s] = xytors(xs,ys);
                %=== build reference element matrices ===%
                V = Vandermonde2D(Norder(elm_rec(ii)), r, s);
                %=== build interpolation matrix ===%
                sensors{ii}.l_iatsensor = Psi*inv(V); %This the local basis function [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)]
                found_rec_id(end+1, 1) = ii;
            end
        end
        elmstofind_rec = setdiff(elmstofind_rec, found_rec_id);
    end
end