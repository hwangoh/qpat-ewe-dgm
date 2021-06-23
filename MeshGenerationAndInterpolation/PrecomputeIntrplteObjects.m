function [PrecomputedIntrplteObjects] = PrecomputeIntrplteObjects(N_Elm,Nodes,x,y,N_points)

% PrecompIntrplteObjects computes all the objects required for
% interpolation that depend only on the mesh.
%
% NOTE: Make sure units of nodes and units of x and y match as inputs!
%
% Inputs:
%   N_Elm - Number of triangulated elements on the mesh to be interpolated over
%   Nodes - N_Nodes by 2 matrix where the ith row contains the x and y coordinates of the ith node on the mesh to be interpolated over
%   x: x coordinates of the points to be interpolated
%   y: x coordinates of the points to be interpolated
%   N_points: number of points to be interpolated
%
% Outputs:
%      PrecompIntrplteObjects: (These are used for when interpolation on the same meshes are reused as well as computing the Jacobian
%                               of the intrplteFEM2DGM function)
%                            DT - Delaunay triangulated mesh required to use the pointLocation function. Therefore corresponds with the
%                                 PointLocs vector.
%                            PrecompCoeffMatricesInv - N_Elm by 1 cell where each cell contains
%                                the required matrix for computing the coefficients for that element's
%                                interpolating function of two variables. These matrices will be required for
%                                computing the Jacobian of this interpolation function
%                            PointLocations - Np by 1 vector containing element numbers (with reference to DT) each grid point lies in. 
%                            InterpMatrix - N_points by N matrix that represents the interpolation operator
%
% Hwan Goh 22/07/2015, University of Auckland, New Zealand

%=== Triangulation ===%
PrecomputedIntrplteObjects.DT=DelaunayTri(Nodes(:,1),Nodes(:,2))'; %if you want to use the pointLocation function, you have to store the topology matrix as a Delaunay type structure; Mesh.Elements won't work

%=== Inverted Matrices ===%
PrecomputedIntrplteObjects.PrecompCoeffMatricesInv = cell(N_Elm,1);
for i=1:N_Elm;
    Node_1 = PrecomputedIntrplteObjects.DT(i,1);
    Node_2 = PrecomputedIntrplteObjects.DT(i,2);
    Node_3 = PrecomputedIntrplteObjects.DT(i,3);
    coordnts = [Nodes(Node_1,:);Nodes(Node_2,:);Nodes(Node_3,:)];
    A = [1, coordnts(1,1), coordnts(1,2); 1, coordnts(2,1), coordnts(2,2);
        1, coordnts(3,1), coordnts(3,2)];
    PrecomputedIntrplteObjects.PrecompCoeffMatricesInv{i} = inv(A);
end

%=== Grid Point Locations ===%
for i = 1:N_points;
    PrecomputedIntrplteObjects.PointLocations(i) = pointLocation(PrecomputedIntrplteObjects.DT,[x(i),y(i)]); %returns the element that the gridpoint lies in.
end

%=== Interpolation Matrix ===%
PrecomputedIntrplteObjects.InterpMatrix = sparse(N_points,size(Nodes,1));
for ii=1:N_points
    ElementIndex = PrecomputedIntrplteObjects.PointLocations(ii);
    Node_1 = PrecomputedIntrplteObjects.DT(ElementIndex,1);
    Node_2 = PrecomputedIntrplteObjects.DT(ElementIndex,2);
    Node_3 = PrecomputedIntrplteObjects.DT(ElementIndex,3);
    B = cell2mat(PrecomputedIntrplteObjects.PrecompCoeffMatricesInv(ElementIndex));
    PrecomputedIntrplteObjects.InterpMatrix(ii,Node_1) = B(1,1) + B(2,1)*x(ii) + B(3,1)*y(ii);
    PrecomputedIntrplteObjects.InterpMatrix(ii,Node_2) = B(1,2) + B(2,2)*x(ii) + B(3,2)*y(ii);
    PrecomputedIntrplteObjects.InterpMatrix(ii,Node_3) = B(1,3) + B(2,3)*x(ii) + B(3,3)*y(ii);
end

%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Backup: Precomputed Coefficients using Normal Vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating Coefficients, Normal Vectors Approach %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if RunOptions.IntrplteNormalVecs == 1;
% for i=1:N_Elm;
%     %=== Obtaining Values ===%
%     Node_1 = DT(i,1);
%     Node_2 = DT(i,2);
%     Node_3 = DT(i,3);
%     coordnts = [Nodes(Node_1,:);Nodes(Node_2,:);Nodes(Node_3,:)];
%     value = [values(Node_1), values(Node_2), values(Node_3)];
%     %=== Setting the Initial Point ===%
%     x0 = coordnts(1,1);
%     y0 = coordnts(1,2);
%     z0 = value(1);
%     iniPt = [x0,y0,z0];
%     %=== Constructing the Normal Vector ===%
%     v1 = [coordnts(2,1),coordnts(2,2),value(2)] - iniPt;
%     v2 = [coordnts(3,1),coordnts(3,2),value(3)] - iniPt;
%     PrecompIntrplteObjects.PrecompCoeff(:,i) = cross(v1,v2);
%     %=== Coefficients Approach with Normal Vector ===%
%     %PrecompCoeff2{i} = [(n(1)*x0 +n(2)*y0 - n(3)*z0)/n(3), -(n(1)/n(3)), -(n(2)/n(3))];    
% end
% end
