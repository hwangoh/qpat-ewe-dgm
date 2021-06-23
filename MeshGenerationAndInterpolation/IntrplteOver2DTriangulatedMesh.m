function [IntpltedValues,PrecomputedIntrplteObjects] = IntrplteOver2DTriangulatedMesh(N_Elm,Nodes,values,x,y,N_points,PrecomputedIntrplteObjects)

% IntrplteOver2DTriangulatedMesh interpolates points onto a given 2D
% triangulated mesh with set values at the vertices of each triangular
% element.
% Note that DT denotes the elements differently from QPATkWaveMesh2D.m,
% however the FEM node numbers are the same. Therefore, the element number
% found by the inbuilt Matlab function pointLocation will be different, but
% the corresponding node number are the same.
%
% NOTE: Make sure units of nodes and units of x and y match as inputs!
%
% Inputs:
%   N_Elm - Number of elements in the triangulated mesh on the mesh to be interpolated over
%   Nodes - N_Nodes by 2 matrix where the ith row contains the x and y coordinates of the ith node on the mesh to be interpolated over
%   values - N_Node by 1 vector containing the values at each node of the Optical forward problem mesh
%   x - x coordinates of nodes to be interpolated
%   y - y coordinates of nodes to be interpolated
%   N_points - Number of points to be interpolated
%   PrecomputedIntrplteObjects - Objects that depend on the Mesh nodes. May have been computed earlier and so can 
%                            be called here to avoid repeating computations. Set to 0 in function call if you
%                            want to compute new objects for a new mesh.
%
% Outputs:
%      intpltedValues - N_Points by 1 vector containing the interpolated
%                       values
%      PrecompIntrplteObjects: (These are used for when interpolation on the same meshes are reused as well as computing the Jacobian
%                               of the intrplteFEM2DGM function)
%                 DT - Delaunay triangulated mesh required to use the pointLocation function. Therefore corresponds with the
%                      PointLocs vector.
%                 PrecompCoeffMatricesInv - N_Elm by 1 cell where each cell contains
%                                the required matrix for computing the coefficients for that element's
%                                interpolating function of two variables. These matrices will be required for
%                                computing the Jacobian of this interpolation function
%                 PointLocations - N_points by 1 vector containing element numbers (with reference to DT) each grid point lies in. 
%                                This will be required for computing the Jacobian of this
%                                interpolation function.
%                 PrecompIntrplteCoeff - 3 by N_Elm matrix where each column contains the three coefficients required for
%                                interpolation over an element
%                 InterpMatrix - N_points by N matrix that represents the interpolation operator
%
% Hwan Goh 19/12/2015, University of Auckland, New Zealand

%Whether to compute only the coefficients, whilst reusing precomputed structures
if isstruct(PrecomputedIntrplteObjects)~=1;
    ComputeNewIntrplteObjects = 1;
else
    ComputeNewIntrplteObjects = 0;
end

%=== Setting Up ===%
IntpltedValues = zeros(N_points,1);
%test = zeros(N_pts,2); %For debugging

%=== Precomputing Mesh Objects for linfunct2D ===%
if ComputeNewIntrplteObjects == 1;
    [PrecomputedIntrplteObjects] = PrecomputeIntrplteObjects(N_Elm,Nodes,x,y,N_points);
end
%=== Shortening Labels Again ===%
PointLocations = PrecomputedIntrplteObjects.PointLocations;
DT = PrecomputedIntrplteObjects.DT;
PrecompCoeffMatricesInv = PrecomputedIntrplteObjects.PrecompCoeffMatricesInv;

%=== Computing New Coefficients for New Nodal Values ===%
PrecomputedIntrplteObjects.PrecompIntrplteCoeff = zeros(3,N_Elm);
for i = 1:N_Elm
    Node_1 = DT(i,1);
    Node_2 = DT(i,2);
    Node_3 = DT(i,3);
    z = [values(Node_1); values(Node_2); values(Node_3)];
    B = cell2mat(PrecompCoeffMatricesInv(i));
    PrecomputedIntrplteObjects.PrecompIntrplteCoeff(:,i) = B*z;
end

%=== Computing New Interpolated Values Using Linear Function of Two Variables ===%
for i=1:N_points
    Coeff = PrecomputedIntrplteObjects.PrecompIntrplteCoeff(:,PointLocations(i));
    IntpltedValues(i) = Coeff(1) + Coeff(2)*x(i) + Coeff(3)*y(i);
end
