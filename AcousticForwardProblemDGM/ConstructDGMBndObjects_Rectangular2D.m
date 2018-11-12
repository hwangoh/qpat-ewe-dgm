function [BCType,TopBndInfo,LeftBndInfo,RightBndInfo,BttmBndInfo] = ConstructDGMBndObjects_Rectangular2D(RunOptions,Nfaces,K,VX,VY,EToV)

% ConstructDGMBndObjects_Rectangular2D constructs a list of the boundary nodes and
% boundary vertices as well as the KxNfaces matrix BCType which are
% required for constructing the DGM mesh using Hesthaven/Warburton's method.
% Specifically, this is for my rectangular mesh. If using Trelis, BCType is
% already generated using Timo's ReadGrid_trelis.m code and this code is not required.
%
% Inputs:
%   RunOptions:
%       - TopBC: Boundary conditions for top edge
%       - LeftBC: Boundary conditions for left edge
%       - RightBC: Boundary conditions for right edge
%       - BottomBC: Boundary conditions for bottom edge
%   Mesh:
%       - Nfaces: Number of faces in each element
%       - K: Number of elements
%       - VX: x-coordinates of vertices [m] same as first column of Mesh.Nodes but relabeled to match textbook codes
%       - VY: y-coordinates of vertices [m] same as first column of Mesh.Nodes but relabeled to match textbook codes
%       - EToV: N_Elm by 3 matrix storing the indices of the nodes corresponding to the vertices of the triangular elements of Mesh. Same as Mesh.Nodes but relabeled to match textbook codes
%
% Outputs:
%   BCType -       N_Elm by 3 matrix where each entry represents what boundary condition the edge is associated with. If not a boundary node, the entry is 0. This is used later for constructing bnodes in BuildBCMaps2D as well as for p-nonconforming meshes as in the script BuildPNonCon2D. Notice that only two elements have two non-zero entries because only two elements have two edges on the boundary
%   TopBndInfo -   N_Elm by 3 matrix where an entry '1' represents whether it is a top edge boundary edge. It not a boundary edge, the entry is 0. This is used later for constructing bnodes in BuildBCMaps2D as well as for p-nonconforming meshes as in the script BuildPNonCon2D.
%   LeftBndInfo -  N_Elm by 3 matrix where an entry '2' represents whether it is a left edge boundary edge. It not a boundary edge, the entry is 0. This is used later for constructing bnodes in BuildBCMaps2D as well as for p-nonconforming meshes as in the script BuildPNonCon2D.
%   RightBndInfo - N_Elm by 3 matrix where an entry '3' represents whether it is a right edge boundary edge. It not a boundary edge, the entry is 0. This is used later for constructing bnodes in BuildBCMaps2D as well as for p-nonconforming meshes as in the script BuildPNonCon2D.
%   BttmBndInfo -  N_Elm by 3 matrix where an entry '4' represents whether it is a bttm edge boundary edge. It not a boundary edge, the entry is 0. This is used later for constructing bnodes in BuildBCMaps2D as well as for p-nonconforming meshes as in the script BuildPNonCon2D.
%
% Hwan Goh, University of Auckland, New Zealand 1/02/2017
% Last Edited: 16/11/2017 - removed reliance on Globals2D and save everything into structures
%
% Notes: Due to rotating elements, constructing list of boundary nodes 
%        using this script doesn't work. CorrectBCTable checks to see if two nodes 
%        on the boundary share the same element edge. However, there are
%        some nodes on the boundary that are the vertex of another element
%        and so these are missed.

Out = 2; Dirichlet = 6; Neumann = 7; %Associating boundary condition type with an integer
epsilon = 1e-14; %tolerance

%=== Boundary Vertices ===%
[DGMMesh_bttmBndVerts,~] = find(VY(:) <= min(VY(:)) + epsilon);
[DGMMesh_leftBndVerts,~] = find(VX(:) <= min(VX(:)) + epsilon);
[DGMMesh_rightBndVerts,~] = find(VX(:) >= max(VX(:)) - epsilon);
[DGMMesh_topBndVerts,~] = find(VY(:) >= max(VY(:)) - epsilon);
DGMMeshBndVerts = [DGMMesh_bttmBndVerts,DGMMesh_leftBndVerts,DGMMesh_rightBndVerts,DGMMesh_topBndVerts]; %removed sort(...), wasn't required?
DGMMeshBndVerts = DGMMeshBndVerts(:);

%=== Construct BCType ===%
BCType = zeros(K,3);
[BCType,TopBndInfo,LeftBndInfo,RightBndInfo,BttmBndInfo] = CorrectBCTable(Nfaces,K,VX,VY,EToV,BCType,DGMMesh_topBndVerts,DGMMesh_leftBndVerts,DGMMesh_rightBndVerts,DGMMesh_bttmBndVerts,...
                                                                          eval(RunOptions.TopBC),eval(RunOptions.LeftBC),eval(RunOptions.RightBC),eval(RunOptions.BottomBC));
