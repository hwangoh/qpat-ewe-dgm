function [M] = ConstructMassMatrixFEM2D(Elements,Nodes)

% ConstructMassMatrixFEM2D constructs the basic mass matrix for FEM where there are no spatially
% dependent coefficients in the integrals 
%
% Inputs:
%   Mesh:
%      Elements - Matrix where row number corresponds to the element number and the entries of the row are the vertices of the element
%      Nodes - Coordinates of Nodes
%
% Outputs:
%      M - Mass Matrix
%
% Hwan Goh 19/09/2017, University of Auckland, New Zealand
% Adapted from Pascal Cheon, I use "adapted" in loosest sense of the word
    
N_Elements = size(Elements,1); % number of elements
N_Nodes = size(Nodes,1); % number of nodes
L = [-1,1,0;-1,0,1]; %Matrix with the gradients of the nodal functions for the reference triangle as its columns

M = zeros(N_Nodes,N_Nodes);
for K = 1:N_Elements
    Vertices = Elements(K,:); %Retrieves the indices of the vertices of elements
    VerticesCoords = Nodes(Vertices,:); %Forms a 4x2 matrix where each row represents a vertex and the columns represent the x and y coordinates
    JacT = L*VerticesCoords;
    detJacT = abs(det(JacT));
    A = (1/24)*[2 1 1; 1 2 1; 1 1 2];
    M(Vertices,Vertices) = M(Vertices,Vertices) + detJacT*A;
end