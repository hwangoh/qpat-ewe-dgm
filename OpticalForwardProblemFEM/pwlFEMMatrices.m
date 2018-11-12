function [S,M,R,b] = pwlFEMMatrices(Mesh,Prmtrs)

% pwcFEMMatrices constructs the blocks needed to form the system matrix for
% FEM where a piecewise linear basis is assumed for the diffusion and scattering coefficients.
%
% Inputs:
%   Mesh:
%      N_Nodes - Number of nodes
%      Nodes - Coordinates of Nodes
%      N_Elm - Number of elements
%      Elements - Matrix where row number corresponds to the element number and the entries of the row are the vertices of the element
%      Lght_ElmtInd - A (lght_Nelm x L) matrix where each column represents the set of elements under one light source
%      Bnd_ElmInd - A column vector listing the element indices for elements with an edge on the boundary
%      I_s - Diffuse boundary condition  
%   Prmtrs:
%      Prmtrs.kappa_elmts - 3 by N_Elm matrix elementwise organisation of kappa values
%      Prmtrs.mu_a_elmts - 3 by N_Elm matrix elementwise organisation of mu_a values
%
% Outputs:
%   S - Stiffness Matrix
%   M - Mass Matrix
%   R - Boundary Matrix
%   b - Lightsource Vector
%
% Hwan Goh 04/06/2013, University of Auckland, New Zealand

%% Shortening labels
N_Nodes = Mesh.N_Nodes;
Nodes = Mesh.Nodes;
N_elm = Mesh.N_Elm;
Elements = Mesh.Elements;
Lght_ElmtInd = Mesh.Lght_ElmtInd;
Bnd_ElmInd = Mesh.Bnd_ElmInd;
I_s = Mesh.I_s;
kappa = Prmtrs.kappa_elmts;
mu_a = Prmtrs.mu_a_elmts;

%% Setting up Matrices
S = sparse(N_Nodes, N_Nodes); %Stiffness Matrix
M = sparse(N_Nodes, N_Nodes); %Mass Matrix
R = sparse(N_Nodes, N_Nodes); %Boundary Matrix
b = sparse(N_Nodes,1)'; %Lightsource Vector

%% Filling the Matrix entries
L = [-1,1,0;-1,0,1]; %Matrix with the gradients of the nodal functions for the reference triangle as its columns.
for K=1:N_elm
    Vertices=Elements(K,:); %Retrieves the indices of the vertices of the element
    VerticesCoords=Nodes(Vertices,:); %Forms a 3x2 matrix where each row represents a vertex and the columns represent the x and y coordinates.
    JacT = L*VerticesCoords;
    detJacT = abs(det(JacT));
    
    S(Vertices,Vertices) = S(Vertices,Vertices) + pwlIntStiffness(VerticesCoords,kappa(:,K),L,JacT,detJacT); %Fills in the entries of S
    M(Vertices,Vertices) = M(Vertices,Vertices) + pwlIntMass(VerticesCoords,mu_a(:,K),detJacT); %Fills in the entries of M
    
    %Check if element i has an edge on the boundary
    if any(find(K==Bnd_ElmInd)) 
        [R(Vertices,Vertices),b(Vertices)] = bndIntgrl(Mesh,VerticesCoords,K,I_s,Lght_ElmtInd); %Fills in the entry of R. Notice that E lists the coordinates of node that's not on the boundary in the 3rd row.
    end
end