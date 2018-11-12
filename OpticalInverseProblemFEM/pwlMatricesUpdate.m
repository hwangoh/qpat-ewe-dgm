function [S,M] = pwlMatricesUpdate(Mesh,Prmtrs)

% pwcMatricesUpdate updates the stiffness and mass system matrices for
% Gauss-Newton Iterations when computing the inverse problem. This ensures
% the efficiency of our computation as the absorption coefficient only
% applies to the stiffness and mass matrices. This code is the same as
% pwlFEMMatrices, just without the computation of the R matrix and b vector
% A piecewise linear basis is assumed for the diffusion and scattering coefficients.
%
% Inputs:
%   Mesh:
%      Mesh.N_Nodes - Number of nodes
%      Mesh.Nodes - Coordinates of Nodes
%      Mesh.N_Elm - Number of elements
%      Mesh.Elements - Matrix where row number corresponds to the element number and the entries of the row are the vertices of the element
%
%   Prmtrs:
%      Prmtrs.kappa - 1 x N_elm vector where each entry represents the value of the diffusion coefficient for that element
%      Prmtrs.mu - 1 x N_elm vector where each entry represents the value of the absorption coefficient for that element 
%      Prmtrs.kappa_elmts - Elementwise organisation of kappa values
%      Prmtrs.mu_a_elmts - Elementwise organisation of mu_a values
%
% Outputs:
%   S - Stiffness Matrix
%   M - Mass Matrix
%
% Hwan Goh 31/08/2013, University of Auckland, New Zealand

%Shortening labels
N_Nodes = Mesh.N_Nodes;
Nodes = Mesh.Nodes;
N_elm = Mesh.N_Elm;
Elements = Mesh.Elements;
kappa = Prmtrs.kappa_elmts;
mu_a = Prmtrs.mu_a_elmts;

%Setting up Matrices
S = sparse(N_Nodes, N_Nodes); %Stiffness Matrix
M = sparse(N_Nodes, N_Nodes); %Mass Matrix

%Filling the Matrix entries
L = [-1,1,0;-1,0,1]; %Matrix with the gradients of the nodal functions for the reference triangle as its columns.
for i=1:N_elm
    Ver=Elements(i,:); %Retrieves the indices of the vertices of the element
    E=Nodes(Ver,:); %Forms a 3x2 matrix where each row represents a vertex and the columns represent the x and y coordinates.
    JacT = L*E;
    detJacT = abs(det(JacT));
    S(Ver,Ver) = S(Ver,Ver) + pwlIntStiffness(E,kappa(:,i),L,JacT,detJacT); %Fills in the entries of S
    M(Ver,Ver) = M(Ver,Ver) + pwlIntMass(E,mu_a(:,i),detJacT); %Fills in the entries of M
end
