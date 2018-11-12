function PosTMassMatrixFEM = ConstructPosTMassMatrixFEM2D(MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,MassCoeffs)

% ConstructPosTMassMatrixFEM2D constructs the FEM mass matrix required for utilizing
% line search with the adjoint state method with a positivity transform.
% A piecewise linear basis is assumed for p(varrho^(m)).
%
% Inputs:
%   MeshIElm: Number of elements by 3 array containing the indices of the nodes in each element
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshINodes: 2 by Number of nodes array for the coordinates of the nodes on the inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   MassCoeffs - POsitivity transformed mass matrix coefficients
%
% Outputs:
%   PosTMassMatrixFEM - Mass matrix M(zeta(p^(m)))
%
% Hwan Goh 3/05/2018, University of Auckland, New Zealand (Almost been 28 for a whole month!)

% Shortening labels
N_Nodes = MeshIN_Nodes;
Nodes = MeshINodes;
N_elm = MeshIN_Elm;
Elements = MeshIElements;

% Mass Matrix
PosTMassMatrixFEM = zeros(N_Nodes,N_Nodes);

% Filling the Matrix entries
L = [-1,1,0;-1,0,1]; %Matrix with the gradients of the nodal functions for the reference triangle as its columns.
for K=1:N_elm
    Vertices=Elements(K,:); %Retrieves the indices of the vertices of the element
    VerticesCoords=Nodes(Vertices,:); %Forms a 3x2 matrix where each row represents a vertex and the columns represent the x and y coordinates.
    JacT = L*VerticesCoords;
    detJacT = abs(det(JacT));  
    PosTMassMatrixFEM(Vertices,Vertices) = PosTMassMatrixFEM(Vertices,Vertices) + pwlIntMass(VerticesCoords,MassCoeffs(:,K),detJacT); %Fills in the entries of ExpTMassMatrixFEM
end