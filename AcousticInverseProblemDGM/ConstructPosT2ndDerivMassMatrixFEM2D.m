function PosT2ndDerivMassMatrixFEM = ConstructPosT2ndDerivMassMatrixFEM2D(p0,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,MassCoeffs,PrecompCoeffMatricesInv)

% ConstructPosTMassMatrixFEM2D constructs the FEM mass matrix required for utilizing
% the Newton direction with the adjoint state method with positivity transform.
% A piecewise linear basis is assumed for zeta(p^(m)).
%
% Inputs:
%   p0 - Input vector of size N_Nodes with respect to FEM Mesh
%   MeshIElm: Number of elements by 3 array containing the indices of the nodes in each element
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshINodes: 2 by Number of nodes array for the coordinates of the nodes on the inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   MassCoeffs - Coefficients for the mass matrix
%   PrecompCoeffMatricesInv - N_Elm by 1 cell where each cell contains
%                             the required matrix for computing the coefficients for that element's
%                             interpolating function of two variables. These matrices will be required for
%                             computing the Jacobian of this interpolation function
%
% Outputs:
%   pTransformedMassMatrixFEM - Mass matrix M(zeta(p^(m)))
%
% Hwan Goh 3/05/2018, University of Auckland, New Zealand (Almost been 28 for a whole month!)

%=== Shortening labels ===%
N_Nodes = MeshIN_Nodes;
Nodes = MeshINodes;
N_Elm = MeshIN_Elm;
Elements = MeshIElements;

%=== Organising into 3 by number of elements array ===%
p0 = FEM_Construct3ByN_ElmArray(MeshINodes,MeshIElements,p0);

%=== Mass Matrix ===%
PosT2ndDerivMassMatrixFEM = zeros(N_Nodes,N_Nodes);

%=== Filling the Matrix entries ===%
L = [-1,1,0;-1,0,1]; %Matrix with the gradients of the nodal functions for the reference triangle as its columns.
for K=1:N_Elm
    Vertices=Elements(K,:); %Retrieves the indices of the vertices of the element
    VerticesCoords=Nodes(Vertices,:); %Forms a 3x2 matrix where each row represents a vertex and the columns represent the x and y coordinates.
    JacT = L*VerticesCoords;
    detJacT = abs(det(JacT));  
    PosT2ndDerivMassMatrixFEM(Vertices,Vertices) = PosT2ndDerivMassMatrixFEM(Vertices,Vertices) + pwlIntPosT2ndDerivMass(p0(:,K),VerticesCoords,MassCoeffs(:,K),detJacT,cell2mat(PrecompCoeffMatricesInv(K))); %Fills in the entries of ExpTMassMatrixFEM
end



