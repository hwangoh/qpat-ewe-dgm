function [DGMMesh,PrecomputedIntrplteObjects] = ConstructDGMMesh(RunOptions,Mesh)

% ConstructDGMMesh constructs the inversion mesh for DGM, given the properties of the
% optical inversion mesh stored in MeshI.
%
% Inputs:
%   RunOptions:
%       - DGMPolyOrder: polynomial order used for approximation
%       - ScalingOptical - Relates to the units the optical mesh is in. Default is metres.
%   Mesh:
%       - N_Nodes: Number of nodes on the optical inversion mesh
%       - Nodes: N_Nodes by 2 matrix storing the x and y coordinates of Mesh's nodes
%       - Mesh.N_Elm: Number of triangular elements on Mesh
%       - Mesh.Elements: N_Elm by 3 matrix storing the indices of the
%         nodes corresponding to the vertices of the triangular elements of Mesh
%
% Outputs:
%   DGMMesh:
%       - N: Polynomial order used for approximation
%       - Norder: K by 1 matrix where each entry represents the polynomial order of each element
%       - Nfaces: Number of faces for each element. Since we are using triangular faces, then Nfaces = 3
%       - Nv: Number of vertices in triangulated grid
%       - VX: x-coordinates of vertices [m] same as first column of Mesh.Nodes but relabeled to match textbook codes
%       - VY: y-coordinates of vertices [m] same as second column of Mesh.Nodes but relabeled to match textbook codes
%       - EToV: N_Elm by 3 matrix storing the indices of the nodes corresponding to the vertices of the triangular elements of Mesh. Same as Mesh.Nodes but relabeled to match textbook codes
%       - K: Number of elements on the DGM mesh
%       - x: x-coordinate of the nodes of the DGM mesh [m]
%       - y: y-coordinate of the nodes of the DGM mesh [m]
%       - Np: Number of nodes per element on the DGM mesh
%       - VertexNodesGlobalIndicesFEMI: indices referring to the vertices of elements. This is for circumventing the gradient inconsistencies arising from the adjoint state method; only the gradient values at these points are considered
%       - VertexNodesGlobalIndicesDGM: indices referring to the vertices of elements; includes indices for repeated degrees of freedom corresponding to the same spatial point
%       - pinfo:
%              - Np: Number of nodes per p-th order element on the DGM mesh
%              - K: Number of p-th order elements on the DGM mesh
%              - x: x-coordinate of the nodes of the DGM inversion mesh [m] corresponding to nodes in p-th order elements
%              - y: y-coordinate of the nodes of the DGM inversion mesh [m] corresponding to nodes in p-th order elements
%              - ids: Np(p)*K(p) indices for nodes corresponding to pth-order elements 
%              - Fmask, r, s , Dr, Ds, LIFT, V, rx, sx, ry, sy, nx, ny, sJ, J, nxny, nxnx, nyny, fmapM, interpP: objects required for computation of PDE discretized by the DGM
%       If rectangular mesh:
%       - mapTop:    indices referring to the entries of vmapP associated with the coordinates of the top edge nodes
%       - mapLeft:   indices referring to the entries of vmapP associated with the coordinates of the left edge nodes
%       - mapRight:  indices referring to the entries of vmapP associated with the coordinates of the right edge nodes
%       - mapBttm:   indices referring to the entries of vmapP associated with the coordinates of the bottom edge nodes
%       - vmapTop:   indices referring to the coordinates of the top edge nodes
%       - vmapLeft:  indices referring to the coordinates of the left edge nodes
%       - vmapRight: indices referring to the coordinates of the right edge nodes
%       - vmapBttm: indices referring to the coordinates of the bottom edge nodes
%   PrecompIntrplteObjects: (These are used for when interpolation on the same meshes are reused as well as computing the Jacobian of the intrplteFEM2DGM function)
%       - DT - Delaunay triangulated mesh required to use the pointLocation function. Therefore corresponds with the
%              PointLocs vector.
%       - PrecompCoeffMatricesInv - N_Elm by 1 cell where each cell contains
%                                the required matrix for computing the coefficients for that element's
%                                interpolating function of two variables. These matrices will be required for
%                                computing the Jacobian of this interpolation function
%       - PointLocations - Np by 1 vector containing element numbers (with reference to DT) each grid point lies in. 
%                                This will be required for computing the Jacobian of this
%                                interpolation function.
%
% Hwan Goh, University of Auckland, New Zealand 01/01/2016 - Happy New Year!!!
% Last Edited: 16/11/2017 - removed reliance on Globals2D and save everything into structures

%=== Renaming variables to match the textbook provided codes ===%
N = RunOptions.DGMPolyOrder; %Polynomial order used for approximation
Nfaces = 3; %Since we are using triangular elements
Nv = Mesh.N_Nodes; %Number of vertices in triangulated grid
VX = Mesh.Nodes(:,1)'*(1/RunOptions.ScalingOptical); %x-coordinates of vertices [m]
VY = Mesh.Nodes(:,2)'*(1/RunOptions.ScalingOptical); %y-coordinates of vertices [m]
K = Mesh.N_Elm; %Number of triangular elements
EToV = Mesh.Elements; %Element to vertex matrix

%=== Timo's Method for ensuring anti-clockwise indexing of element boundary nodes ===%
%Note: this is also added at the end of TriangulateRect when the optical mesh is generated for QPAT
for ii = 1:K
  iv = EToV(ii,:); 
  x1 = VX(iv);
  y1 = VY(iv);
  rpx = mean(x1);
  rpy = mean(y1);
  for jj = 1:3
    angles(jj) = atan2((y1(jj)-rpy),(x1(jj)-rpx))/pi*180;
  end
  [~, id] = sort(angles);
  EToV(ii,:) = EToV(ii,id);
end

%=== Initialize solver and construct grid and metric ===%
Norder = RunOptions.DGMPolyOrder*ones(K,1);
if RunOptions.UseTrelisMesh == 1 %Trelis generated circular mesh, BCType is already provided by Timo's ReadGrid_trelis.m
    Out = 2; Dirichlet = 6; Neumann = 7; %Associating boundary condition type with an integer
    Mesh.BCType = eval(RunOptions.BoundaryCondition)*Mesh.BCType;
    [pinfo,x,y,Np] = BuildPNonCon2D(K,VX,VY,EToV,Norder,Mesh.BCType);
end
if RunOptions.UseMyRectangularMesh == 1 %My generated rectangular mesh
    [BCType,TopBndInfo,LeftBndInfo,RightBndInfo,BttmBndInfo] = ConstructDGMBndObjects_Rectangular2D(RunOptions,Nfaces,K,VX,VY,EToV);
    [pinfo,x,y,Np,mapTop,mapLeft,mapRight,mapBttm,vmapTop,vmapLeft,vmapRight,vmapBttm] = BuildPNonCon_Rectangular2D(K,VX,VY,EToV,Norder,BCType,TopBndInfo,LeftBndInfo,RightBndInfo,BttmBndInfo);
end

%=== Global Node Indices ===%
VertexNodesxy = [VX(:),VY(:)];
Nodesxy = [x(:),y(:)];
% VertexNodesGlobalIndicesDGM = [];

% No longer needed since I have the transpose of the interpolation matrix
% for i = 1:size(VX(:),1)
%     Coords = VertexNodesxy(i,:);
%     rowstemp = find(ismember(Nodesxy,Coords,'rows'));
%     VertexNodesGlobalIndicesFEM(i) = rowstemp(1); %Element vertices for FEM mesh
%     VertexNodesGlobalIndicesDGM = [VertexNodesGlobalIndicesDGM;rowstemp]; %Element vertices for DGM mesh
% end
% VertexNodesGlobalIndicesFEM = VertexNodesGlobalIndicesFEM(:);

%=== Construct Interpolation Objects ===%
[PrecomputedIntrplteObjects] = PrecomputeIntrplteObjects(Mesh.N_Elm,Mesh.Nodes*(1/RunOptions.ScalingOptical),x,y,Np*K);

%=== Saving Outputs in a Structure ===%
DGMMesh.N = N;
DGMMesh.Norder = Norder;
DGMMesh.Nfaces = Nfaces;
DGMMesh.Nv = Nv;
DGMMesh.VX = VX;
DGMMesh.VY = VY;
DGMMesh.EToV = EToV;
DGMMesh.K = K;
DGMMesh.x = x;
DGMMesh.y = y;
DGMMesh.Np = Np;
% DGMMesh.VertexNodesGlobalIndicesFEM = VertexNodesGlobalIndicesFEM;
% DGMMesh.VertexNodesGlobalIndicesDGM = VertexNodesGlobalIndicesDGM;
DGMMesh.pinfo = pinfo;

if RunOptions.UseMyRectangularMesh == 1
    DGMMesh.mapTop = mapTop;
    DGMMesh.mapLeft = mapLeft;
    DGMMesh.mapRight = mapRight;
    DGMMesh.mapBttm = mapBttm;
    DGMMesh.vmapTop = vmapTop;
    DGMMesh.vmapLeft = vmapLeft;
    DGMMesh.vmapRight = vmapRight;
    DGMMesh.vmapBttm = vmapBttm;
end








