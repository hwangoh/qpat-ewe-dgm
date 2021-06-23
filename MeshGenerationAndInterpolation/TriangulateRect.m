function [Mesh]=TriangulateRect(RunOptions,Mesh)

% TriangulateRect triangulates the rectangular pixelated grid 
%
% Inputs:
%   RunOptions:
%      ScalingOptical - Relates to the units the optical mesh is in.
%                       Default is metres.
%      RotateCornerElements - Rotate the corner elements so that every corner
%                             is shared by two elements. Note that this
%                             will create elements with different
%                             orientations.
%   Mesh:
%      Dimensns - Dimensions of the mesh [length width] [mm]
%      ExtndBndx - Extension of boundary in the x direction [mm]
%      ExtndBndy - Extension of boundary in the y direction [mm]
%
% Outputs:
%   Mesh:
%        Mesh.Dimesns - Dimensions of mesh [length, width] [mm]
%        Mesh.Nodes - Coordinates of nodes
%        Mesh.Elements - Matrix where row number corresponds to the element number and the entries of the row are the vertices of the element
%        Mesh.N_Nodes - Number of nodes
%        Mesh.N_Elm - Number of elements
%        Mesh.N_Lght - Number of light sources
%        Mesh.Bnd_ElmInd - A column vector listing the element indices for elements with an edge on the boundary
%        Mesh.Bnd_NodeInd - A column vector listing the node indices for the nodes on the boundary
%        Mesh.Lght_ElmtInd - A (lght_Nelm x L) matrix where each column represents the set of elements under one light source
%        Mesh.Lght_Nelm - Maximum number of elements for one light source
%        Mesh.I_s - Diffuse boundary condition
%        Mesh.minDist - Minimum distance between nodes
%        Mesh.maxDist - Maximum distance between nodes
%        Mesh.Maxsize - Maximum length of an edge of an element [mm]
%        Mesh.Minsize - Minimum length of an edge of an element [mm]
%        Mesh.ExtndBndx - Extension of boundary in the x direction [mm]
%        Mesh.ExtndBndy - Extension of boundary in the y direction [mm]
%        Mesh.CornerNodes - [BttmLftNode;BttmRghtNode;TopLftNode;TopRghtNode];
%        Mesh.CornerNodesCoord - [Nodes(BttmLftNode,:);Nodes(BttmRghtNode,:);Nodes(TopLftNode,:);Nodes(TopRghtNode,:)];
%        Mesh.DomainComplmnt - Indices of the grid points on kgrid that are not part of the optical computational domain 
% 
% Hwan Goh, University of Auckland, New Zealand - 23/10/2013
% Last Edited - 20/1/2018: Added normals to boundary elements

%Shortening labels
Nx = Mesh.Nx;
Ny = Mesh.Ny;
dx = Mesh.dx;
dy = Mesh.dy;
Dimensns = Mesh.Dimensns;
I_s = RunOptions.OptI_s;
ExtndBndx = Mesh.ExtndBndx;
ExtndBndy = Mesh.ExtndBndy;

%Checking boundary extensions
if (mod(Mesh.ExtndBndx,5)~=0 || mod(Mesh.ExtndBndy,5)~=0)
    printf(['Boundary extensions must be multiples of 5'])
    return
end

%Defining x_vec and y_vec
x_vec = -Dimensns:dx:Dimensns;
y_vec = -Dimensns:dy:Dimensns;
x_vec = x_vec';
y_vec = y_vec';
%Constructing grid point coordinates
N_GridPts = Nx*Ny;
GridPts = zeros(N_GridPts,2);
if RunOptions.DefineNodesVertically == 1
    for i=1:Nx;
        GridPts(1+(i-1)*Ny:i*Ny,:) = [x_vec(i)*ones(Ny,1),y_vec(1:Ny)];
    end
end
if RunOptions.DefineNodesHorizontally == 1
    for i=1:Ny;
        GridPts(1+(i-1)*Nx:i*Nx,:) = [x_vec(1:Nx),y_vec(i)*ones(Nx,1)];
    end
end
Nodes = GridPts*RunOptions.ScalingOptical;

%Removing nodes on extended boundary
if (ExtndBndx>0 && ExtndBndy>0)
    [minx,~] = find(Nodes(:,1)<-Dimensns(1));
    [maxx,~] = find(Nodes(:,1)>Dimensns(1));
    [miny,~] = find(Nodes(:,2)<-Dimensns(2));
    [maxy,~] = find(Nodes(:,2)>Dimensns(2));
    DomainComplmnt = [minx;maxx;miny;maxy];
    DomainComplmnt = unique(sort(DomainComplmnt));
    Nodes(DomainComplmnt',:)=[];
end

%=========================================================================%
%                           Generating Elements
%=========================================================================%
Elements=unique(sort(delaunay(Nodes(:,1),Nodes(:,2))')','rows');

%Obtaining the node index of the corner nodes
BttmLft = ismember(Nodes,[min(Nodes(:,1)),min(Nodes(:,2))],'rows');
BttmLftNode = find(BttmLft);

BttmRght = ismember(Nodes,[max(Nodes(:,1)),min(Nodes(:,2))],'rows');
BttmRghtNode = find(BttmRght);

TopLft = ismember(Nodes,[min(Nodes(:,1)),max(Nodes(:,2))],'rows');
TopLftNode = find(TopLft);

TopRght = ismember(Nodes,[max(Nodes(:,1)),max(Nodes(:,2))],'rows');
TopRghtNode = find(TopRght);

%=== Timo's Method for ensuring anti-clockwise indexing of element boundary nodes ===%
for ii = 1:size(Elements,1)
  iv = Elements(ii,:); 
  x1 = Nodes(iv,1);
  y1 = Nodes(iv,2);
  rpx = mean(x1);
  rpy = mean(y1);
  for jj = 1:3
    angles(jj) = atan2((y1(jj)-rpy),(x1(jj)-rpx))/pi*180;
  end
  [~, id] = sort(angles);
  Elements(ii,:) = Elements(ii,id);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rotating Corner Elements %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if RunOptions.RotateCornerElementsToOneElement == 1
%=== One Element in Each Corner ===%
[BLElmt,~] = find(Elements==BttmLftNode);
if size(BLElmt,1)==2;
    %Finding the element that shares 2 nodes with the single corner element
    %Finding the elements that share the corner node
    BLElmtNodes = Elements(BLElmt,:);
    [~,BLCornerNodeInElmt1Col] = find(BLElmtNodes(1,:) == BttmLftNode);
    [~,BLCornerNodeInElmt2Col] = find(BLElmtNodes(2,:) == BttmLftNode);
    [~,BLDiffNodeInElmt1Col] = find(ismember(BLElmtNodes(1,:),BLElmtNodes(2,:))==0);
    [~,BLDiffNodeInElmt2Col] = find(ismember(BLElmtNodes(2,:),BLElmtNodes(1,:))==0); 
    [~,BLSharedNodesInElmt2Cols] = find(ismember(BLElmtNodes(2,:),BLElmtNodes(1,:))==1); 
    [~,BLIndex] = find(BLSharedNodesInElmt2Cols ~= BLCornerNodeInElmt2Col);
    BLDiffNodeInElmt1 = BLElmtNodes(1,BLDiffNodeInElmt1Col);
    BLDiffNodeInElmt2 = BLElmtNodes(2,BLDiffNodeInElmt2Col);
    %Rotating BLElmt1 and BLElmt2
    Elements(BLElmt(1),BLCornerNodeInElmt1Col) = BLDiffNodeInElmt2;
    Elements(BLElmt(2),BLSharedNodesInElmt2Cols(BLIndex)) = BLDiffNodeInElmt1;
end

[BRElmt,~] = find(Elements==BttmRghtNode);
if size(BRElmt,1)==2;
    %Finding the elements that share the corner node
    BRElmtNodes = Elements(BRElmt,:);
    [~,BRCornerNodeInElmt1Col] = find(BRElmtNodes(1,:) == BttmRghtNode);
    [~,BRCornerNodeInElmt2Col] = find(BRElmtNodes(2,:) == BttmRghtNode);
    [~,BRDiffNodeInElmt1Col] = find(ismember(BRElmtNodes(1,:),BRElmtNodes(2,:))==0);
    [~,BRDiffNodeInElmt2Col] = find(ismember(BRElmtNodes(2,:),BRElmtNodes(1,:))==0); 
    [~,BRSharedNodesInElmt2Cols] = find(ismember(BRElmtNodes(2,:),BRElmtNodes(1,:))==1); 
    [~,BRIndex] = find(BRSharedNodesInElmt2Cols ~= BRCornerNodeInElmt2Col);
    BRDiffNodeInElmt1 = BRElmtNodes(1,BRDiffNodeInElmt1Col);
    BRDiffNodeInElmt2 = BRElmtNodes(2,BRDiffNodeInElmt2Col);
    %Rotating BRElmt1 and BRElmt2
    Elements(BRElmt(1),BRCornerNodeInElmt1Col) = BRDiffNodeInElmt2;
    Elements(BRElmt(2),BRSharedNodesInElmt2Cols(BRIndex)) = BRDiffNodeInElmt1;
end    

[TLElmt,~] = find(Elements==TopLftNode);
if size(TLElmt,1)==2;
    %Finding the elements that share the corner node
    TLElmtNodes = Elements(TLElmt,:);
    [~,TLCornerNodeInElmt1Col] = find(TLElmtNodes(1,:) == TopLftNode);
    [~,TLCornerNodeInElmt2Col] = find(TLElmtNodes(2,:) == TopLftNode);
    [~,TLDiffNodeInElmt1Col] = find(ismember(TLElmtNodes(1,:),TLElmtNodes(2,:))==0);
    [~,TLDiffNodeInElmt2Col] = find(ismember(TLElmtNodes(2,:),TLElmtNodes(1,:))==0); 
    [~,TLSharedNodesInElmt2Cols] = find(ismember(TLElmtNodes(2,:),TLElmtNodes(1,:))==1); 
    [~,TLIndex] = find(TLSharedNodesInElmt2Cols ~= TLCornerNodeInElmt2Col);
    TLDiffNodeInElmt1 = TLElmtNodes(1,TLDiffNodeInElmt1Col);
    TLDiffNodeInElmt2 = TLElmtNodes(2,TLDiffNodeInElmt2Col);
    %Rotating TLElmt1 and TLElmt2
    Elements(TLElmt(1),TLCornerNodeInElmt1Col) = TLDiffNodeInElmt2;
    Elements(TLElmt(2),TLSharedNodesInElmt2Cols(TLIndex)) = TLDiffNodeInElmt1;
end    

[TRElmt,~] = find(Elements==TopRghtNode);
if size(TRElmt,1)==2;
    %Finding the elements that share the corner node
    TRElmtNodes = Elements(TRElmt,:);
    [~,TRCornerNodeInElmt1Col] = find(TRElmtNodes(1,:) == TopRghtNode);
    [~,TRCornerNodeInElmt2Col] = find(TRElmtNodes(2,:) == TopRghtNode);
    [~,TRDiffNodeInElmt1Col] = find(ismember(TRElmtNodes(1,:),TRElmtNodes(2,:))==0);
    [~,TRDiffNodeInElmt2Col] = find(ismember(TRElmtNodes(2,:),TRElmtNodes(1,:))==0); 
    [~,TRSharedNodesInElmt2Cols] = find(ismember(TRElmtNodes(2,:),TRElmtNodes(1,:))==1); 
    [~,TRIndex] = find(TRSharedNodesInElmt2Cols ~= TRCornerNodeInElmt2Col);
    TRDiffNodeInElmt1 = TRElmtNodes(1,TRDiffNodeInElmt1Col);
    TRDiffNodeInElmt2 = TRElmtNodes(2,TRDiffNodeInElmt2Col);
    %Rotating TRElmt1 and TRElmt2
    Elements(TRElmt(1),TRCornerNodeInElmt1Col) = TRDiffNodeInElmt2;
    Elements(TRElmt(2),TRSharedNodesInElmt2Cols(TRIndex)) = TRDiffNodeInElmt1;
end    
end

%=== Two Elements in Each Corner ===%
if RunOptions.RotateCornerElementsToTwoElements == 1
[BLElmt,~] = find(Elements==BttmLftNode);
if size(BLElmt,1)==1;
    %Finding the element that shares 2 nodes with the single corner element
    BLElmtNodes = Elements(BLElmt,:);
    [BLr BLc] = find(BLElmtNodes~=BttmLftNode);
    [BLShrdElmt1,~] = find(Elements==BLElmtNodes(BLr(1),BLc(1)));
    [BLShrdElmt2,~] = find(Elements==BLElmtNodes(BLr(2),BLc(2)));
    BLShrdElmts = intersect(BLShrdElmt1,BLShrdElmt2);
    H = find(BLShrdElmts ~= BLElmt);
    BLShrdElmt = BLShrdElmts(H);
    %Rotating BLElmt and BLShrdElmt
    BLShrdElmtNodes = Elements(BLShrdElmt,:);
    Elements(BLElmt,:)=[BttmLftNode,BLElmtNodes(2),BLShrdElmtNodes(3)];
    Elements(BLShrdElmt,:)=[BttmLftNode,BLShrdElmtNodes(2),BLShrdElmtNodes(3)];
end

[BRElmt,~] = find(Elements==BttmRghtNode);
if size(BRElmt,1)==1;
    %Finding the element that shares 2 nodes with the single corner element
    BRElmtNodes = Elements(BRElmt,:);
    [BRr BRc] = find(BRElmtNodes~=BttmRghtNode);
    [BRShrdElmt1,~] = find(Elements==BRElmtNodes(BRr(1),BRc(1)));
    [BRShrdElmt2,~] = find(Elements==BRElmtNodes(BRr(2),BRc(2)));
    BRShrdElmts = intersect(BRShrdElmt1,BRShrdElmt2);
    H = find(BRShrdElmts ~= BRElmt);
    BRShrdElmt = BRShrdElmts(H);
    %Rotating BRElmt and BRShrdElmt
    BRShrdElmtNodes = Elements(BRShrdElmt,:);
    Elements(BRElmt,:)=[BRElmtNodes(1),BttmRghtNode,BRShrdElmtNodes(2)];
    Elements(BRShrdElmt,:)=[BttmRghtNode,BRShrdElmtNodes(2),BRShrdElmtNodes(3)];
end    

[TLElmt,~] = find(Elements==TopLftNode);
if size(TLElmt,1)==1;
    %Finding the element that shares 2 nodes with the single corner element
    TLElmtNodes = Elements(TLElmt,:);
    [TLr TLc] = find(TLElmtNodes~=TopLftNode);
    [TLShrdElmt1,~] = find(Elements==TLElmtNodes(TLr(1),TLc(1)));
    [TLShrdElmt2,~] = find(Elements==TLElmtNodes(TLr(2),TLc(2)));
    TLShrdElmts = intersect(TLShrdElmt1,TLShrdElmt2);
    H = find(TLShrdElmts ~= TLElmt);
    TLShrdElmt = TLShrdElmts(H);
    %Rotating TLElmt and TLShrdElmt
    TLShrdElmtNodes = Elements(TLShrdElmt,:);
    Elements(TLElmt,:)=[TLShrdElmtNodes(2),TopLftNode,TLElmtNodes(3)];
    Elements(TLShrdElmt,:)=[TLShrdElmtNodes(1),TLShrdElmtNodes(2),TopLftNode];
end    

[TRElmt,~] = find(Elements==TopRghtNode);
if size(TRElmt,1)==1;
    %Finding the element that shares 2 nodes with the single corner element
    TRElmtNodes = Elements(TRElmt,:);
    [TRr TRc] = find(TRElmtNodes~=TopRghtNode);
    [TRShrdElmt1,~] = find(Elements==TRElmtNodes(TRr(1),TRc(1)));
    [TRShrdElmt2,~] = find(Elements==TRElmtNodes(TRr(2),TRc(2)));
    TRShrdElmts = intersect(TRShrdElmt1,TRShrdElmt2);
    H = find(TRShrdElmts ~= TRElmt);
    TRShrdElmt = TRShrdElmts(H);
    %Rotating TRElmt and TRShrdElmt
    TRShrdElmtNodes = Elements(TRShrdElmt,:);
    Elements(TRElmt,:)=[TRShrdElmtNodes(1),TRElmtNodes(2),TopRghtNode];
    Elements(TRShrdElmt,:)=[TRShrdElmtNodes(1),TRShrdElmtNodes(2),TopRghtNode];
end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generating Boundary Elements %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
[bttmEdg]=EdgElmntsGenerator(Elements,Nodes,min(Nodes(:,2)),2);
[rghtEdg]=EdgElmntsGenerator(Elements,Nodes,max(Nodes(:,1)),1);
[topEdg]=EdgElmntsGenerator(Elements,Nodes,max(Nodes(:,2)),2);
[lftEdg]=EdgElmntsGenerator(Elements,Nodes,min(Nodes(:,1)),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Boundary Structures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bnd_ElmInd = unique(sort([bttmEdg';rghtEdg';topEdg';lftEdg']));
N_Bnd = size(Bnd_ElmInd,1);
N_bttmEdg = size(bttmEdg',1);
N_lftEdg = size(lftEdg',1);
N_rghtEdg = size(rghtEdg',1);
N_topEdg = size(topEdg',1);
[bttmBndNodes,~] = find(Nodes(:,2) == min(Nodes(:,2)));
[leftBndNodes,~] = find(Nodes(:,1) == min(Nodes(:,1)));
[rightBndNodes,~] = find(Nodes(:,1) == max(Nodes(:,1)));
[topBndNodes,~] = find(Nodes(:,2) == max(Nodes(:,2)));
Bnd_NodeInd = unique(sort([bttmBndNodes,leftBndNodes,rightBndNodes,topBndNodes]));
Lght_ElmtInd = Bnd_ElmInd;
Lght_Nelm = size(Lght_ElmtInd,1);

%=========================================================================%
%                                Outputs
%=========================================================================%
Mesh.Dimensns = Dimensns;
Mesh.Nodes = Nodes;
Mesh.Elements=Elements;
Mesh.N_Nodes=size(Nodes,1);
Mesh.N_Elm=size(Mesh.Elements,1);
Mesh.N_Lght=4;
Mesh.Bnd_ElmInd=Bnd_ElmInd;
Mesh.Bnd_NodeInd = Bnd_NodeInd;
Mesh.Lght_ElmtInd=Lght_ElmtInd;
Mesh.Lght_Nelm = Lght_Nelm;
Mesh.I_s = I_s;
Mesh.ExtndBndx = ExtndBndx;
Mesh.ExtndBndy = ExtndBndy;
Mesh.CornerNodes = [BttmLftNode;BttmRghtNode;TopLftNode;TopRghtNode];
Mesh.CornerNodesCoord = [Nodes(BttmLftNode,:);Nodes(BttmRghtNode,:);Nodes(TopLftNode,:);Nodes(TopRghtNode,:)];
if (ExtndBndx>0 && ExtndBndy>0)
    Mesh.DomainComplmnt = DomainComplmnt;
else
    Mesh.DomainComplmnt = [];
end

[Mesh.MinDist Mesh.MaxDist Mesh.Minsize Mesh.Maxsize] = ...
    CheckNodeDistances(Mesh.Elements,Mesh.Nodes,Mesh.N_Elm,Mesh.N_Nodes);

%=========================================================================%
%                      Function: EdgElmntsGenerator
%=========================================================================%

function [EdgElmnts]=EdgElmntsGenerator(Elements,Nodes,CmmnCoordnt,axs)

% EdgElmntsGenerator generates the boundary elements of each edge of a
% rectangular domain
%
% Inputs:
%    Elements - Matrix where row number corresponds to the element number 
%               and the entries of the row are the vertices of the element
%    Nodes - Coordinates of the node
%    CmmnCoordnt - The coordinate that all edge nodes share
%    axs - The axis on which the shared coordinate lies; X-axis=1,Y-axis=2
%
% Outputs:
%    EdgElmnts - Vector listing the elements relevant to the edge]
%
% Hwan Goh, University of Auckland, New Zealand - 26/10/2013

        %===================================================%
        %   Finding elements with nodes on the boundary     %
        %===================================================%

[row col] = find(Nodes == CmmnCoordnt);
A=[row col];
EdgNodes = [];
for i=1:size(A,1);
    if A(i,2)==axs;
        EdgNodes(i)=A(i,1); %will have entries of zero since we exclude 
                             %when the x-coordinate equals min(Nodes,2)
    end
end
G=find(EdgNodes); %vector of the location of non-zero entries of bNodes
ValidNodes=zeros(size(G,2),1);
for i=1:size(G,2);
    ValidNodes(i) = EdgNodes(G(i)); %obtains the values of the non-zero entries
end
H = {};
for i=1:size(ValidNodes);
    [row ~] = find(Elements == ValidNodes(i)); %
    H{i} = row;
end
Elmnts = [];
for i=1:size(H,2);
    Elmnts = [Elmnts; cell2mat(H(i))];
end
Elmnts = unique(Elmnts);
N_Bnd = size(Elmnts,1);

      %==========================================================%
      %   %Removing elements without an edge on the boundary     %
      %==========================================================%

for i=1:N_Bnd;
    ElmntNodes = zeros(3,2);
    for j=1:3;
        ElmntNodes(j,:) = Nodes(Elements(Elmnts(i),j),:);
    end
    [rows col] = find(ElmntNodes(:,axs) == CmmnCoordnt);
    if size(rows,1) < 2;
        Elmnts(i)=0;
    end
end

validElmnts = find(Elmnts);
EdgElmnts = [];
for i=1:size(validElmnts)
    EdgElmnts(i) = Elmnts(validElmnts(i));
end

%=========================================================================%
%                      Function: CheckNodeDistances
%=========================================================================%

function [minDist maxDist minEDist maxEDist]=CheckNodeDistances(Elements,Nodes,N_Elm,N_Nodes)

% CheckNodeDistances calculates some distance values for a given mesh. Main
% use of this is to check for needle elements
%
% Inputs:
%   Elements - Matrix where row number corresponds to the elementnumber and the entries of the row are the vertices of the element
%   Nodes - Coordinates of nodes
%   N_Elm - Number of elements
%   N_Nodes - Number of Nodes
%
% Outputs:
%   minDist - minimum distance between nodes
%   maxDist - maximum distance between nodes
%   minEDist - minimum size of an edge
%   maxEDist - maximum distance of an edge
% 
% Hwan Goh, University of Auckland, New Zealand - 14/11/2013

distances = zeros(N_Nodes,N_Nodes);

for i=1:N_Nodes;
    for j=1:N_Nodes;
        distances(j,i) = sqrt((Nodes(i,1)-Nodes(j,1))^2 + (Nodes(i,2)-Nodes(j,2))^2);
    end
end

minDist = min(distances(distances~=0));
maxDist = max(distances(:));
    
Edistances = zeros(N_Elm,3);

for i=1:N_Elm;
    Edistances(i,1) = sqrt((Nodes(Elements(i,1),1)-Nodes(Elements(i,2),1))^2 + (Nodes(Elements(i,1),2)-(Nodes(Elements(i,2),2)))^2);
    Edistances(i,2) = sqrt((Nodes(Elements(i,1),1)-Nodes(Elements(i,3),1))^2 + (Nodes(Elements(i,1),2)-(Nodes(Elements(i,3),2)))^2);
    Edistances(i,3) = sqrt((Nodes(Elements(i,2),1)-Nodes(Elements(i,3),1))^2 + (Nodes(Elements(i,2),2)-(Nodes(Elements(i,3),2)))^2); 
end

minEDist = min(Edistances(:));
maxEDist = max(Edistances(:));
