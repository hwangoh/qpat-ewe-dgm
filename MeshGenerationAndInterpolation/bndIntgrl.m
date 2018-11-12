function [R_int,b_int] = bndIntgrl(Mesh,VerticesCoords,i,I_s,Lght_ElmtInd)

% bndIntgrl computes the integral of the edge of an element lying on the boundary
% as well as the integral of the edge of an element lying under a light source on the boundary.
%
% Inputs:
%    Mesh:
%      Dimensns - Dimensions of the Mesh, 2x1 for rectangular and 1x1 for
%                 circular
%    E - a 3 by 2 matrix where each row contains coordinates of the
%        vertices for the element
%    i - Current element number
%    I_s - The diffuse boundary current
%    Lght_ElmtInd - A (lght_Nelm x L) matrix where each column represents the set of elements under one light source
%
% Outputs:
%    R_int - Matrix representing the values of the integrals for the R matrix
%    b_int - Value of the integral for the b vector
%
% Hwan Goh 06/06/2013, University of Auckland, New Zealand
% adapted from P.J. Hadwin 07/07/2010, University of Auckland, New Zealand
% Hwan Goh 27/01/2017, Almost four years later! Merged this with bndIntgrlLight.m for efficiency
%
% Notes: "i==Lght_ElmtInd" creates a matrix of the same dimensions as Lght_ElmtInd, where if an entry of Lght_ElmtInd = i then the corresponding entry of i==Lght_ElmtInd is a 1, else a zero. There should only be one, and "find" will extract the linear index of the entry. Then "any" returns one if non-zero 
%        light_Ind = find(any(i==Lght_ElmtInd)); %Returns the light source index the element is under. "any(i==Lght_ElmtInd)" returns a (1xN_lghtElm) matrix of ones and zeroes with a one if a column of Lght_ElmtInd contains an entry equal to i. Then "find" returns the column number of "any(i==Lght_ElmtInd)".

%=== Rectangular Mesh ===%
if size(Mesh.Dimensns,2) == 2;
    H = {};
    H{1} = find(VerticesCoords(:,2)==min(Mesh.Nodes(:,2)));
    H{2} = find(VerticesCoords(:,1)==min(Mesh.Nodes(:,1)));
    H{3} = find(VerticesCoords(:,1)==max(Mesh.Nodes(:,1)));
    H{4} = find(VerticesCoords(:,2)==max(Mesh.Nodes(:,2)));

    for j=1:4;
        if size(cell2mat(H(j)),1)==2
            EdgNode = cell2mat(H(j));
        end
    end
    %=== R matrix ==%
    Jac = norm(VerticesCoords(EdgNode(1),:) - VerticesCoords(EdgNode(2),:),2);
    if ~any(find(EdgNode==3));
        R_int = Jac*(2/pi)*[1/3 1/6 0; 1/6 1/3 0; 0 0 0]; %including zeroes to represent the integrals of the nodes not on the boundary.
    end
    if ~any(find(EdgNode==2));
        R_int = Jac*(2/pi)*[1/3 0 1/6; 0 0 0; 1/6 0 1/3];
    end
    if ~any(find(EdgNode==1));
        R_int = Jac*(2/pi)*[0 0 0; 0 1/3 1/6; 0 1/6 1/3];
    end 
    %=== b vector ==%
    if any(find(i==Lght_ElmtInd)) %Check if element i is also under a light source
        if ~any(find(EdgNode==3));
            b_int = Jac*I_s*[1;1;0]; %including zeroes to represent the integrals of the nodes not on the boundary.
        end
        if ~any(find(EdgNode==2));
            b_int = Jac*I_s*[1;0;1];
        end
        if ~any(find(EdgNode==1));
            b_int = Jac*I_s*[0;1;1];
        end
    end
%=== Circular Mesh ===%
elseif size(Mesh.Dimensns,2) == 1; %If mesh is circular
    %=== R matrix ==%
    Jac = norm(VerticesCoords(1,:) - VerticesCoords(2,:),2);
    R_int = Jac*(2/pi)*[1/3 1/6 0; 1/6 1/3 0; 0 0 0]; %including zeroes to represent the integrals of the nodes not on the boundary.
    %=== b vector ==%
    if any(find(i==Lght_ElmtInd))
        b_int = Jac*I_s*[1;1;0];
    end
end