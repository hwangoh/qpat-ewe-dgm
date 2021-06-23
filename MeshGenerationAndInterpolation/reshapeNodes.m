function C = reshapeNodes(Mesh,Prmtr)

% reshapeNodes reorganizes the Prmtr vector into an elementwise
% representation of its values.
%
% Inputs:
%   Mesh - FEM Mesh Structure
%        Nodes - Coordinates of Nodes [x y z] needs to be at least 2-D
%        Elements - Nodes of Elements, assumed to be linear Basis
%   Prmtr - 1 by N_Nodes vector where each entry reperesents the
%        coefficient's value at that node
%
% Outputs:
%   C = 3 by N_elm matrix where each column represents an element and each
%       entry of the column represents the value of the parameter at the
%       node
%
% Hwan Goh, University of Auckland, New Zealand 01/09/2013
% Adapted from P. J. Hadwin, University of Auckland, New Zealand
%    24/11/11 - Original

N=size(Prmtr,2);
g=Mesh.Nodes; 
H=Mesh.Elements;
[I1 I2]=size(H);

X=reshape(g(H(:,1:I2),1),I1,I2)';
Y=reshape(g(H(:,1:I2),2),I1,I2)';
C=zeros([size(X) N]);

for ii=1:N
    C(:,:,ii)=reshape(Prmtr(H(:,1:I2),ii),I1,I2)';
end