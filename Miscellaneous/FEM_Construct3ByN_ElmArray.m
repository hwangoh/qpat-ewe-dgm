function ThreebyN_ElmArray = FEM_Construct3ByN_ElmArray(Nodes,Elements,Prmtr)


% FEM_ConstructN_ElmBy3Array constructs a 3 by N_Elm array
%
%Inputs: 
%   Nodes - Coordinates of Nodes [x y z] needs to be at least 2-D
%   Elements - Nodes of Elements, assumed to be linear Basis
%   Prmtr - 1 by N_Nodes vector where each entry reperesents the
%           coefficient's value at that node
%
%Outputs:
%   ThreebyN_ElmArray - The N_Elm by 3 array
%
%Hwan Goh, University of Auckland 12/10/2017 
%Adapted from P. J. Hadwin, University of Auckland, New Zealand, 24/11/2011

%Input Manipulation 
N=size(Prmtr,2);
g=Nodes; 
H=Elements;
[I1, I2]=size(H);

X=reshape(g(H(:,1:I2),1),I1,I2)';
Y=reshape(g(H(:,1:I2),2),I1,I2)';
ThreebyN_ElmArray=zeros([size(X) N]);

for ii=1:N
    ThreebyN_ElmArray(:,:,ii)=reshape(Prmtr(H(:,1:I2),ii),I1,I2)';
end