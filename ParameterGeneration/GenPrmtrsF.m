function [prmtr, R]=GenPrmtrsF(Mesh,bckgrnd,scalar,maximum)

% GenPrmtrsF generates the scattering and absorption coefficients for QPAT forward problem in a circular domain. 
% Specifically, Gaussian blobs are generated.
%
%  Inputs:
%    Mesh - The linear mesh
%    bckgrnd - Background values for the parameter.
%    scalar - Controls the size of the parameter to reflect experimental values
%    maximum - Creates a maximum size of the parameter to reflect experimental values 
%
%  Outputs:
%    prmtr - N_nodes  by 1 vector of draws
%    R - [min(prmtr) max(prmtr)]
%
% Hwan Goh, University of Auckland, New Zealand 15/06/2013
% Adapted from P.J. Hadwin, University of Auckland, New Zealand 25/09/2012

N=randi([2 3]);
Ix=[-7;5.5;3.5]; Iy=[0.5;6;-6.75];
Center=[Ix(1:N)+0.5*randn(N,1) Iy(1:N)-abs(randn(N,1))];
Cor=0.5*(randi([0 100],N,2)+randi([0 10],N,2))/5;
R(1)=max([1.2*rand(1)/5 0.0075]);
R(2)=R(1)+rand(1);

N=Mesh.N_Nodes;
prmtr=zeros(N,1)+bckgrnd;
for ii=1:size(Center,1)
    Mu=Center(ii,:)';
    Gam=diag(Cor(ii,:));
    Bubb=gaussian(Mesh.Nodes(1:N,:)',Mu,Gam)';
    BuB = 1/(2*max(Bubb))*Bubb;
    BuB = scalar*BuB;
    prmtr = prmtr + BuB;
end




