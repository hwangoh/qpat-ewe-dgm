function [prmtrD, R]=GenPrmtrsTest(Mesh,scalar)

% GenPrmtrsITest generates the initial pressure. 
%
%  Inputs:
%   RunOptions:
%      ScalingOptical - Relates to the units the optical mesh is in.
%                       Default is metres.
%    Mesh - The data mesh
%    MeshI - The inverse mesh
%    bckgrnd - Background values for the parameter.
%    scalar - Controls the size of the parameter to reflect experimental values
%    PrmtrType - Whether generating absorption or scattering coefficient
%
%  Outputs:
%    prmtr - N_nodes  by 1 vector of draws
%    R - [min(prmtr) max(prmtr)]
%
% Hwan Goh, University of Auckland, New Zealand 26/03/2018
% Adapted from P.J. Hadwin, University of Auckland, New Zealand 25/09/2012

Center = [0.005,0.005];
Cor = [3.700000000000000,3.700000000000000];
R = [0.117420633216005,0.455140043037382];

N=Mesh.N_Nodes;
prmtrD=zeros(N,1);

for ii=1:size(Center,1)
    Mu = Center(ii,:)';
    Gam = diag(Cor(ii,:))*10^-6;
    Bubb = GenPrmtrsGaussian(Mesh.Nodes(1:N,:)',Mu,Gam)';
    BuB = 1/(2*max(Bubb))*Bubb;
    BuB = scalar*BuB;
    prmtrD = prmtrD + BuB;
end