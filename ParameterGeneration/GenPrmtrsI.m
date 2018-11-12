function [prmtrD, prmtrI, R]=GenPrmtrsI(RunOptions,MeshD,MeshI,bckgrnd,scalar,PrmtrType)

% GenPrmtrsI generates the scattering and absorption coefficients for QPAT inverse problem in a circular domain. 
% Specifically, Gaussian blobs are generated. Note that this works the same
% as GenPrmtrs for the forward problem, only it is run twice; once for the
% forward mesh and once for the inverse mesh.
%
%  Inputs:
%   RunOptions:
%      ScalingOptical - Relates to the units the optical mesh is in.
%                       Default is metres.
%    MeshD - The data mesh
%    MeshI - The inverse mesh
%    bckgrnd - Background values for the parameter.
%    scalar - Controls the size of the parameter to reflect experimental values
%    PrmtrType - Whether generating absorption or scattering coefficient
%
%  Outputs:
%    prmtr - N_nodes  by 1 vector of draws
%    R - [min(prmtr) max(prmtr)]
%
% Hwan Goh, University of Auckland, New Zealand 22/08/2013
% Adapted from P.J. Hadwin, University of Auckland, New Zealand 25/09/2012

if RunOptions.GeneratePrmtrs == 1 %Generate parameters
    N=randi([2 3]);
    Ix=[-7;5.5;3.5]; Iy=[0.5;6;-6.75]; %Choose original coordinates of center of parameters
    Cor=0.5*(randi([0 100],N,2)+randi([0 10],N,2))/5; %Controls the size of the parameters
    Center=[Ix(1:N)+0.5*randn(N,1) Iy(1:N)-abs(randn(N,1))]*10^-3*RunOptions.ScalingOptical; %*10^-3 because the grid is in metres by default
    R(1)=max([1.2*rand(1)/5 0.0075]);
    R(2)=R(1)+rand(1);
end
if RunOptions.GeneratePrmtrs == 0 %Load pre-generated parameters
    if strcmp(PrmtrType,'absorption') == 1
        if  RunOptions.GenerateAwayFromBoundary == 1
%             Center = [-0.007031527436595,3.758556517836881e-04;0.005857371451913,0.004510302392215;0.003397516970850,-0.004159034489800];
%             Cor = [1.00000000000000,1.00000000000000;2.200000000000000,1.900000000000000;2.00000000000000,3.800000000000000];
%             R = [0.230278182334297,0.886018881490884];
            Center = [0.003,0.003];
            Cor = [2.00000000000000,2.00000000000000];
            R = [0.230278182334297,0.886018881490884];
        else
            Center = [-0.007578200827832,-4.642294226316275e-04;0.005233221445342,0.005479939898545;0.002498682132058,-0.006770027851643];
            Cor = [6.300000000000000,8.800000000000000;5.900000000000000,6.300000000000000;1.400000000000000,3.700000000000000];
%             Center = [-0.004078200827832,-4.642294226316275e-04;0.005233221445342,0.005479939898545;0.002498682132058,-0.006770027851643];
%             Cor = [4.300000000000000,4.800000000000000;5.900000000000000,6.300000000000000;1.400000000000000,3.700000000000000];
            R = [0.117420633216005,0.455140043037382];
        end
    end
    if strcmp(PrmtrType,'scattering') == 1
        Center = [-0.006533135918664,4.709942362912737e-04;0.005675160500678,0.005817547832494];
        Cor = [2.900000000000000,5;3.400000000000000,5.100000000000000];
        R = [0.056948059145165,0.515796887325096];
    end
end

%=========================================================================%
%                             Forward Mesh
%=========================================================================%
N=MeshD.N_Nodes;
prmtrD=zeros(N,1)+bckgrnd;
blobDiameter = RunOptions.CentralGaussianSize;

%=== Randomly Generated ===%
if RunOptions.GenerateCentralGaussian == 0;
for ii=1:size(Center,1)
    Mu = Center(ii,:)';
    Gam = diag(Cor(ii,:))*10^-6*(RunOptions.ScalingOptical^2); %remember the covariance matrix is E[(X - E[X])(X - E[X])'] and since (X - E[X]) is scaled by 10^-3 then we should have 10^-6
    Bubb = GenPrmtrsGaussian(MeshD.Nodes(1:N,:)',Mu,Gam)';
    BuB = 1/(2*max(Bubb))*Bubb;
    BuB = scalar*BuB;
    prmtrD = prmtrD + BuB;
end
end

%=== Central Gaussian ===%
if RunOptions.GenerateCentralGaussian == 1;
    xAxisNodes = find(MeshD.Nodes(:,2)==0);
    CenterNodeNum = xAxisNodes(1+floor((size(xAxisNodes,1)/2)));
    Center = MeshD.Nodes(CenterNodeNum,:);
    Mu = Center';
    Gam = [blobDiameter,0;0,blobDiameter]*(RunOptions.ScalingOptical^2);
    Bubb = GenPrmtrsGaussian(MeshD.Nodes(1:N,:)',Mu,Gam)';
    BuB = 1/(2*max(Bubb))*Bubb;
    BuB = scalar*BuB;
    prmtrD = prmtrD + BuB;
end

%=========================================================================%
%                             Inverse Mesh
%=========================================================================%
N=MeshI.N_Nodes;
prmtrI=zeros(N,1)+bckgrnd;

%=== Randomly Generated ===%
if RunOptions.GenerateCentralGaussian == 0;
for ii=1:size(Center,1)
    Mu=Center(ii,:)';
    Gam=diag(Cor(ii,:))*10^-6*(RunOptions.ScalingOptical^2);
    Bubb=GenPrmtrsGaussian(MeshI.Nodes(1:N,:)',Mu,Gam)';
    BuB = 1/(2*max(Bubb))*Bubb;
    BuB = scalar*BuB;
    prmtrI = prmtrI + BuB;
end
end

%=== Central Gaussian ===%
if RunOptions.GenerateCentralGaussian == 1;
    xAxisNodes = find(MeshI.Nodes(:,2)==0);
    CenterNodeNum = xAxisNodes(1+floor((size(xAxisNodes,1)/2)));
    Center = MeshD.Nodes(CenterNodeNum,:);
    Mu = Center';
    Gam = [blobDiameter,0;0,blobDiameter]*(RunOptions.ScalingOptical^2);
    Bubb = GenPrmtrsGaussian(MeshI.Nodes(1:N,:)',Mu,Gam)';
    BuB = 1/(2*max(Bubb))*Bubb;
    BuB = scalar*BuB;
    prmtrI = prmtrI + BuB;
end    
