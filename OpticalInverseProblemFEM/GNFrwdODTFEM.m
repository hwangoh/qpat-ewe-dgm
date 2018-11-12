function [p0,PrmtrsI,phi,A,R,b_pattern] = GNFrwdODTFEM(RunOptions,DataVrbls,MeshI,mu_a,mu_s,PrmtrsPrp,R,b_pattern)

% GNFrwdODTFEM computes the optical forward problem of QPAT using FEM to
% obtain the initial pressure p0 to be used for Gauss-Newton iterations.
%
% Inputs:
%   DataVrbls:
%      Gruneisen - Gruneisen parameter
%   Mesh:
%      Mesh.N_Nodes - Number of nodes  
%   PrmtrsI:
%      Prmtrs.kappa_elmts - Elementwise organisation of kappa values
%      Prmtrs.mu_a_elmts - Elementwise organisation of mu_a values
%   PrmtrsPrp:
%      g - mean of the cosine of the scattering angle, used for computing the diffusion coefficient
%   R - Boundary Matrix, doesn't depend on mu_a so can be called and reused
%       here if already computed previously
%   b_pattern - Lightsource Vector for current illumination pattern, 
%               so can be called and reused here if already computed previously
%
% Outputs:
%   p0 - Initial Pressure [mm]
%   PrmtrsI:
%          kappa - diffusion coefficient
%   phi - Light fluence
%   A - S+M+R
%   R - Boundary Matrix, doesn't depend on mu_a so only needs to be
%       computed here once
%   b_pattern - Lightsource Vector for current illumination pattern, 
%               doesn't depend on mu_a so only needs to be computed here once
%
% Hwan Goh 29/12/2015, University of Auckland, New Zealand

%absorption := mu_a
mu_a_elmts = reshapeNodes(MeshI,mu_a); %reshaping required for ODTForward
%diffusion := kappa
kappa = zeros(MeshI.N_Nodes,1);
for ii=1:MeshI.N_Nodes;
    kappa(ii) = 1/(2*(mu_a(ii) + (1-PrmtrsPrp.g)*mu_s(ii)));
end
kappa_elmts = reshapeNodes(MeshI,kappa); %reshaping required for ODTForward. Also no need to scale by 10^-3 since mu_a and mu_s are already in m^-1 when kappa is formed

PrmtrsI.mu_a = mu_a;
PrmtrsI.mu_a_elmts = mu_a_elmts;
PrmtrsI.mu_s = mu_s;
PrmtrsI.kappa = kappa;
PrmtrsI.kappa_elmts = kappa_elmts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Computations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Storing Outputs of Each Iteration ===%
H = sparse(MeshI.N_Nodes,1);
p0 = sparse(MeshI.N_Nodes,1);
phi = sparse(MeshI.N_Nodes,1);

%=== ODT Solution and System Matrices ===%
if ischar(R)==1 && ischar(b_pattern)==1; %If first time generating R and b_pattern
    clear R b_pattern
    [phi,S,M,R,~,b_pattern] = DA_FEM2D_OpticalForward(RunOptions,MeshI,PrmtrsI);
    A = S+M+R;
end
if ischar(R)~=1 && ischar(b_pattern)~=1; %If using previously generated R and b_pattern
    [S,M] = pwlMatricesUpdate(MeshI,PrmtrsI);
    A = S+M+R;
%     L = chol(A);
%     phi = L\(L'\b_pattern');
    phi = A\b_pattern';
end

%=== Computing fmu_a ===%
H = (PrmtrsI.mu_a.*phi);
p0 = DataVrbls.Gruneisen*H;
    
       