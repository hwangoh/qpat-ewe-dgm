function [fmu_a,PrmtrsI,phi,A,R,b_pattern,p0] = DA_FEM2D_GNFwrd(RunOptions,DataVrblsOptical,MeshI,mu_a,mu_s,PrmtrsPrp,R,b_pattern)

% DA_FEM2D_GNFwrd computes the optical forward problem of ODT using FEM to
% obtain the initial pressure p0.
%
% Inputs:
%   RunOptions: Not used here, but required as an input so that
%               GNIterationsQPAT is a general function
%   DataVrblsOptical:
%      NumberofSensors - Number of sensors
%   MeshI:
%      MeshI.N_Nodes - Number of nodes
%   mu_a: current absorbtion coefficient. Note that these fields were deliberately separated out from the 
%                                         structure PrmtrsI for the purpose of using fminsearch when line search is implemented
%   mu_s: selected scattering coefficient
%   PrmtrsPrp:
%      g - mean of the cosine of the scattering angle, used for computing the diffusion coefficient
%   R: Boundary Matrix
%   b: Lightsource Vector
%
% Outputs:
%   fmu_a: sensory data, which is a function of mu_a  
%   PrmtrsI:
%      PrmtrsI.kappa_elmts - Elementwise organisation of kappa values
%      PrmtrsI.mu_a_elmts - Elementwise organisation of mu_a values
%   phi: light fluence
%   A: S+M+R which depends on the current mu_a
%   R - Boundary Matrix, doesn't depend on mu_a so is an ouput so that it can be used later
%   b_pattern - Lightsource Vector for current illumination pattern, so is an ouput so that it can be rused later 
%
% Hwan Goh 23/12/2015, University of Auckland, New Zealand

[p0,PrmtrsI,phi,A,R,b_pattern] = GNFrwdODTFEM(RunOptions,DataVrblsOptical,MeshI,mu_a,mu_s,PrmtrsPrp,R,b_pattern); %How is this different from the overall function???
fmu_a = p0(:);