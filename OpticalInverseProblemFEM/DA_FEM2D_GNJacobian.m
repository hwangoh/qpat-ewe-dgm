function JacODTFrwd = DA_FEM2D_GNJacobian(MeshI,PrmtrsI,phi,A)

% DA_FEM2D_GNJacobian constructs the Jacobian of the ODT optical
% forward function to be used for Gauss-Newton iterations.
%
% Inputs:
%   MeshI:
%        MeshI.N_Nodes - Number of nodes
%   N_Illum: Number of illumination patterns
%   PrmtrsI:
%          Prmtrs.kappa - N_Nodes by 1 vector for the diffusion coefficient
%          Prmtrs.mu_a - N_Nodes by 1 vector for the absorption coefficient
%   phi - current light fluence depending on the current mu_a
%   A - current matrix A = S + M + R depending on the current mu_a
%   ~ - JacIntrplte isn't required here
%
% Outputs:
%   JacFrwdODT - Jacobian of Jacobian of the DGM acoustic forward function
%
% Hwan Goh 29/12/2015, University of Auckland, New Zealand

JacODTFrwd = ConstructJacmuaphiAbs(MeshI,PrmtrsI,phi,A);


