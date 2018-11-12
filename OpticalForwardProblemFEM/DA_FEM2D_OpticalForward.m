function [phi,S,M,R,b,b_pattern]=DA_FEM2D_OpticalForward(RunOptions,Mesh,Prmtrs)

% DA_FEM2D_OpticalForward.m is the main script file for the optical forward problem
% using the diffusion approximation and the finite element method.
%
% Inputs:
%   Mesh:
%        Mesh.Nodes - Coordinates of nodes
%        Mesh.Elements - Matrix where row number corresponds to the element number and the entries of the row are the vertices of the element
%        Mesh.N_Nodes - Number of nodes
%        Mesh.N_Elm - Number of elements
%        Mesh.N_Lght - Number of light sources
%        Mesh.Bnd_ElmInd - A column vector listing the element indices for elements with an edge on the boundary
%        Mesh.Bnd_NodeInd - A column vector listing the node indices for the nodes on the boundary
%        Mesh.Lght_ElmtInd - A (lght_Nelm x L) matrix where each column represents the set of elements under one light source
%        Mesh.Lght_Width - Vector of computational light source widths
%        Mesh.Lght_Nelm - Number of elements under a light source
%        Mesh.I_s - Diffuse boundary condition   
%   Prmtrs:
%        Prmtrs.kappa_elmts - Elementwise organisation of kappa values
%        Prmtrs.mu_a_elmts - Elementwise organisation of mu_a values
%   PLOTOptical - To plot or not to plot, that is the question
%
%  Outputs:
%   S - Stiffness Matrix
%   M - Mass Matrix
%   R - Boundary Matrix
%   b - Lightsource Vector
%   phi - Fluence vector
%   b_pattern - N_Patterns by N_Nodes Matrix where each row is a lightsource vector with deactivated light sources
%
% Hwan Goh, 04/05/2013, University of Auckland, New Zealand

%% =======================================================================%
%                       Piecewise Linear Case
%=========================================================================%
%Constructing the system matrices for the piecewise linear case
[S, M, R, b] = pwlFEMMatrices(Mesh,Prmtrs);

%% =======================================================================%
%                              Solution
%=========================================================================%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Implementing Illumination Patterns %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ActiveSources_Nodes] = ActivateSources(RunOptions,Mesh); %1 by N_Nodes vector containing a 1 when a node corresponds with a node of an element under an active light source.
b_pattern = zeros(1,Mesh.N_Nodes);
b_pattern = ActiveSources_Nodes.*b;   

%%%%%%%%%%%%%%%
%%% Solving %%%
%%%%%%%%%%%%%%%
A=S+M+R;
% L=chol(A);
phi=A\b_pattern';
