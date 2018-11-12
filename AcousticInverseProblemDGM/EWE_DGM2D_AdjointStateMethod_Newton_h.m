function [Dirctn] = EWE_DGM2D_AdjointStateMethod_Newton_h(RunOptions,SensorsI,FrwdFunction,hRecon,alphaRecon,Gradient,Gradient_alpha,Gradient_InitialGuess,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11AdjTime0,T22AdjTime0,Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,CurrentIteration,hRelError,PLOT)

% EWE_DGM2D_AdjointStateMethod_Newton_h.m computes the Newton or Gauss-Newton direction
% when constrained by the elastic wave equation discretized by the
% discontinuous Galerkin method
%
% Inputs:
%   RunOptions:
%      NumberofTimeSteps - Number of time steps to be computed
%      BC - Boundary conditions
%      EWE_LS_ExponentialTransform - Use exponential transform
%   SensorsI: 1 by number of sensors array containing 3 by 1 cells 
%      id - Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%      xy - coordinates of the sensor
%      l_iatsensor - 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%   FrwdFunction - Script for the forward function of h
%   hRecon - Reconstructed initial pressure. This is required when positivity transform is implemented, because the mapping is linearized around the MAP estimate
%   alphaRecon - Reconstructed hyperprior parameter
%   Gradient - Gradient of functional computed with adjoint-state method
%   Gradient_alpha - Gradient of alpha computed with adjoint-state method
%   Gradient_InitialGuess - Gradient of the functional at the initial guess
%   L_ET11 - Error Model for T11
%   L_ET22 - Error Model for T22
%   L_ET12 - Error Model for T12
%   L_Evx - Error Model for vx
%   L_Evy - Error Model for vy
%   TiiAdjTime0 - Adjoint variables, to be used for computing Newton direction
%   Prior:
%      L_pr - Regularization operator
%      Exp_h - Expected value
%   MeshIElm: Number of elements by 3 array containing the indices of the nodes in each element
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshINodes: 2 by Number of nodes array for the coordinates of the nodes on the inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   pinfoI - Information regarding p-refinement for p-nonconforming meshes on DGM inversion mesh
%   PrecompIntrplteObjectsI - Objects that depend on the inverse Mesh nodes. May have been computed earlier and so can 
%                             be called here to avoid repeating computations. Set to 0 in function call if you
%                             want to compute new objects for a new mesh.
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI - Number of nodes per element on DGM inversion mesh
%   KI - Number of elements on DGM inversion mesh
%   dt - Size of time step, need negative of this to do reverse time stepping
%   PLOT: To plot or not to plot, that is the question
%
% Outputs:
%   Dirctn - Steepest descent seach direction
%
% Hwan Goh 27/6/2018, University of Auckland, New Zealand (In Malaysia after one month trip but going home soon!)


%% =======================================================================%
%                       Conjugate Gradient Inputs
%=========================================================================%
Norm_Gradient_IniGuess = norm(Gradient_InitialGuess,2);
Norm_Gradient_hRecon = norm(Gradient,2);
tol = min(0.5,sqrt(Norm_Gradient_IniGuess/Norm_Gradient_hRecon)*Norm_Gradient_hRecon);
maxit = 100;

if RunOptions.EWE_LS_SigmoidTransform == 1 || RunOptions.EWE_LS_ExponentialTransform == 1
    error('EWE_DGM2D_AdjointStateMethod_Newton_h has not been set up for these positivity constraints')
end

%% =======================================================================%
%                     Conjugate Gradient Iterations
%=========================================================================%
disp(' ')
disp('------------------------------------------------------')
disp('Using bicg to compute Newton search direction ')
disp('------------------------------------------------------')


% Dirctn = bicg(@(hhat,transp_flag) EWE_DGM2D_AdjStateHessianAction(RunOptions,SensorsI,FrwdFunction,hhat,hRecon,...
%                                                                    L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
%                                                                    T11AdjTime0,T22AdjTime0,...
%                                                                    Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,...
%                                                                    pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,...
%                                                                    CurrentIteration,hRelError,PLOT,transp_flag),...
%               -Gradient,tol,maxit,Prior.L_pr',Prior.L_pr,InitialGuess);

if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 0;
alphaRecon = 1;
Dirctn = pcg(@(hhat) EWE_DGM2D_AdjStateHessianAction(RunOptions,SensorsI,FrwdFunction,hhat,hRecon,alphaRecon,...
                                                                   L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                                   T11AdjTime0,T22AdjTime0,...
                                                                   Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,...
                                                                   pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,...
                                                                   CurrentIteration,hRelError,PLOT,'notransp'),...
              -Gradient,tol,maxit,Prior.L_pr',Prior.L_pr);
end
if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 1;
HyperL_pr = sparse(MeshIN_Nodes+1,MeshIN_Nodes+1);
HyperL_pr(MeshIN_Nodes+1,MeshIN_Nodes+1) = Prior.Hyperprior_Std_alpha^(-1);
HyperL_pr(1:MeshIN_Nodes,1:MeshIN_Nodes) = alphaRecon*Prior.L_pr;
Dirctn = pcg(@(halphahat) EWE_DGM2D_AdjStateHessianAction_Hyperprior(RunOptions,SensorsI,FrwdFunction,halphahat,hRecon,alphaRecon,...
                                                                      L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                                      T11AdjTime0,T22AdjTime0,...
                                                                      Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,...
                                                                      pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,...
                                                                      CurrentIteration,hRelError,PLOT,'notransp'),...
              -[Gradient;Gradient_alpha],tol,maxit,HyperL_pr',HyperL_pr);
end
          
          