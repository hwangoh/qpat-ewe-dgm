function HAOutput = EWE_DGM2D_AdjStateHessianAction_Hyperprior(RunOptions,SensorsI,FrwdFunction,halphahat,hRecon,alphaRecon,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11AdjTime0,T22AdjTime0,Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,CurrentIteration,hRelError,PLOT,transp_flag)

% AdjStateHessianAction_Hyperprior compute the action of the Hessian on the vector via
% the adjoint-state method with a hyperprior implemented. This circumvents explicitly constructing the
% Hessian.
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
%   hhat - Input vector of size N_Nodes with respect to FEM Mesh
%   hRecon - Reconstructed initial pressure. This is required when positivity transform is implemented, because the mapping is linearized around the current MAP estimate
%   alphahat - Input scalar representing hyperprior parameter
%   alphaRecon - Reconstructed hyperprior parameter alpha
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
%   PrecompIntrplteObjects: (These are used for when interpolation on the same meshes are reused as well as computing the Jacobian
%                               of the intrplteFEM2DGM function)
%                 DT - Delaunay triangulated mesh required to use the pointLocation function. Therefore corresponds with the
%                      PointLocs vector.
%                 PrecompCoeffMatricesInv - N_Elm by 1 cell where each cell contains
%                                the required matrix for computing the coefficients for that element's
%                                interpolating function of two variables. These matrices will be required for
%                                computing the Jacobian of this interpolation function
%                 PointLocations - Np by 1 vector containing element numbers (with reference to DT) each grid point lies in. 
%                                This will be required for computing the Jacobian of this
%                                interpolation function.
%                 PrecompIntrplteCoeff - 3 by N_Elm matrix where each column contains the three coefficients required for
%                                interpolation over an element
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI - Number of nodes per element on DGM inversion mesh
%   KI - Number of elements on DGM inversion mesh
%   dt - Size of time step, need negative of this to do reverse time stepping
%   PLOT: To plot or not to plot, that is the question
%   CurrentIteration: Current iteration of whatever you're using the Hessian to compute. Just set as 'Not Required' to not display
%   hRelError: Indicator for when the Hessian is used to compute line search iterations. Just set as 'Not Required' to not display
%   transp_flag: Required for using bicg to indicate whether operator is transposed or not
%
% Outputs:
%   HAOutput - Vector resulting from the action of the Hessian on a vector
%
% Hwan Goh 27/6/2018, University of Auckland, New Zealand (In Malaysia after one month trip but going home soon!)
%
% Notes: While there's alot of transforming back and forth 

if strcmp(transp_flag,'transp') || strcmp(transp_flag,'notransp') %For using bicg, this operator is symmetric
    
disp(' ')
disp('------------------------------------------------------')
disp('Computing Action of Hessian on a Vector ')
disp('------------------------------------------------------')
printf(['For the Case ' RunOptions.SaveFileName]);
if isnumeric(CurrentIteration)
    printf(['Current iteration number: ' num2str(CurrentIteration)]);
end
if isnumeric(hRelError)
    printf(['Computing line search, current relative error: ' num2str(hRelError)]);
end

hhat = halphahat(1:end-1);
alphahat = halphahat(end);

%=== If Positivity Constraint is Implemented ===%
if RunOptions.EWE_LS_LogExpTransform == 1
    preLogExpThRecon = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*hRecon)-1);
    pDerivhRecon = exp(Prior.LogExpTkappa*preLogExpThRecon)./(exp(Prior.LogExpTkappa*preLogExpThRecon)+1);
    hhat = pDerivhRecon.*hhat; %Note that h is actually varrho when considering positivity transform
end
if RunOptions.EWE_LS_SigmoidTransform == 1 || RunOptions.EWE_LS_ExponentialTransform == 1
    error('EWE_DGM2D_AdjStateHessianAction has not been set up for these positivity constraints')
end
[T11AdjTime0hat,T22AdjTime0hat,~,~,~,~,~] = EWE_DGM2D_HessianMisfit(RunOptions,hhat,FrwdFunction,SensorsI,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,MeshIN_Elm,MeshINodes,pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,PLOT);
                                                     
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second Order Total Derivative %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Interpolating To FEM Mesh ===%
T11AdjTime0hat = PrecomputedIntrplteObjectsI.InterpMatrix'*T11AdjTime0hat(:);
T22AdjTime0hat = PrecomputedIntrplteObjectsI.InterpMatrix'*T22AdjTime0hat(:);

%=== Total Derivative ===%
if RunOptions.EWE_LS_LogExpTransform == 1 
    preLogExpT1stDeriv = exp(Prior.LogExpTkappa*preLogExpThRecon)./(exp(Prior.LogExpTkappa*preLogExpThRecon)+1); %Derivative of positivity transform at prePosThRecon
    Direction_Prior = Prior.L_pr'*Prior.L_pr*hhat;
    Direction_Data = preLogExpT1stDeriv.*(T11AdjTime0hat + T22AdjTime0hat);
    Direction_alphapart = 2*alphaRecon*alphahat*Prior.L_pr'*Prior.L_pr*(Prior.Exp_h - preLogExpThRecon);
    Dirction_alpha = -2*alphaRecon*dot(Prior.L_pr'*Prior.L_pr*(Prior.Exp_h - preLogExpThRecon),hhat) + alphahat*norm(Prior.L_pr*(Prior.Exp_h - preLogExpThRecon),2)^2 + 2*Prior.Hyperprior_Std_alpha^(-2)*alphahat - MeshIN_Nodes*(1/alphaRecon^2)*alphahat;
    if RunOptions.EWE_DGM2D_UseGaussNewton == 1
        Direction_hpart = alphaRecon^2*Direction_Prior - Direction_Data - Direction_alphapart;
    end
    if RunOptions.EWE_DGM2D_UseNewton == 1
        preLogExpT2ndDeriv = Prior.LogExpTkappa*exp(Prior.LogExpTkappa*preLogExpThRecon)./((exp(Prior.LogExpTkappa*preLogExpThRecon)+1).^2); %Second Derivative of positivity transform at prePosThRecon
        Direction_Data2 = preLogExpT2ndDeriv.*(T11AdjTime0 + T22AdjTime0);
        Direction_hpart = alphaRecon^2*Direction_Prior - Direction_Data - Direction_Data2 - Direction_alphapart;
    end
    HAOutput = [Direction_hpart;Dirction_alpha];
end
end