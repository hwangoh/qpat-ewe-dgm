function [Gradient,Gradient_alpha,T11AdjTime0,T22AdjTime0] = EWE_DGM2D_AdjointStateMethod_Gradient_h(RunOptions,SensorsI,alphaRecon,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps,hRecon,Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,pinfoI,xI,yI,NpI,KI,InterpMatrix,dt,PLOT,ii)

% EWE_DGM2D_AdjointStateMethod_Gradient_h.m computes the derivative of the objective functional
% when constrained by the elastic wave equation discretized by the
% discontinuous Galerkin method
%
% Inputs:
%   RunOptions:
%      NumberofTimeSteps - Number of time steps to be computed
%      BC - Boundary conditions
%      EWE_SD_ExponentialTransform - Use exponential transform
%   SensorsI: 1 by number of sensors array containing 3 by 1 cells 
%      id - Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%      xy - coordinates of the sensor
%      l_iatsensor - 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%   alphaRecon - reconstruction of hyperprior parameter
%   L_ET11 - Error Model for T11
%   L_ET22 - Error Model for T22
%   L_ET12 - Error Model for T12
%   L_Evx - Error Model for vx
%   L_Evy - Error Model for vy
%   L_ET11ErrorTimeSteps - L_E(T11_d - T11(hRecon) - T11ErrorMean)
%   L_ET22ErrorTimeSteps - L_E(T22_d - T22(hRecon) - T22ErrorMean)
%   L_ET12ErrorTimeSteps - L_E(T12_d - T12(hRecon) - T12ErrorMean)
%   L_EvxErrorTimeSteps - L_E(vx_d - vx(hRecon) - vxErrorMean)
%   L_EvyErrorTimeSteps - L_E(vy_d - vy(hRecon) - vyErrorMean)
%   hRecon - Initial pressure at current reconstruction
%   Prior:
%      L_pr - Regularization operator
%      Exp_h - Expected value
%   MeshIElm: Number of elements by 3 array containing the indices of the nodes in each element
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshINodes: 2 by Number of nodes array for the coordinates of the nodes on the inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   pinfoI - Information regarding p-refinement for p-nonconforming meshes on DGM inversion mesh
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI - Number of nodes per element on DGM inversion mesh
%   KI - Number of elements on DGM inversion mesh%   dt - Size of time step, need negative of this to do reverse time stepping
%   PLOT: To plot or not to plot, that is the question
%
% Outputs:
%   Gradient - Steepest descent seach direction for p_0
%   Gradient_alpha - Steepest descent search direction for alpha
%   TiiAdjTime0 - Adjoint variable, to be used for computing Newton direction. 
%
% Hwan Goh, University of Auckland, New Zealand 13/09/2017

NumberofTimeSteps = RunOptions.NumberofTimeSteps;
NumberofSensors = size(SensorsI,2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Derivative of Objective Functional %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Computing L_E'L_E(q_d - q(hRecon)) ===%
L_EAdjL_ET11ErrorTimeSteps = L_ET11'*L_ET11ErrorTimeSteps(:);
L_EAdjL_ET22ErrorTimeSteps = L_ET22'*L_ET22ErrorTimeSteps(:);
L_EAdjL_ET12ErrorTimeSteps = L_ET12'*L_ET12ErrorTimeSteps(:);
L_EAdjL_EvxErrorTimeSteps = L_Evx'*L_EvxErrorTimeSteps(:);
L_EAdjL_EvyErrorTimeSteps = L_Evy'*L_EvyErrorTimeSteps(:);

L_EAdjL_ET11ErrorTimeSteps = reshape(L_EAdjL_ET11ErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);
L_EAdjL_ET22ErrorTimeSteps = reshape(L_EAdjL_ET22ErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);
L_EAdjL_ET12ErrorTimeSteps = reshape(L_EAdjL_ET12ErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);
L_EAdjL_EvxErrorTimeSteps = reshape(L_EAdjL_EvxErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);
L_EAdjL_EvyErrorTimeSteps = reshape(L_EAdjL_EvyErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);

%=== Constructing Term Representing Derivative of Objective Functional ===%
T11DerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
T22DerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
T12DerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
vxDerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
vyDerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
for t=0:NumberofTimeSteps-2 %Must multiply by [l_1(s),l_2(s),...,l_Np(s)]' due to integration with the delta function and then considering for all i \in {1,...,Np}
    %=== Full q Vector Data ===% 
    if RunOptions.FullqVectorData == 1;
        for s=1:size(SensorsI,2)
            T11DerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = T11DerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_ET11ErrorTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
            T22DerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = T22DerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_ET22ErrorTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
            T12DerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = T12DerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_ET12ErrorTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
        end
    end
    for s=1:size(SensorsI,2)
        vxDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = vxDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_EvxErrorTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
        vyDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = vyDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_EvyErrorTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
    end
end

%% %%%%%%%%%%%%%%%%%%%%%
%%% Adjoint Equation %%%
%%%%%%%%%%%%%%%%%%%%%%%%
[T11AdjTime0,T22AdjTime0,~,~,~,~,~,~,~,~] = EWE_DGM2D_LSExpRK4Reverse(RunOptions,...
                                                            T11DerivObjFunctTimeSteps,T22DerivObjFunctTimeSteps,T12DerivObjFunctTimeSteps,vxDerivObjFunctTimeSteps,vyDerivObjFunctTimeSteps,...
                                                            pinfoI,xI,yI,NpI,KI,dt,PLOT);                                                                                                                                                  
%% %%%%%%%%%%%%%%%%%%%%%
%%% Total Derivative %%%
%%%%%%%%%%%%%%%%%%%%%%%%
%=== Interpolating To FEM Mesh ===%
T11AdjTime0 = InterpMatrix'*T11AdjTime0(:);
T22AdjTime0 = InterpMatrix'*T22AdjTime0(:);

%=== Total Derivative ===%
if RunOptions.EWE_LS_LogExpTransform == 1
    preLogExpThRecon = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*hRecon)-1); %Inverse mapping of positivity transform     
    preLogExpT1stDeriv = exp(Prior.LogExpTkappa*preLogExpThRecon)./(exp(Prior.LogExpTkappa*preLogExpThRecon)+1); %Derivative of positivity transform at prePosThRecon
    Gradient_Data = preLogExpT1stDeriv.*(T11AdjTime0 + T22AdjTime0);
    Gradient_Prior = Prior.L_pr'*Prior.L_pr*(Prior.Exp_h - preLogExpThRecon);
    Gradient_alpha = alphaRecon*norm(Prior.L_pr*(Prior.Exp_h - preLogExpThRecon),2)^2 - 2*Prior.Hyperprior_Std_alpha^(-2)*(Prior.Hyperprior_Exp_alpha - alphaRecon) - MeshIN_Nodes*(1/alphaRecon);
end
if RunOptions.EWE_LS_SigmoidTransform == 1
    preSigThRecon = (atanh((hRecon - Prior.SigTalpha)/Prior.SigTbeta))/Prior.SigTkappa + Prior.SigTdelta; %Inverse mapping of positivity transform
    preSigT1stDeriv = Prior.SigTbeta*((sech(Prior.SigTkappa*(preSigThRecon - Prior.SigTdelta))).^2)*Prior.SigTkappa; %Derivative of positivity transform at prePosThRecon
    Gradient_Data = preSigT1stDeriv.*(T11AdjTime0 + T22AdjTime0);
    Gradient_Prior = Prior.L_pr'*Prior.L_pr*(Prior.Exp_h - preSigThRecon);
    Gradient_alpha = alphaRecon*norm(Prior.L_pr*(Prior.Exp_h - preSigThRecon),2)^2 - 2*Prior.Hyperprior_Std_alpha^(-2)*(Prior.Hyperprior_Exp_alpha - alphaRecon) - MeshIN_Nodes*(1/alphaRecon);
end
if RunOptions.EWE_LS_ExponentialTransform == 1
    preExpThRecon = log(hRecon); %Inverse mapping of positivity transform
    preExpT1stDeriv = exp(preExpThRecon); %Derivative of positivity transform at prePosThRecon
    Gradient_Data = preExpT1stDeriv.*(T11AdjTime0 + T22AdjTime0);
    Gradient_Prior = Prior.L_pr'*Prior.L_pr*(Prior.Exp_h - preExpThRecon);
    Gradient_alpha = alphaRecon*norm(Prior.L_pr*(Prior.Exp_h - preExpThRecon),2)^2 - 2*Prior.Hyperprior_Std_alpha^(-2)*(Prior.Hyperprior_Exp_alpha - alphaRecon) - MeshIN_Nodes*(1/alphaRecon);
end
if RunOptions.EWE_LS_LogExpTransform == 0 && RunOptions.EWE_LS_SigmoidTransform == 0 && RunOptions.EWE_LS_ExponentialTransform == 0
    Gradient_Data = T11AdjTime0 + T22AdjTime0;
    Gradient_Prior = Prior.L_pr'*Prior.L_pr*(Prior.Exp_h - hRecon);
    Gradient_alpha = alphaRecon*norm(Prior.L_pr*(Prior.Exp_h - hRecon),2)^2 - 2*Prior.Hyperprior_Std_alpha^(-2)*(Prior.Hyperprior_Exp_alpha - alphaRecon) - MeshIN_Nodes*(1/alphaRecon);
end

%=== Hyperprior Gradient ===%
if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 0;
    alphaRecon = 1;
    Gradient_alpha = 0;
end
Gradient = -Gradient_Data - alphaRecon^2*Gradient_Prior;


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%==== Debugging Gradient ===%
% DataFigure = figure;
% drawnow
% trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),Gradient_Data);
% shading interp %thanks Ru!
% colormap(jet(256))
% 
% PriorFigure = figure;
% drawnow
% trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),alphaRecon^2*Gradient_Prior)
% shading interp %thanks Ru!
% colormap(jet(256))
% 
% GradientFigure = figure;
% drawnow
% trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2), -Gradient)
% shading interp %thanks Ru!
% colormap(jet(256))
% 
% keyboard
% 
% close(DataFigure) 
% close(PriorFigure) 
% close(GradientFigure)
