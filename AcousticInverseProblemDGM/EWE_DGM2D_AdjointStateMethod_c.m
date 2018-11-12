function [gradient] = EWE_DGM2D_AdjointStateMethod_c(RunOptions,SensorsI,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps,cpRecon,csRecon,Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,pinfoI,xI,yI,NpI,KI,VertexNodesGlobalIndicesFEMI,dt,PLOT)

% EWE_DGM2D_AdjointStateMethod_c.m computes the derivative of the objective functional
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
%   L_ET11 - Error Model for T11
%   L_ET22 - Error Model for T22
%   L_ET12 - Error Model for T12
%   L_Evx - Error Model for vx
%   L_Evy - Error Model for vy
%   L_ET11ErrorTimeSteps - L_E(T11_d - T11(p0Recon) - T11ErrorMean)
%   L_ET22ErrorTimeSteps - L_E(T22_d - T22(p0Recon) - T22ErrorMean)
%   L_ET12ErrorTimeSteps - L_E(T12_d - T12(p0Recon) - T12ErrorMean)
%   L_EvxErrorTimeSteps - L_E(vx_d - vx(p0Recon) - vxErrorMean)
%   L_EvyErrorTimeSteps - L_E(vy_d - vy(p0Recon) - vyErrorMean)
%   Prior:
%      L_pr - Regularization operator
%      Exp_p0 - Expected value
%   MeshIElm: Number of elements by 3 array containing the indices of the nodes in each element
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshINodes: 2 by Number of nodes array for the coordinates of the nodes on the inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   pinfoI - Information regarding p-refinement for p-nonconforming meshes on DGM inversion mesh
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI - Number of nodes per element on DGM inversion mesh
%   KI - Number of elements on DGM inversion mesh
%   VertexNodesGlobalIndicesFEMI: indices referring to the vertices of elements. This is for circumventing the gradient inconsistencies arising from the adjoint state method; only the gradient values at these points are considered
%   dt - Size of time step, need negative of this to do reverse time stepping
%   PLOT: To plot or not to plot, that is the question
%
% Outputs:
%   gradient - Steepest descent seach direction
%
% Hwan Goh, University of Auckland, New Zealand 21/04/2018

NumberofTimeSteps = RunOptions.NumberofTimeSteps;
NumberofSensors = size(SensorsI,2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Derivative of Objective Functional %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Computing L_E'L_E(q_d - q(p0Recon)) ===%
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
[~,~,~,~,~,T11AdjTimeStep,T22AdjTimeStep,T12AdjTimeStep,vxAdjTimeStep,vyAdjTimeStep] = EWE_DGM2D_LSExpRK4Reverse(RunOptions,...
                                                                                                                 T11DerivObjFunctTimeSteps,T22DerivObjFunctTimeSteps,T12DerivObjFunctTimeSteps,vxDerivObjFunctTimeSteps,vyDerivObjFunctTimeSteps,...
                                                                                                                 pinfoI,xI,yI,NpI,KI,dt,PLOT);                                                                                                                                                  
%% %%%%%%%%%%%%%%%%%%%%%
%%% Total Derivative %%%
%%%%%%%%%%%%%%%%%%%%%%%%
%=== Interpolating To FEM Mesh ===%
T11AdjTime0 = T11AdjTime0(VertexNodesGlobalIndicesFEMI);
T22AdjTime0 = T22AdjTime0(VertexNodesGlobalIndicesFEMI);

%=== Total Derivative ===%
if RunOptions.EWE_SD_ExponentialTransform == 1
    MassMatrixFEM = ConstructExpTMassMatrixFEM2D(MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,p0DGM);
end
gradient = -MassMatrixFEM*(T11AdjTime0 + T22AdjTime0 + Prior.L_pr'*Prior.L_pr*(Prior.Exp_p0 - p0DGM));
gradient = gradient/(norm(gradient,2));

%=== Total Derivative ===%
V = pinfoI(2).V;
MassMatrixDGM = inv(V*V');

RegTerm_cp = -MassMatrixDGM*(Prior.L_cp'*Prior.L_cp*(Prior.Exp_cp - cpRecon));
RegTerm_cs = -MassMatrixDGM*(Prior.L_cs'*Prior.L_cs*(Prior.Exp_cs - csRecon));

