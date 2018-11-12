function [T11AdjTime0,T22AdjTime0,T11AdjTimeStep,T22AdjTimeStep,T12AdjTimeStep,vxAdjTimeStep,vyAdjTimeStep] = EWE_DGM2D_HessianMisfit(RunOptions,h,FrwdFunction,SensorsI,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,MeshIN_Elm,MeshINodes,pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,PLOT)

% EWE_DGM2D_HessianMisfit computes the action of the Hessian misfit
%
% Inputs:
%   RunOptions:
%      NumberofTimeSteps - Number of time steps to be computed
%   p0 - Initial Pressure
%   SensorsI: 1 by number of sensors array containing 3 by 1 cells 
%      id - Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%      xy - coordinates of the sensor
%      l_iatsensor - 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%   L_ET11 - Error Model for T11
%   L_ET22 - Error Model for T22
%   L_ET12 - Error Model for T12
%   L_Evx - Error Model for vx
%   L_Evy - Error Model for vy
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   pinfoI - Information regarding p-refinement for p-nonconforming meshes on DGM inversion mesh
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI - Number of nodes per element on DGM inversion mesh
%   KI - Number of elements on DGM inversion mesh
%   dt - Size of time step, need negative of this to do reverse time stepping
%   PLOT: To plot or not to plot, that is the question
%
% Outputs:
%   T11AdjTime0 - Adjoint variable of T11 at initial time
%   T22AdjTime0 - Adjoint variable of T22 at initial time
%   T12AdjTime0 - Adjoint variable of T12 at initial time
%   vxAdjTime0 - Adjoint variable of vx at initial time
%   vxAdjTime0 - Adjoint variable of vx at initial time
%   T11AdjTimeStep - Adjoint variable of T11
%   T22AdjTimeStep - Adjoint variable of T22
%   T12AdjTimeStep - Adjoint variable of T12
%   vxAdjTimeStep - Adjoint variable of vx
%   vyAdjTimeStep - Adjoint variable of vy
%
% Hwan Goh 3/7/2018, University of Auckland, New Zealand (Back in NZ!)
%
% Note that this is just the incremental adjoint equation! Still need to
% multiply by mass matrix and take minuses (data part of the L_p'' term).

PLOT.DGMForward = 0; %Suppress plotting of wave propagation
PLOT.AdjWaveField = 0; %Suppress plotting of adjoint wave field propagation
PLOT.DGMForwardSensorData = 0; %Suppress plotting of sensory data
PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshI; %This is for plotting forward elastic wave propagation on inversion mesh. FwrdFunction is called which is EWE_DGM2D_LSExpRK4. However, that function calls PLOT.TRI_DGMMeshD when plotting, so here we replace it with DGMMeshI but re-use the same name
NumberofSensors = size(SensorsI,2);
NumberofTimeSteps = RunOptions.NumberofTimeSteps;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forward Wave Propagation %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hDGM = IntrplteOver2DTriangulatedMesh(MeshIN_Elm,MeshINodes,h,xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
IniCond = reshape(hDGM,NpI,KI)./(rhoI*RunOptions.SpecificHeatCoeff);
[T11p0TimeSteps,T22p0TimeSteps,T12p0TimeSteps,vxp0TimeSteps,vyp0TimeSteps] = FrwdFunction(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating Current f(p0) to Sensor Locations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T11InterpolatedTimeSteps = zeros(NumberofSensors,NumberofTimeSteps);
T22InterpolatedTimeSteps = zeros(NumberofSensors,NumberofTimeSteps);
T12InterpolatedTimeSteps = zeros(NumberofSensors,NumberofTimeSteps);
vxInterpolatedTimeSteps = zeros(NumberofSensors,NumberofTimeSteps);
vyInterpolatedTimeSteps = zeros(NumberofSensors,NumberofTimeSteps);
for t=1:NumberofTimeSteps
    %=== Full q Vector Data ===%
    if RunOptions.FullqVectorData == 1;
        for s=1:NumberofSensors
            T11InterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T11p0TimeSteps(SensorsI{s}.id,t);
            T22InterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T22p0TimeSteps(SensorsI{s}.id,t);
            T12InterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T12p0TimeSteps(SensorsI{s}.id,t);
            vxInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxp0TimeSteps(SensorsI{s}.id,t);
            vyInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyp0TimeSteps(SensorsI{s}.id,t);
        end
    end
    %=== Orthogonal Velocities Data ===%
    if RunOptions.VelocitiesData == 1;
        for s=1:NumberofSensors
            vxInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxp0TimeSteps(SensorsI{s}.id,t);
            vyInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyp0TimeSteps(SensorsI{s}.id,t);
        end
    end
end

%=== Computing L_E(q(p0)) ===%
L_ET11InterpolatedTimeSteps = L_ET11*T11InterpolatedTimeSteps(:);
L_ET22InterpolatedTimeSteps = L_ET22*T22InterpolatedTimeSteps(:);
L_ET12InterpolatedTimeSteps = L_ET12*T12InterpolatedTimeSteps(:);
L_EvxInterpolatedTimeSteps = L_Evx*vxInterpolatedTimeSteps(:);
L_EvyInterpolatedTimeSteps = L_Evy*vyInterpolatedTimeSteps(:);

L_ET11InterpolatedTimeSteps = reshape(L_ET11InterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);
L_ET22InterpolatedTimeSteps = reshape(L_ET22InterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);
L_ET12InterpolatedTimeSteps = reshape(L_ET12InterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);
L_EvxInterpolatedTimeSteps = reshape(L_EvxInterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);
L_EvyInterpolatedTimeSteps = reshape(L_EvyInterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);

%=== Computing L_E'L_E(q(p0)) ===%
L_EAdjL_ET11InterpolatedTimeSteps = L_ET11'*L_ET11InterpolatedTimeSteps(:);
L_EAdjL_ET22InterpolatedTimeSteps = L_ET22'*L_ET22InterpolatedTimeSteps(:);
L_EAdjL_ET12InterpolatedTimeSteps = L_ET12'*L_ET12InterpolatedTimeSteps(:);
L_EAdjL_EvxInterpolatedTimeSteps = L_Evx'*L_EvxInterpolatedTimeSteps(:);
L_EAdjL_EvyInterpolatedTimeSteps = L_Evy'*L_EvyInterpolatedTimeSteps(:);

L_EAdjL_ET11InterpolatedTimeSteps = reshape(L_EAdjL_ET11InterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);
L_EAdjL_ET22InterpolatedTimeSteps = reshape(L_EAdjL_ET22InterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);
L_EAdjL_ET12InterpolatedTimeSteps = reshape(L_EAdjL_ET12InterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);
L_EAdjL_EvxInterpolatedTimeSteps = reshape(L_EAdjL_EvxInterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);
L_EAdjL_EvyInterpolatedTimeSteps = reshape(L_EAdjL_EvyInterpolatedTimeSteps,NumberofSensors,NumberofTimeSteps);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Adjoint Wave Propagation %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Constructing Term Representing Second Derivative of Objective Functional ===%
T11SecondDerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
T22SecondDerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
T12SecondDerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
vxSecondDerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
vySecondDerivObjFunctTimeSteps = zeros(NpI*KI,NumberofTimeSteps);
for t=0:NumberofTimeSteps-2 %Must multiply by [l_1(s),l_2(s),...,l_Np(s)]' due to integration with the delta function and then considering for all i \in {1,...,Np}
    %=== Full q Vector Data ===% 
    if RunOptions.FullqVectorData == 1;
        for s=1:size(SensorsI,2)
            T11SecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = T11SecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_ET11InterpolatedTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
            T22SecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = T22SecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_ET22InterpolatedTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
            T12SecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = T12SecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_ET12InterpolatedTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
        end
    end
    for s=1:size(SensorsI,2)
        vxSecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = vxSecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_EvxInterpolatedTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
        vySecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) = vySecondDerivObjFunctTimeSteps(SensorsI{s}.id,NumberofTimeSteps - t) + L_EAdjL_EvyInterpolatedTimeSteps(s,NumberofTimeSteps - t)*SensorsI{s}.l_iatsensor';
    end
end
%=== Adjoint Wave Propagation ===%
[T11AdjTime0,T22AdjTime0,T11AdjTimeStep,T22AdjTimeStep,T12AdjTimeStep,vxAdjTimeStep,vyAdjTimeStep] = EWE_DGM2D_LSExpRK4Reverse(RunOptions,...
                                                                                                                              -T11SecondDerivObjFunctTimeSteps,-T22SecondDerivObjFunctTimeSteps,-T12SecondDerivObjFunctTimeSteps,-vxSecondDerivObjFunctTimeSteps,-vySecondDerivObjFunctTimeSteps,...
                                                                                                                               pinfoI,xI,yI,NpI,KI,dt,PLOT);     