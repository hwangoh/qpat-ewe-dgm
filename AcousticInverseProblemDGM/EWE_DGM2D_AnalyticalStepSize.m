function StepSize = EWE_DGM2D_AnalyticalStepSize(RunOptions,FrwdFunction,Dirctn,hRecon,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps,Prior,N_ElmI,NodesI,PrecomputedIntrplteObjectsI,SensorsI,pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT)

% EWE_DGM2D_AnalyticalStepSize uses the linearity of the dependence of the state
% variable q on the control parameter h to compute the quadratic
% corresponding to the regularized error functional.
%
% Inputs:
%   RunOptions:
%   FrwdFunction - Script for the forward function of h
%   Dirctn - Current Step Direction
%   hRecon - Current estimate of h
%   L_ET11 - Error Model for T11
%   L_ET22 - Error Model for T22
%   L_ET12 - Error Model for T12
%   L_Evx - Error Model for vx
%   L_Evy - Error Model for vy
%   L_ET11ErrorTimeSteps - L_E(T11_d - T11(p^m) - T11ErrorMean)
%   L_ET22ErrorTimeSteps - L_E(T22_d - T22(p^m) - T22ErrorMean)
%   L_ET12ErrorTimeSteps - L_E(T12_d - T12(p^m) - T12ErrorMean)
%   L_EvxErrorTimeSteps - L_E(vx_d - vx(p^m) - vxErrorMean)
%   L_EvyErrorTimeSteps - L_E(vy_d - vy(p^m) - vyErrorMean)
%   Prior:
%      L_pr - Regularization operator
%      Exp_h - Expected value
%   N_ElmI - Number of elements in optical mesh
%   NodesI - Nx2 matrix containing coordinates of nodes on inverse optical mesh
%   PrecompIntrplteObjectsI - Objects that depend on the inverse Mesh nodes. May have been computed earlier and so can 
%                             be called here to avoid repeating computations. Set to 0 in function call if you
%                             want to compute new objects for a new mesh.
%   SensorsI: 1 by number of sensors array containing 3 by 1 cells 
%      id - Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%      xy - coordinates of the sensor
%      l_iatsensor - 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%   pinfoI - Invesion mesh information regarding p-refinement for p-nonconforming meshes
%   xI: x coordinates of the DGM inverse mesh nodes, for interpolation onto optical inverse mesh
%   yI: y coordinates of the DGM inverse mesh nodes, for interpolation onto optical inverse mesh
%   NpI: Number of grid points in one element of the DGM inverse mesh
%   KI: Number of elements of the DGM inverse mesh
%   rhoI: Medium density of inverse mesh for computing initial condition from h
%   dt: Time step size
%   PLOT - To plot or not to plot, that is the question
%
% Outputs:
%   StepSize - the optimum step size
%
% Notes: Does not match analytical stepsize on document! L_11ETijErrorTimeSteps includes mean already!
%
% Hwan Goh, University of Auckland, New Zealand 12/10/2017

disp('Computing Analytical Step Size ')
NumberofSensors = size(SensorsI,2);
PLOT.DGMForward = 0;
PLOT.DGMPlotzAxis = [0 1];
PLOT.ColourAxis = [0 1];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forward Function With Step Direction As Initial Condition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DirctnDGM = IntrplteOver2DTriangulatedMesh(N_ElmI,NodesI,Dirctn,xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
IniCond = reshape(DirctnDGM,NpI,KI)./(rhoI*RunOptions.SpecificHeatCoeff);

[T11LSTimeSteps,T22LSTimeSteps,T12LSTimeSteps,vxLSTimeSteps,vyLSTimeSteps] = FrwdFunction(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT); 

%=== Interpolating Current fh to Sensor Locations ===%
T11LSInterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
T22LSInterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
T12LSInterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
vxLSInterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
vyLSInterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
for t=1:RunOptions.NumberofTimeSteps
    %=== Full q Vector Data ===%
    if RunOptions.FullqVectorData == 1;
        for s=1:NumberofSensors
            T11LSInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T11LSTimeSteps(SensorsI{s}.id,t);
            T22LSInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T22LSTimeSteps(SensorsI{s}.id,t);
            T12LSInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T12LSTimeSteps(SensorsI{s}.id,t);
            vxLSInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxLSTimeSteps(SensorsI{s}.id,t);
            vyLSInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyLSTimeSteps(SensorsI{s}.id,t);
        end
    end
    %=== Velocities Data ===%
    if RunOptions.VelocitiesData == 1;
        for s=1:NumberofSensors
            vxLSInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxLSTimeSteps(SensorsI{s}.id,t);
            vyLSInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyLSTimeSteps(SensorsI{s}.id,t);
        end
    end
end

%=== Applying Error Model ===%
L_ET11LSInterpolatedTimeSteps = L_ET11*T11LSInterpolatedTimeSteps(:);
L_ET22LSInterpolatedTimeSteps = L_ET22*T22LSInterpolatedTimeSteps(:);
L_ET12LSInterpolatedTimeSteps = L_ET12*T12LSInterpolatedTimeSteps(:);
L_EvxLSInterpolatedTimeSteps = L_Evx*vxLSInterpolatedTimeSteps(:);
L_EvyLSInterpolatedTimeSteps = L_Evy*vyLSInterpolatedTimeSteps(:);

L_EfLS = [L_ET11LSInterpolatedTimeSteps;L_ET22LSInterpolatedTimeSteps;L_ET12LSInterpolatedTimeSteps;L_EvxLSInterpolatedTimeSteps;L_EvyLSInterpolatedTimeSteps];
L_EError = [L_ET11ErrorTimeSteps(:);L_ET22ErrorTimeSteps(:);L_ET12ErrorTimeSteps(:);L_EvxErrorTimeSteps(:);L_EvyErrorTimeSteps(:)];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forming Quadratic and Solving for Step Size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
a = (1/2)*(norm(L_EfLS,2)^2 + norm(Prior.L_pr*Dirctn,2)^2);
b = -dot(L_EfLS,L_EError) - dot(Prior.L_pr*Dirctn,Prior.L_pr*(Prior.Exp_h - hRecon));
c = (1/2)*(norm(L_EError,2)^2 + norm(Prior.L_pr*(Prior.Exp_h - hRecon),2)^2);

s = -10:0.01:1e4;
ErrorFunctionalQuadratic = zeros(1,size(s,2));
for i=1:size(s,2)
    ErrorFunctionalQuadratic(i) = a*(s(i)^2) + b*s(i) + c;
end
% figure
% plot(s,ErrorFunctionalQuadratic);
[~,index] = min(ErrorFunctionalQuadratic);

StepSize = s(index);