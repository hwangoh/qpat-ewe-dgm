function [T11InterpolatedTimeSteps,T22InterpolatedTimeSteps,T12InterpolatedTimeSteps,vxInterpolatedTimeSteps,vyInterpolatedTimeSteps,Error_qh,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps] = EWE_DGM2D_ObjectiveFunctionalErrorTerm(RunOptions,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,T11NewTimeSteps,T22NewTimeSteps,T12NewTimeSteps,vxNewTimeSteps,vyNewTimeSteps,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,SensorsI)

% EWE_DGM2D_ObjectiveFunctionalErrorTerm takes data and computed time steps from an initial condition to form the error functional L_E(q_d - q(hRecon) - ErrorMean)
%
% Inputs:
%   RunOptions
%   L_ET11 - Error Model for T11
%   L_ET22 - Error Model for T22
%   L_ET12 - Error Model for T12
%   L_Evx - Error Model for vx
%   L_Evy - Error Model for vy
%   T11DataTimeSteps - Time step data for T11, dimensions are NumberofSensors by NumberofTimeSteps
%   T22DataTimeSteps - Time step data for T22, dimensions are NumberofSensors by NumberofTimeSteps
%   T12DataTimeSteps - Time step data for T12, dimensions are NumberofSensors by NumberofTimeSteps
%   vxDataTimeSteps -  Time step data for vx, dimensions are NumberofSensors by NumberofTimeSteps
%   vyDataTimeSteps -  Time step data for vy, dimensions are NumberofSensors by NumberofTimeSteps
%   T11NewTimeSteps - New time steps for T11, dimensions are NumberofSensors by NumberofTimeSteps
%   T22NewTimeSteps - New time steps for T22, dimensions are NumberofSensors by NumberofTimeSteps
%   T12NewTimeSteps - New time steps for T12, dimensions are NumberofSensors by NumberofTimeSteps
%   vxNewTimeSteps -  New time steps for vx, dimensions are NumberofSensors by NumberofTimeSteps
%   vyNewTimeSteps -  New time steps for vy, dimensions are NumberofSensors by NumberofTimeSteps
%   T11ErrorMean - Mean of error model for T11, dimensions are NumberofSensors by NumberofTimeSteps
%   T22ErrorMean - Mean of error model for T22, dimensions are NumberofSensors by NumberofTimeSteps
%   T12ErrorMean - Mean of error model for T12, dimensions are NumberofSensors by NumberofTimeSteps
%   vxErrorMean - Mean of error model for vx, dimensions are NumberofSensors by NumberofTimeSteps
%   vyErrorMean - Mean of error model for vy, dimensions are NumberofSensors by NumberofTimeSteps
%   SensorsI: 1 by number of sensors array containing 3 by 1 cells 
%      id - Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%      xy - coordinates of the sensor
%      l_iatsensor - 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%
% Outputs: Total Error of q ||q_d - q(hRecon)||
%   T11InterpolatedTimeSteps - T11NewTimeSteps interpolated to sensors
%   T22InterpolatedTimeSteps - T22NewTimeSteps interpolated to sensors
%   T12InterpolatedTimeSteps - T12NewTimeSteps interpolated to sensors
%   vxInterpolatedTimeSteps - vxNewTimeSteps interpolated to sensors
%   vyInterpolatedTimeSteps - vyNewTimeSteps interpolated to sensors
%   Error_qh - ||q_d - q(hRecon)||
%   L_ET11ErrorTimeSteps - L_E(T11_d - T11(p^m) - T11ErrorMean)
%   L_ET22ErrorTimeSteps - L_E(T22_d - T22(p^m) - T22ErrorMean)
%   L_ET12ErrorTimeSteps - L_E(T12_d - T12(p^m) - T12ErrorMean)
%   L_EvxErrorTimeSteps - L_E(vx_d - vx(p^m) - vxErrorMean)
%   L_EvyErrorTimeSteps - L_E(vy_d - vy(p^m) - vyErrorMean) 
%
% Hwan Goh, University of Auckland, New Zealand 17/02/2018

NumberofSensors = size(SensorsI,2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating Current f(hRecon) to Sensor Locations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T11InterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
T22InterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
T12InterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
vxInterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
vyInterpolatedTimeSteps = zeros(NumberofSensors,RunOptions.NumberofTimeSteps);
for t=1:RunOptions.NumberofTimeSteps
    %=== Full q Vector Data ===%
    if RunOptions.FullqVectorData == 1;
        for s=1:NumberofSensors
            T11InterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T11NewTimeSteps(SensorsI{s}.id,t);
            T22InterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T22NewTimeSteps(SensorsI{s}.id,t);
            T12InterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T12NewTimeSteps(SensorsI{s}.id,t);
            vxInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxNewTimeSteps(SensorsI{s}.id,t);
            vyInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyNewTimeSteps(SensorsI{s}.id,t);
        end
    end
    %=== Orthogonal Velocities Data ===%
    if RunOptions.VelocitiesData == 1;
        for s=1:NumberofSensors
            vxInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxNewTimeSteps(SensorsI{s}.id,t);
            vyInterpolatedTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyNewTimeSteps(SensorsI{s}.id,t);
        end
    end
end

%% %%%%%%%%%%
%%% Error %%%
%%%%%%%%%%%%%
%=== Computing q_d - q(hRecon) - ErrorMean ===%
T11ErrorTimeSteps = T11DataTimeSteps - T11InterpolatedTimeSteps - T11ErrorMean;
T22ErrorTimeSteps = T22DataTimeSteps - T22InterpolatedTimeSteps - T22ErrorMean;
T12ErrorTimeSteps = T12DataTimeSteps - T12InterpolatedTimeSteps - T12ErrorMean;
vxErrorTimeSteps = vxDataTimeSteps - vxInterpolatedTimeSteps - vxErrorMean;
vyErrorTimeSteps = vyDataTimeSteps - vyInterpolatedTimeSteps - vyErrorMean;

if RunOptions.EWE_LS_RemoveTimesSteps ~= 0;
    T11ErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
    T22ErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
    T12ErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
    vxErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
    vyErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
end

%=== Computing L_E(q_d - q(hRecon) - ErrorMean) ===% %This precomputed because it is used later in EWE_DGM2D_AdjointStateMethod and EWE_DGM2D_PNonConf_StrainAdjoint
L_ET11ErrorTimeSteps = L_ET11*T11ErrorTimeSteps(:);
L_ET22ErrorTimeSteps = L_ET22*T22ErrorTimeSteps(:);
L_ET12ErrorTimeSteps = L_ET12*T12ErrorTimeSteps(:);
L_EvxErrorTimeSteps = L_Evx*vxErrorTimeSteps(:);
L_EvyErrorTimeSteps = L_Evy*vyErrorTimeSteps(:);

L_ET11ErrorTimeSteps = reshape(L_ET11ErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);
L_ET22ErrorTimeSteps = reshape(L_ET22ErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);
L_ET12ErrorTimeSteps = reshape(L_ET12ErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);
L_EvxErrorTimeSteps = reshape(L_EvxErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);
L_EvyErrorTimeSteps = reshape(L_EvyErrorTimeSteps,NumberofSensors,RunOptions.NumberofTimeSteps);

%=== Total Error of q ||q_d - q(hRecon)|| ===%
Error_qh = norm([T11ErrorTimeSteps(:);T22ErrorTimeSteps(:);T12ErrorTimeSteps(:);vxErrorTimeSteps(:);vyErrorTimeSteps(:)],2);
