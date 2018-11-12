function [SampleMeanT11,SampleMeanT22,SampleMeanT12,SampleMeanvx,SampleMeanvy,SampleCovT11,SampleCovT22,SampleCovT12,SampleCovvx,SampleCovvy]=EWE_DGM2D_AE_ConstructStats(RunOptions,MeshANodes,MeshAElements,MeshADimensns,MeshINodes,MeshIElements,xA,yA,NpA,KA,PrecomputedIntrplteObjectsA,pinfoA,xI,yI,NpI,KI,PrecomputedIntrplteObjectsI,pinfoI,SensorsA,SensorsI,EdgeSensorInd,dt,Prior,PLOT)

% EWE_DGM2D_AE_ConstructStats calculates the approximation error mean and approximation error 
% covariance for the enhanced error model. The priors samples are drawn
% from a prior model of 
%
% Inputs:
%   RunOptions:
%              NumberofTimeSteps
%              AEN_Samples - number of samples we wish to use to calculate approximation error
%   MeshANodes,MeshINodes - Number of FEM Nodes by 2 array storing the coordinates of the nodes
%   MeshAElements,MeshIElements - Number of Elements by 3 array where the ith row contains the indices of the nodes in the ith element
%   MeshADimensns - Width and height of the accurate mesh
%   xA, xI - x-coordinates of nodes on the DGM Mesh
%   yA, yI - y-coordinates of nodes on the DGM Mesh
%   NpA, NpI - Number of nodes per element on the DGM Mesh
%   KA, KI - Number of elements
%   PrecomputedIntrplteObjectsA,PrecomputedIntrplteObjectsI:
%                             Objects that depend on the inverse Mesh nodes. May have been computed earlier and so can 
%                             be called here to avoid repeating computations. Set to 0 in function call if you
%                             want to compute new objects for a new mesh.
%   pinfoA,pinfoI - Invesion mesh information regarding p-refinement for p-nonconforming meshes
%   SensorsA, SensorsI: 1 by number of sensors array containing 3 by 1 cells 
%                       - id: Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%                       - xy: coordinates of the sensor
%                       - l_iatsensor: 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%   EdgeSensorInd - Sensor indices along one edge, for plotting sensor data
%   dt - size of time steps
%   Prior:
%      corr_p0 - 2*1 vector = [correlation_x, correlation_y]. The larger the number,
%          the more points near to a marginalisation point are correlated
%          in the x and y direction.
%      bounds_p0 - approximately [p0_min, p0_max]
%      Exp_p0 - initial guess of p0
%   PLOT - To plot or not to plot, that is the question
%
% Outputs:
%      SampleMeanTij, SampleMeanvi - approximation error mean; dimensions are NumberofSensors by NumberofTimeSteps after using reshape
%      SampleCovTij, SampleCovvi - approximation error covariance; dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%
% Hwan Goh 15/01/2018, University of Auckland, New Zealand

PLOT.DGMForward = 0; %Suppress plotting of wave propagation
PLOT.PriorMargPoints = 0; %Suppress plotting of marginalisation points
PLOT.PriorSamples = 0; %Suppress plotting of prior samples
PLOT.DGMForwardSensorData = 0; %Suppress plotting of sensory data
N_Samples = RunOptions.AEN_Samples;
NumberofSensors = size(SensorsI,2);
NumberofTimeSteps = RunOptions.NumberofTimeSteps;
printf([num2str(N_Samples) ' number of samples to be computed for approximation error ']);

%=========================================================================%
%                         Generating Samples
%=========================================================================%
%% %%%%%%%%%%%%%%%%%%%
%%% Accurate Model %%%
%%%%%%%%%%%%%%%%%%%%%%
printf('Computing samples of accurate model');

%=== Initial Pressure Samples ===%
[~,~,~,~,p0AS_FEM]=SmoothnessPrior_Informative(MeshANodes,MeshAElements,2*MeshADimensns,Prior.InformSmooth_Bounds_p0,Prior.InformSmooth_Corr_p0,Prior.Exp_p0,N_Samples,PLOT);
% [~,~,p0AS_FEM] = SmoothnessPrior_AutoCorr(MeshANodes,Prior.Exp_p0,Prior.AC_Var_p0,Prior.AC_Corr_p0,N_Samples,PLOT);

%=== Forward Data Samples ===%
T11ASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
T22ASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
T12ASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
vxASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
vyASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);

PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshA; %This is for plotting forward elastic wave propagation on inversion mesh. FwrdFunction is called which is EWE_DGM2D_LSExpRK4. However, that function calls PLOT.TRI_DGMMeshD when plotting, so here we replace it with DGMMeshA but re-use the same name

%=== Accurate Forward Data ===%
for n=1:N_Samples
    printf(['\nComputing Accurate Model Sample ' num2str(n) ' of ' num2str(N_Samples)]);
    printf(['For the Case ' RunOptions.SaveFileName]);
    [p0AS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),p0AS_FEM(:,n),xA,yA,NpA*KA,PrecomputedIntrplteObjectsA);
    IniCond = p0AS_DGM/(RunOptions.ElasticMediumDensity*RunOptions.SpecificHeatCoeff);
    if RunOptions.TimeLSERK4 == 1;
        [T11ASTimeSteps,T22ASTimeSteps,T12ASTimeSteps,vxASTimeSteps,vyASTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,xA,yA,NpA,KA,pinfoA,dt,PLOT);
    end
    T11ASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    T22ASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    T12ASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    vxASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    vyASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    %=== Interpolation to Form Sensory Data ===%
    for t=1:NumberofTimeSteps
        for s=1:NumberofSensors
            if RunOptions.FullqVectorData == 1;
                T11ASDataTimeSteps(s,t) = SensorsA{s}.l_iatsensor*T11ASTimeSteps(SensorsA{s}.id,t);
                T22ASDataTimeSteps(s,t) = SensorsA{s}.l_iatsensor*T22ASTimeSteps(SensorsA{s}.id,t);
                T12ASDataTimeSteps(s,t) = SensorsA{s}.l_iatsensor*T12ASTimeSteps(SensorsA{s}.id,t);
            end
            vxASDataTimeSteps(s,t) = SensorsA{s}.l_iatsensor*vxASTimeSteps(SensorsA{s}.id,t);
            vyASDataTimeSteps(s,t) = SensorsA{s}.l_iatsensor*vyASTimeSteps(SensorsA{s}.id,t);
        end
    end
    %=== Plotting Sensory Data ===%
    if RunOptions.UseFullDomainData ~= 1 && PLOT.DGMForwardSensorData == 1;
        figure
        for ii=1:RunOptions.NumberofTimeSteps
            plot(1:1:length(EdgeSensorInd),sqrt(vxASDataTimeSteps(EdgeSensorInd,ii).^2 + vyASDataTimeSteps(EdgeSensorInd,ii).^2),'o')
            ylim([0,600])
            pause(0.1)
        end
    end
    T11ASamplesDataTimeSteps(:,n) = T11ASDataTimeSteps(:);
    T22ASamplesDataTimeSteps(:,n) = T22ASDataTimeSteps(:);
    T12ASamplesDataTimeSteps(:,n) = T12ASDataTimeSteps(:);
    vxASamplesDataTimeSteps(:,n) = vxASDataTimeSteps(:);
    vyASamplesDataTimeSteps(:,n) = vyASDataTimeSteps(:);
end
if RunOptions.FullqVectorData == 1;
clear T11ASTimeSteps T22ASTimeSteps T12ASTimeSteps T11ASDataTimeSteps T22ASDataTimeSteps T12ASDataTimeSteps
end
clear vxASTimeSteps vyASTimeSteps vxASDataTimeSteps vyASDataTimeSteps

%=== Saving Samples and Forward Data ===% Will only be using 00034 Mesh for accurate model
SaveFileNameAESamplesForward = sprintf('AESamplesForward-%s',RunOptions.SaveFileName)
save(SaveFileNameAESamplesForward,'p0AS_FEM','T11ASamplesDataTimeSteps','T22ASamplesDataTimeSteps','T12ASamplesDataTimeSteps','vxASamplesDataTimeSteps','vyASamplesDataTimeSteps','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%
%%% Less Accurate Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf('Computing samples of less accurate model ');
%=== Initial Pressure Samples ===%
[p0IS_FEM(:,1), PrecomputedIntrplteObjectsDtoIFEM] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),p0AS_FEM(:,1),MeshINodes(:,1),MeshINodes(:,2),size(MeshINodes,1),0);
for n=2:N_Samples %Just so we do not need to compute the interpolation objects over and over again
    [p0IS_FEM(:,n), ~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),p0AS_FEM(:,n),MeshINodes(:,1),MeshINodes(:,2),size(MeshINodes,1),PrecomputedIntrplteObjectsDtoIFEM);
end
if PLOT.PriorSamples==1;
    figure
    FemPlot2D(MeshINodes,MeshIElements,p0IS_FEM(:,1:min(N_Samples,10)))
    title(PLOT.Figure_Prior_Title,'FontWeight','bold')
end

%=== Forward Data Samples ===%
T11ISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
T22ISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
T12ISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
vxISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
vyISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);

PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshI; %This is for plotting forward elastic wave propagation on inversion mesh. FwrdFunction is called which is EWE_DGM2D_LSExpRK4. However, that function calls PLOT.TRI_DGMMeshD when plotting, so here we replace it with DGMMeshI but re-use the same name

%=== Less Accurate Forward Data ===%
NumberofSensors = size(SensorsI,2);
for n=1:N_Samples
    printf(['\nComputing Less Accurate Model Sample ' num2str(n) ' of ' num2str(N_Samples)]);
    printf(['For the Case ' RunOptions.SaveFileName]);
    [p0IS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshIElements,1),MeshINodes*(1/RunOptions.ScalingOptical),p0IS_FEM(:,n),xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
    IniCond = p0IS_DGM/(RunOptions.ElasticMediumDensity*RunOptions.SpecificHeatCoeff);
    if RunOptions.TimeLSERK4 == 1;
        [T11ISTimeSteps,T22ISTimeSteps,T12ISTimeSteps,vxISTimeSteps,vyISTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT);
    end
    T11ISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    T22ISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    T12ISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    vxISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    vyISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
    %=== Interpolation to Form Sensory Data ===%
    for t=1:NumberofTimeSteps     
        for s=1:NumberofSensors
            if RunOptions.FullqVectorData == 1;
                T11ISDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T11ISTimeSteps(SensorsI{s}.id,t);
                T22ISDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T22ISTimeSteps(SensorsI{s}.id,t);
                T12ISDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T12ISTimeSteps(SensorsI{s}.id,t);
            end
            vxISDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxISTimeSteps(SensorsI{s}.id,t);
            vyISDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyISTimeSteps(SensorsI{s}.id,t);
        end
    end
    %=== Plotting Sensor Data ===%
    if RunOptions.UseFullDomainData ~= 1 && PLOT.DGMForwardSensorData == 1;
        figure
        for ii=1:RunOptions.NumberofTimeSteps
            plot(1:1:length(EdgeSensorInd),sqrt(vxISDataTimeSteps(EdgeSensorInd,ii).^2 + vyISDataTimeSteps(EdgeSensorInd,ii).^2),'o')
            ylim([0,600])
            pause(0.1)
        end
    end
    T11ISamplesDataTimeSteps(:,n) = T11ISDataTimeSteps(:);
    T22ISamplesDataTimeSteps(:,n) = T22ISDataTimeSteps(:);
    T12ISamplesDataTimeSteps(:,n) = T12ISDataTimeSteps(:);
    vxISamplesDataTimeSteps(:,n) = vxISDataTimeSteps(:);
    vyISamplesDataTimeSteps(:,n) = vyISDataTimeSteps(:);
end
if RunOptions.FullqVectorData == 1;
    clear T11ISTimeSteps T22ISTimeSteps T12ISTimeSteps T11ISDataTimeSteps T22ISDataTimeSteps T12ISDataTimeSteps
end
clear vxISTimeSteps vyISTimeSteps  vxISDataTimeSteps vyISDataTimeSteps

%% =======================================================================%
%                         Approximation Error
%=========================================================================%
errorT11 = T11ASamplesDataTimeSteps - T11ISamplesDataTimeSteps;
errorT22 = T22ASamplesDataTimeSteps - T22ISamplesDataTimeSteps;
errorT12 = T12ASamplesDataTimeSteps - T12ISamplesDataTimeSteps;
errorvx = vxASamplesDataTimeSteps - vxISamplesDataTimeSteps;
errorvy = vyASamplesDataTimeSteps - vyISamplesDataTimeSteps;

clear T11ASamplesDataTimeSteps T11ISamplesDataTimeSteps T22ASamplesDataTimeSteps T22ISamplesDataTimeSteps T12ASamplesDataTimeSteps T12ISamplesDataTimeSteps vxASamplesDataTimeSteps vxISamplesDataTimeSteps vyASamplesDataTimeSteps vyISamplesDataTimeSteps

printf('Computing sample means');
SampleMeanT11 = (1/N_Samples)*sum(errorT11,2);
SampleMeanT22 = (1/N_Samples)*sum(errorT22,2);
SampleMeanT12 = (1/N_Samples)*sum(errorT12,2);
SampleMeanvx = (1/N_Samples)*sum(errorvx,2);
SampleMeanvy = (1/N_Samples)*sum(errorvy,2);

printf('Computing sample covariances');
if RunOptions.FullqVectorData == 1;
    SampleCovT11 = 1/(N_Samples - 1)*(errorT11*errorT11') - (N_Samples/(N_Samples - 1))*(SampleMeanT11*SampleMeanT11');
    SampleCovT22 = 1/(N_Samples - 1)*(errorT22*errorT22') - (N_Samples/(N_Samples - 1))*(SampleMeanT22*SampleMeanT22');
    SampleCovT12 = 1/(N_Samples - 1)*(errorT12*errorT12') - (N_Samples/(N_Samples - 1))*(SampleMeanT12*SampleMeanT12');
else
    SampleCovT11 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
    SampleCovT22 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
    SampleCovT12 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
end
SampleCovvx = 1/(N_Samples - 1)*(errorvx*errorvx') - (N_Samples/(N_Samples - 1))*(SampleMeanvx*SampleMeanvx');
SampleCovvy = 1/(N_Samples - 1)*(errorvy*errorvy') - (N_Samples/(N_Samples - 1))*(SampleMeanvy*SampleMeanvy');

%=== Reshaping sample means to be used later ===%
SampleMeanT11 = reshape(SampleMeanT11,NumberofSensors,NumberofTimeSteps);
SampleMeanT22 = reshape(SampleMeanT22,NumberofSensors,NumberofTimeSteps);
SampleMeanT12 = reshape(SampleMeanT12,NumberofSensors,NumberofTimeSteps);
SampleMeanvx = reshape(SampleMeanvx,NumberofSensors,NumberofTimeSteps);
SampleMeanvy = reshape(SampleMeanvy,NumberofSensors,NumberofTimeSteps);

printf('Computation of sample means and covariances complete');

