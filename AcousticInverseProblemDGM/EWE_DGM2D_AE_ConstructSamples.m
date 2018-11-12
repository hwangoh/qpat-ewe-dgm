function [hAS_FEM,hIS_FEM,errorT11,errorT22,errorT12,errorvx,errorvy]=EWE_DGM2D_AE_ConstructSamples(RunOptions,MeshANodes,MeshAElements,MeshADimensns,MeshADomainIndices,MeshINodes,MeshIElements,xA,yA,NpA,KA,PrecomputedIntrplteObjectsA,pinfoA,DGMMeshANorder,DGMMeshAlambda,DGMMeshArho,xI,yI,NpI,KI,rhoI,PrecomputedIntrplteObjectsI,pinfoI,SensorsA,SensorsI,dt,Prior,PLOT)

% EWE_DGM2D_AE_ConstructSamples computes the samples required for calculating 
% the approximation error mean and approximation error covariance
%
% Inputs:
%   RunOptions:
%              NumberofTimeSteps
%              AEN_Samples - number of samples we wish to use to calculate approximation error
%   MeshANodes, MeshINodes - Number of FEM Nodes by 2 array storing the coordinates of the nodes
%   MeshAElements, MeshIElements - Number of Elements by 3 array where the ith row contains the indices of the nodes in the ith element
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
%   dt - size of time steps
%   Prior:
%      corr_h - 2*1 vector = [correlation_x, correlation_y]. The larger the number,
%          the more points near to a marginalisation point are correlated
%          in the x and y direction.
%      bounds_h - approximately [h_min, h_max]
%      Exp_h - initial guess of h
%   PLOT - To plot or not to plot, that is the question
%
% Outputs:
%      hAS_FEM, hIS_FEM - Prior model sample draws, mainly used for QR
%                           construction of approximation error statistics
%      errorTij, errorvi - Sample errors
%
% Hwan Goh 15/01/2018, University of Auckland, New Zealand
%          27/07/2018 - huge overhaul of approximation error codes

PLOT.DGMForward = 0; %Suppress plotting of wave propagation
PLOT.PriorMargPoints = 0; %Suppress plotting of marginalisation points
PLOT.PriorSamples = 0; %Suppress plotting of prior samples
PLOT.DGMForwardSensorData = 0; %Suppress plotting of sensory data
N_Samples = RunOptions.AEN_Samples;
NumberofSensors = size(SensorsI,2);
NumberofTimeSteps = RunOptions.NumberofTimeSteps;
printf([num2str(N_Samples) ' samples to be computed for approximation error']);

if RunOptions.GenerateAndSave_AESamples == 1
%=========================================================================%
%                         Accurate Model Samples
%=========================================================================%
printf('Computing samples of accurate model');

%=== Absorbed Energy Density Samples ===%
if RunOptions.UseInformSmoothPrior == 1;
    [~,~,~,~,hAS_FEM]=SmoothnessPrior_Informative(MeshANodes,MeshAElements,2*MeshADimensns,Prior.InformSmooth_Bounds_h,Prior.InformSmooth_Corr_h,Prior.Exp_h,Prior.InformSmooth_Normalize,N_Samples,PLOT);
end
if RunOptions.UseACPrior == 1;
    [~,~,~,hAS_FEM] = SmoothnessPrior_AutoCorr(MeshANodes,Prior.Exp_h,Prior.AC_Var_h,Prior.AC_Corr_h,N_Samples,PLOT);
end

%=== Solid Layer Samples ===%
if RunOptions.FluidDomainWithSolidLayerMeshD == 1
    if RunOptions.AE_VarySolidLayer == 1 
%         [muAS_FEM,rhoAS_FEM] = SmoothnessPrior_AutoCorr_EL(RunOptions,MeshAElements,MeshANodes,MeshADomainIndices,Prior.Exp_rho_e,Prior.AC_Var_rho_e,Prior.AC_Corr_rho_e,Prior.AC_ExpShiftExp_rho_e,Prior.AC_ExpShiftVar_rho_e,Prior.Exp_cs,Prior.AC_Var_cs,Prior.AC_Corr_cs,Prior.AC_ExpShiftExp_cs,Prior.AC_ExpShiftVar_cs,N_Samples,PLOT);
    end
end

%=== Forward Data Samples ===%
if RunOptions.FullqVectorData == 1;
    T11ASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
    T22ASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
    T12ASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
else
    T11ASamplesDataTimeSteps = sparse(NumberofSensors*NumberofTimeSteps,N_Samples);
    T22ASamplesDataTimeSteps = sparse(NumberofSensors*NumberofTimeSteps,N_Samples);
    T12ASamplesDataTimeSteps = sparse(NumberofSensors*NumberofTimeSteps,N_Samples);
end
vxASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
vyASamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);

PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshA; %This is for plotting forward elastic wave propagation on inversion mesh. FwrdFunction is called which is EWE_DGM2D_LSExpRK4. However, that function calls PLOT.TRI_DGMMeshD when plotting, so here we replace it with DGMMeshA but re-use the same name

%=== Accurate Forward Data ===%
for n=1:N_Samples
    printf(['\nComputing Accurate Model Sample ' num2str(n) ' of ' num2str(N_Samples)]);
    printf(['For the Case ' RunOptions.SaveFileNameAESamples]);
    if RunOptions.FluidDomainWithSolidLayerMeshD == 1
        if RunOptions.AE_VarySolidLayer == 1
%             [DGMMeshAmu, ~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),muAS_FEM(:,n),xA,yA,NpA*KA,PrecomputedIntrplteObjectsA);
%             [DGMMeshArho, ~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),rhoAS_FEM(:,n),xA,yA,NpA*KA,PrecomputedIntrplteObjectsA);
            ElasticMediumDensity = 1850 + 20*rand(1);
            ElasticMediumWaveSpeed_cs = 1500 + 20*rand(1);
            ElasticLamemu = ElasticMediumDensity*ElasticMediumWaveSpeed_cs^2;
            DGMMeshAmu = zeros(NpA,KA);
            DGMMeshArho = zeros(NpA,KA);
            for ii=1:KA
                if MeshADomainIndices(ii) == 1 || MeshADomainIndices(ii) == 3
                    DGMMeshAmu(:,ii) = 0;
                    DGMMeshArho(:,ii) = RunOptions.AcousticMediumDensity;
                else
                    DGMMeshAmu(:,ii) = ElasticLamemu;
                    DGMMeshArho(:,ii) = ElasticMediumDensity;
                end
            end
            pinfoA = EWE_DGM2D_PrecomputeUpwindFluxPNonConf(RunOptions,pinfoA,DGMMeshANorder,DGMMeshArho,DGMMeshAlambda,DGMMeshAmu,MeshADomainIndices,RunOptions.FluidDomainWithSolidLayerMeshD,RunOptions.SolidMeshD);
        end
    end
    [hAS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),hAS_FEM(:,n),xA,yA,NpA*KA,PrecomputedIntrplteObjectsA);
    p0AS_DGM = RunOptions.LinearHeatExpansionCoeff*reshape(hAS_DGM,NpA,KA)./(DGMMeshArho*RunOptions.SpecificHeatCoeff);
    if RunOptions.TimeLSERK4 == 1;
        [T11ASTimeSteps,T22ASTimeSteps,T12ASTimeSteps,vxASTimeSteps,vyASTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,p0AS_DGM,xA,yA,NpA,KA,pinfoA,dt,PLOT);
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
    if RunOptions.FullqVectorData == 1;
        T11ASamplesDataTimeSteps(:,n) = T11ASDataTimeSteps(:);
        T22ASamplesDataTimeSteps(:,n) = T22ASDataTimeSteps(:);
        T12ASamplesDataTimeSteps(:,n) = T12ASDataTimeSteps(:);
    end
    vxASamplesDataTimeSteps(:,n) = vxASDataTimeSteps(:);
    vyASamplesDataTimeSteps(:,n) = vyASDataTimeSteps(:);
end
if RunOptions.FullqVectorData == 1;
clear T11ASTimeSteps T22ASTimeSteps T12ASTimeSteps T11ASDataTimeSteps T22ASDataTimeSteps T12ASDataTimeSteps
end
clear vxASTimeSteps vyASTimeSteps vxASDataTimeSteps vyASDataTimeSteps

%=== Saving and Loading Accurate Model Samples ===% this is reusable because most of the time we consider different inaccurate models; so now we don't need to regenerate accurate samples. But it's simpler to comment the save/load in/out than add a boolean
save(RunOptions.SaveFileNameAESamples,'hAS_FEM','T11ASamplesDataTimeSteps','T22ASamplesDataTimeSteps','T12ASamplesDataTimeSteps','vxASamplesDataTimeSteps','vyASamplesDataTimeSteps','-v7.3')
save(RunOptions.SaveFileNameAESamplesAccurate,'hAS_FEM','T11ASamplesDataTimeSteps','T22ASamplesDataTimeSteps','T12ASamplesDataTimeSteps','vxASamplesDataTimeSteps','vyASamplesDataTimeSteps','-v7.3')
% load(RunOptions.SaveFileNameAESamplesAccurate,'hAS_FEM','T11ASamplesDataTimeSteps','T22ASamplesDataTimeSteps','T12ASamplesDataTimeSteps','vxASamplesDataTimeSteps','vyASamplesDataTimeSteps')

%=========================================================================%
%                       Less Accurate Model Samples
%=========================================================================%
printf('Computing samples of less accurate model');
%=== Initial Pressure Samples ===%
[hIS_FEM(:,1), PrecomputedIntrplteObjectsAtoIFEM] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),hAS_FEM(:,1),MeshINodes(:,1),MeshINodes(:,2),size(MeshINodes,1),0);
for n=2:N_Samples %Just so we do not need to compute the interpolation objects over and over again
    [hIS_FEM(:,n), ~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),hAS_FEM(:,n),MeshINodes(:,1),MeshINodes(:,2),size(MeshINodes,1),PrecomputedIntrplteObjectsAtoIFEM);
end
if PLOT.PriorSamples==1;
    figure
    FemPlot2D(MeshINodes,MeshIElements,hIS_FEM(:,1:min(N_Samples,10)))
    title(PLOT.Figure_Prior_Title,'FontWeight','bold')
end

%=== Forward Data Samples ===%
if RunOptions.FullqVectorData == 1;
    T11ISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
    T22ISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
    T12ISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
else
    T11ISamplesDataTimeSteps = sparse(NumberofSensors*NumberofTimeSteps,N_Samples);
    T22ISamplesDataTimeSteps = sparse(NumberofSensors*NumberofTimeSteps,N_Samples);
    T12ISamplesDataTimeSteps = sparse(NumberofSensors*NumberofTimeSteps,N_Samples);
end
vxISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);
vyISamplesDataTimeSteps = zeros(NumberofSensors*NumberofTimeSteps,N_Samples);

PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshI; %This is for plotting forward elastic wave propagation on inversion mesh. FwrdFunction is called which is EWE_DGM2D_LSExpRK4. However, that function calls PLOT.TRI_DGMMeshD when plotting, so here we replace it with DGMMeshI but re-use the same name

%=== Less Accurate Forward Data ===%
NumberofSensors = size(SensorsI,2);
for n=1:N_Samples
    printf(['\nComputing (SVD) Less Accurate Model Sample ' num2str(n) ' of ' num2str(N_Samples)]);
    printf(['For the Case ' RunOptions.SaveFileNameAESamples]);
    [hIS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshIElements,1),MeshINodes*(1/RunOptions.ScalingOptical),hIS_FEM(:,n),xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
    p0IS_DGM = RunOptions.LinearHeatExpansionCoeff*reshape(hIS_DGM,NpI,KI)./(rhoI*RunOptions.SpecificHeatCoeff);
    if RunOptions.TimeLSERK4 == 1;
        [T11ISTimeSteps,T22ISTimeSteps,T12ISTimeSteps,vxISTimeSteps,vyISTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,p0IS_DGM,xI,yI,NpI,KI,pinfoI,dt,PLOT);
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
    if RunOptions.FullqVectorData == 1;
        T11ISamplesDataTimeSteps(:,n) = T11ISDataTimeSteps(:);
        T22ISamplesDataTimeSteps(:,n) = T22ISDataTimeSteps(:);
        T12ISamplesDataTimeSteps(:,n) = T12ISDataTimeSteps(:);
    end
    vxISamplesDataTimeSteps(:,n) = vxISDataTimeSteps(:);
    vyISamplesDataTimeSteps(:,n) = vyISDataTimeSteps(:);
end
if RunOptions.FullqVectorData == 1;
    clear T11ISTimeSteps T22ISTimeSteps T12ISTimeSteps T11ISDataTimeSteps T22ISDataTimeSteps T12ISDataTimeSteps
end
clear vxISTimeSteps vyISTimeSteps  vxISDataTimeSteps vyISDataTimeSteps

%=========================================================================%
%                               Error Samples
%=========================================================================%
errorT11 = T11ASamplesDataTimeSteps(:,1:N_Samples) - T11ISamplesDataTimeSteps(:,1:N_Samples);
errorT22 = T22ASamplesDataTimeSteps(:,1:N_Samples) - T22ISamplesDataTimeSteps(:,1:N_Samples);
errorT12 = T12ASamplesDataTimeSteps(:,1:N_Samples) - T12ISamplesDataTimeSteps(:,1:N_Samples);
errorvx = vxASamplesDataTimeSteps(:,1:N_Samples) - vxISamplesDataTimeSteps(:,1:N_Samples);
errorvy = vyASamplesDataTimeSteps(:,1:N_Samples) - vyISamplesDataTimeSteps(:,1:N_Samples);

printf('Saving AE Samples');
save(RunOptions.SaveFileNameAESamples,'hAS_FEM','T11ASamplesDataTimeSteps','T22ASamplesDataTimeSteps','T12ASamplesDataTimeSteps','vxASamplesDataTimeSteps','vyASamplesDataTimeSteps',...
                           'hIS_FEM','T11ISamplesDataTimeSteps','T22ISamplesDataTimeSteps','T12ISamplesDataTimeSteps','vxISamplesDataTimeSteps','vyISamplesDataTimeSteps',...
                           'errorT11','errorT22','errorT12','errorvx','errorvy',...
                           'MeshANodes','MeshAElements','MeshADimensns','PrecomputedIntrplteObjectsA','pinfoA','SensorsA','xA','yA','NpA','KA',...
                           '-v7.3')
else
printf('Loading AE Samples');
load(RunOptions.SaveFileNameAESamples,'hAS_FEM','hIS_FEM',...
                                      'errorT11','errorT22','errorT12','errorvx','errorvy')
end


