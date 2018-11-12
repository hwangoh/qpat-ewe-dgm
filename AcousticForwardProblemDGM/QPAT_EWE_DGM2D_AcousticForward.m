% QPAT_EWE_DGM2D_AcousticForward solves the acoustic forward problem using the
% discontinuous Galerkin method
%
% Hwan Goh, University of Auckland, New Zealand - 8/7/2015
% Last Edited: 16/11/2017 - removed reliance on Globals2D and save everything into structures

disp(' ')
disp('-----------------')
disp('Generating Data')
disp('-----------------')

%% =======================================================================%
%                            Interpolation
%=========================================================================%
DataVrblsWave.h = DataVrblsOptical.h;

%=== Histogram of True p0 Values ===%
% [Nbins,Nc]=hist(DataVrblsWave.p0,50);
% Ndx=Nc(2)-Nc(1);
% figure
% bar(Nc, Nbins)
% 
% keyboard

%=== Interpolate Absorbed Energy Density ===%
[hDGM, PrecomputedIntrplteObjectsD] = IntrplteOver2DTriangulatedMesh(MeshD.N_Elm,MeshD.Nodes*(1/RunOptions.ScalingOptical),DataVrblsWave.h,DGMMeshD.x,DGMMeshD.y,DGMMeshD.Np*DGMMeshD.K,PrecomputedIntrplteObjectsD);
if RunOptions.InverseCrime == 1;
    PrecomputedIntrplteObjectsI = PrecomputedIntrplteObjectsD; %In inverse crime case, PrecompIntrplteObjectsI was not previously computed in AcoustForwardDGMSetup
end

%=== Plotting Interpolated Absorbed Energy Density ===%
PLOT.TRI_DGMMeshD=delaunay(DGMMeshD.x,DGMMeshD.y); %To be used later when plotting wave propagation
if PLOT.IntrplteAbsorbedEnergyDensity==1;
    figure(PLOT.Figure_AbsorbedEnergyDensityInterpolated)
    trisurf(PLOT.TRI_DGMMeshD,DGMMeshD.x,DGMMeshD.y,hDGM);
    if PLOT.DGMPlotBirdsEyeView == 1;
        view(2)
    end
    shading interp %thanks Ru!
    if PLOT.DGMPlotUsezLim == 1;
        zlim(PLOT.AbsorbedEnergyDensityzAxis);
    end
    colorbar
    colormap(jet(256))
    title(PLOT.Figure_AbsorbedEnergyDensityInterpolated_Title,'FontWeight','bold')
end

%=== Setting Up Plotting on Inversion Mesh ===%
if RunOptions.InverseCrime == 1
    PLOT.TRI_DGMMeshI = PLOT.TRI_DGMMeshD;
else
    PLOT.TRI_DGMMeshI = delaunay(DGMMeshI.x,DGMMeshI.y);
end

%=== Interpolate Absorbed Energy Density ===%
if RunOptions.StrainVelocityForm == 1
    DataVrblsWave.Gruneisen = RunOptions.LinearHeatExpansionCoeff;
end
if RunOptions.StressVelocityForm == 1
    DataVrblsWave.Gruneisen = -(2*RunOptions.ElasticLamemu_e + 3*RunOptions.ElasticLamelambda_e);
end
DataVrblsWave.p0DGM = DataVrblsWave.Gruneisen*reshape(hDGM,DGMMeshD.Np,DGMMeshD.K)./(DGMMeshD.rho*RunOptions.SpecificHeatCoeff);

%=== Set Initial Pressure on Optical Forward Mesh, to be Interpolated ===%
if RunOptions.InterplteGaussianp0 == 1;
    DataVrblsWave.p0 = exp(-(RunOptions.GaussianWidth)*((MeshD.Nodes(:,1)*(1/RunOptions.ScalingOptical)).^2+(MeshD.Nodes(:,2)*(1/RunOptions.ScalingOptical)).^2));
    [DataVrblsWave.p0DGM, PrecomputedIntrplteObjectsD] = IntrplteOver2DTriangulatedMesh(MeshD.N_Elm,MeshD.Nodes*(1/RunOptions.ScalingOptical),DataVrblsWave.p0,DGMMeshD.x,DGMMeshD.y,DGMMeshD.Np*DGMMeshD.K,PrecomputedIntrplteObjectsD);
    PLOT.AbsorbedEnergyDensityzAxis = [0 1]; %z axis limits for plotting QPAT elastic wave propogation
    PLOT.DGMPlotzAxis = [0 1]; %z axis limits for plotting QPAT elastic wave propogation
    PLOT.DGMAdjointPlotzAxis = [0 1]; %z axis limits for plotting QPAT elastic wave propogation
    PLOT.HoldColourAxis = 1; %Hold colour axis
    PLOT.ColourAxis = [0 1]; %Colour axis
end

%=== Test Case on DGM Mesh ===%
if RunOptions.p0Test1 == 1 || RunOptions.p0Test2 == 1;
    AcousticForwardTestCases
end


%% =======================================================================%
%                       Precomputing Flux Terms
%=========================================================================%
%=== Precompute Some Upwind Flux Terms for Non-Conforming Mesh ===%
if RunOptions.UseMyRectangularMesh == 1
    MeshD.DomainIndices = 'Not constructed'; %Domain Indices
end
DGMMeshD.pinfo = EWE_DGM2D_PrecomputeUpwindFluxPNonConf(RunOptions,DGMMeshD.pinfo,DGMMeshD.Norder,DGMMeshD.rho,DGMMeshD.lambda,DGMMeshD.mu,MeshD.DomainIndices,RunOptions.FluidDomainWithSolidLayerMeshD,RunOptions.SolidMeshD);             

%% =======================================================================%
%                            Time Stepping
%=========================================================================%
if RunOptions.LoadNoisyData == 0
%=== Low-Storage Five-Stage Explicit 4th Order Runge-Kutta Method ===%    
if RunOptions.TimeLSERK4 == 1;
    [T11TimeSteps,T22TimeSteps,T12TimeSteps,vxTimeSteps,vyTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,DataVrblsWave.p0DGM,DGMMeshD.x,DGMMeshD.y,DGMMeshD.Np,DGMMeshD.K,DGMMeshD.pinfo,dt,PLOT);
end

%% =======================================================================%
%                            Sensory Data
%=========================================================================%
DataVrblsWave.T11DataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
DataVrblsWave.T22DataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
DataVrblsWave.T12DataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
DataVrblsWave.vxDataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
DataVrblsWave.vyDataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
%=== Interpolation to Form Sensory Data ===%
for t=1:RunOptions.NumberofTimeSteps
    for s=1:DataVrblsWave.NumberofSensors
        if RunOptions.FullqVectorData == 1;
            DataVrblsWave.T11DataTimeSteps(s,t) = DataVrblsWave.SensorsD{s}.l_iatsensor*T11TimeSteps(DataVrblsWave.SensorsD{s}.id,t);
            DataVrblsWave.T22DataTimeSteps(s,t) = DataVrblsWave.SensorsD{s}.l_iatsensor*T22TimeSteps(DataVrblsWave.SensorsD{s}.id,t);
            DataVrblsWave.T12DataTimeSteps(s,t) = DataVrblsWave.SensorsD{s}.l_iatsensor*T12TimeSteps(DataVrblsWave.SensorsD{s}.id,t);
        end
        DataVrblsWave.vxDataTimeSteps(s,t) = DataVrblsWave.SensorsD{s}.l_iatsensor*vxTimeSteps(DataVrblsWave.SensorsD{s}.id,t);
        DataVrblsWave.vyDataTimeSteps(s,t) = DataVrblsWave.SensorsD{s}.l_iatsensor*vyTimeSteps(DataVrblsWave.SensorsD{s}.id,t);
    end
end

clear T11TimeSteps T22TimeSteps T12TimeSteps vxTimeSteps vyTimeSteps

if RunOptions.AddNoise == 1
    %=== Adding Noise - Max minus Min of Data ===%
    if RunOptions.NoiseMinMax == 1;
        if RunOptions.FullqVectorData == 1;
            MaxMinusMinT11 = max(DataVrblsWave.T11DataTimeSteps(:)) - min(DataVrblsWave.T11DataTimeSteps(:));
            MaxMinusMinT22 = max(DataVrblsWave.T22DataTimeSteps(:)) - min(DataVrblsWave.T22DataTimeSteps(:));
            MaxMinusMinT12 = max(DataVrblsWave.T12DataTimeSteps(:)) - min(DataVrblsWave.T12DataTimeSteps(:));
        end
        MaxMinusMinvx = max(DataVrblsWave.vxDataTimeSteps(:)) - min(DataVrblsWave.vxDataTimeSteps(:));
        MaxMinusMinvy = max(DataVrblsWave.vyDataTimeSteps(:)) - min(DataVrblsWave.vyDataTimeSteps(:));
        for s = 1:DataVrblsWave.NumberofSensors
            for t = 1:RunOptions.NumberofTimeSteps
                if RunOptions.FullqVectorData == 1;
                    DataVrblsWave.T11DataTimeSteps(s,t) = DataVrblsWave.T11DataTimeSteps(s,t) + RunOptions.NoiseLevel*MaxMinusMinT11*randn(1);
                    DataVrblsWave.T22DataTimeSteps(s,t) = DataVrblsWave.T22DataTimeSteps(s,t) + RunOptions.NoiseLevel*MaxMinusMinT22*randn(1);
                    DataVrblsWave.T12DataTimeSteps(s,t) = DataVrblsWave.T12DataTimeSteps(s,t) + RunOptions.NoiseLevel*MaxMinusMinT12*randn(1);
                end
                DataVrblsWave.vxDataTimeSteps(s,t) = DataVrblsWave.vxDataTimeSteps(s,t) + RunOptions.NoiseLevel*MaxMinusMinvx*randn(1);
                DataVrblsWave.vyDataTimeSteps(s,t) = DataVrblsWave.vyDataTimeSteps(s,t) + RunOptions.NoiseLevel*MaxMinusMinvy*randn(1);
            end
        end
    end
    %=== Adding Noise - Max minus Min of Data at Each Sensor ===%
    if RunOptions.NoiseMinMaxS == 1;
        DataVrblsWave.MaxMinusMinST11 = zeros(1,DataVrblsWave.NumberofSensors);
        DataVrblsWave.MaxMinusMinST22 = zeros(1,DataVrblsWave.NumberofSensors);
        DataVrblsWave.MaxMinusMinST12 = zeros(1,DataVrblsWave.NumberofSensors);
        DataVrblsWave.MaxMinusMinSvx = zeros(1,DataVrblsWave.NumberofSensors);
        DataVrblsWave.MaxMinusMinSvy = zeros(1,DataVrblsWave.NumberofSensors);
        for s = 1:DataVrblsWave.NumberofSensors
            if RunOptions.FullqVectorData == 1;
                DataVrblsWave.MaxMinusMinST11(s) = max(DataVrblsWave.T11DataTimeSteps(s,:)) - min(DataVrblsWave.T11DataTimeSteps(s,:));
                DataVrblsWave.MaxMinusMinST22(s) = max(DataVrblsWave.T22DataTimeSteps(s,:)) - min(DataVrblsWave.T22DataTimeSteps(s,:));
                DataVrblsWave.MaxMinusMinST12(s) = max(DataVrblsWave.T12DataTimeSteps(s,:)) - min(DataVrblsWave.T12DataTimeSteps(s,:));
            end
            DataVrblsWave.MaxMinusMinSvx(s) = max(DataVrblsWave.vxDataTimeSteps(s,:)) - min(DataVrblsWave.vxDataTimeSteps(s,:));
            DataVrblsWave.MaxMinusMinSvy(s) = max(DataVrblsWave.vyDataTimeSteps(s,:)) - min(DataVrblsWave.vyDataTimeSteps(s,:));           
            for t = 1:RunOptions.NumberofTimeSteps
                if RunOptions.FullqVectorData == 1;
                    DataVrblsWave.T11DataTimeSteps(s,t) = DataVrblsWave.T11DataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxMinusMinST11(s)*randn(1);
                    DataVrblsWave.T22DataTimeSteps(s,t) = DataVrblsWave.T22DataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxMinusMinST22(s)*randn(1);
                    DataVrblsWave.T12DataTimeSteps(s,t) = DataVrblsWave.T12DataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxMinusMinST12(s)*randn(1);
                end
                DataVrblsWave.vxDataTimeSteps(s,t) = DataVrblsWave.vxDataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxMinusMinSvx(s)*randn(1);
                DataVrblsWave.vyDataTimeSteps(s,t) = DataVrblsWave.vyDataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxMinusMinSvy(s)*randn(1);
            end
        end
    end
    %=== Adding Noise - Max of Data ===%
    if RunOptions.NoiseMax == 1;
        if RunOptions.FullqVectorData == 1;
            MaxT11 = max(abs(DataVrblsWave.T11DataTimeSteps(:)));
            MaxT22 = max(abs(DataVrblsWave.T22DataTimeSteps(:)));
            MaxT12 = max(abs(DataVrblsWave.T12DataTimeSteps(:)));
        end
        Maxvx = max(abs(DataVrblsWave.vxDataTimeSteps(:)));
        Maxvy = max(abs(DataVrblsWave.vyDataTimeSteps(:)));
        for s = 1:DataVrblsWave.NumberofSensors
            for t = 1:RunOptions.NumberofTimeSteps
                if RunOptions.FullqVectorData == 1;
                    DataVrblsWave.T11DataTimeSteps(s,t) = DataVrblsWave.T11DataTimeSteps(s,t) + RunOptions.NoiseLevel*MaxT11*randn(1);
                    DataVrblsWave.T22DataTimeSteps(s,t) = DataVrblsWave.T22DataTimeSteps(s,t) + RunOptions.NoiseLevel*MaxT22*randn(1);
                    DataVrblsWave.T12DataTimeSteps(s,t) = DataVrblsWave.T12DataTimeSteps(s,t) + RunOptions.NoiseLevel*MaxT12*randn(1);
                end
                DataVrblsWave.vxDataTimeSteps(s,t) = DataVrblsWave.vxDataTimeSteps(s,t) + RunOptions.NoiseLevel*Maxvx*randn(1);
                DataVrblsWave.vyDataTimeSteps(s,t) = DataVrblsWave.vyDataTimeSteps(s,t) + RunOptions.NoiseLevel*Maxvy*randn(1);
            end
        end
    end
    %=== Adding Noise - Max of Data At Each Sensor ===%
    if RunOptions.NoiseMaxS == 1;
        DataVrblsWave.MaxST11 = zeros(1,DataVrblsWave.NumberofSensors);
        DataVrblsWave.MaxST22 = zeros(1,DataVrblsWave.NumberofSensors);
        DataVrblsWave.MaxST12 = zeros(1,DataVrblsWave.NumberofSensors);
        DataVrblsWave.MaxSvx = zeros(1,DataVrblsWave.NumberofSensors);
        DataVrblsWave.MaxSvy = zeros(1,DataVrblsWave.NumberofSensors);
        for s = 1:DataVrblsWave.NumberofSensors
            if RunOptions.FullqVectorData == 1;
                DataVrblsWave.MaxST11(s) = max(abs(DataVrblsWave.T11DataTimeSteps(s,:)));
                DataVrblsWave.MaxST22(s) = max(abs(DataVrblsWave.T22DataTimeSteps(s,:)));
                DataVrblsWave.MaxST12(s) = max(abs(DataVrblsWave.T12DataTimeSteps(s,:)));
            end
            DataVrblsWave.MaxSvx(s) = max(abs(DataVrblsWave.vxDataTimeSteps(s,:)));
            DataVrblsWave.MaxSvy(s) = max(abs(DataVrblsWave.vyDataTimeSteps(s,:)));
            for t = 1:RunOptions.NumberofTimeSteps
                if RunOptions.FullqVectorData == 1;
                    DataVrblsWave.T11DataTimeSteps(s,t) = DataVrblsWave.T11DataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxST11(s)*randn(1);
                    DataVrblsWave.T22DataTimeSteps(s,t) = DataVrblsWave.T22DataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxST22(s)*randn(1);
                    DataVrblsWave.T12DataTimeSteps(s,t) = DataVrblsWave.T12DataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxST12(s)*randn(1);
                end
                DataVrblsWave.vxDataTimeSteps(s,t) = DataVrblsWave.vxDataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxSvx(s)*randn(1);
                DataVrblsWave.vyDataTimeSteps(s,t) = DataVrblsWave.vyDataTimeSteps(s,t) + RunOptions.NoiseLevel*DataVrblsWave.MaxSvy(s)*randn(1);
            end
        end
    end
end

%% =======================================================================%
%                  Loading and Saving Sensory Data
%=========================================================================%
%=== Saving Sensor Data ===%
T11DataTimeSteps = DataVrblsWave.T11DataTimeSteps;
T22DataTimeSteps = DataVrblsWave.T22DataTimeSteps;
T12DataTimeSteps = DataVrblsWave.T12DataTimeSteps;
vxDataTimeSteps = DataVrblsWave.vxDataTimeSteps;
vyDataTimeSteps = DataVrblsWave.vyDataTimeSteps;
RunOptions.SaveFileNameData = sprintf('%s-%s-%sNoise%s-%sD-%dSensors-%sFinalTime',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDataDomain,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.TrelisMeshDElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime); 
% save(RunOptions.SaveFileNameData,'T11DataTimeSteps','T22DataTimeSteps','T12DataTimeSteps','vxDataTimeSteps','vyDataTimeSteps','-v7.3')
end

if RunOptions.LoadNoisyData == 1
%=== Loading Sensor Data ===%
disp('Loading Data')
RunOptions.SaveFileNameData = sprintf('%s-%s-%sNoise%s-%sD-%dSensors-%sFinalTime',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDataDomain,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.TrelisMeshDElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime); 
load(RunOptions.SaveFileNameData,'T11DataTimeSteps','T22DataTimeSteps','T12DataTimeSteps','vxDataTimeSteps','vyDataTimeSteps')
DataVrblsWave.T11DataTimeSteps = T11DataTimeSteps;
DataVrblsWave.T22DataTimeSteps = T22DataTimeSteps;
DataVrblsWave.T12DataTimeSteps = T12DataTimeSteps;
DataVrblsWave.vxDataTimeSteps = vxDataTimeSteps;
DataVrblsWave.vyDataTimeSteps = vyDataTimeSteps;
end
clear T11DataTimeSteps T22DataTimeSteps T12DataTimeSteps vxDataTimeSteps vyDataTimeSteps

%% =======================================================================%
%                           Plotting Sensor Data
%=========================================================================%
DataVrblsWave.EdgeSensorInd = find(DataVrblsWave.SensorCoords(:,1)==-0.01);
if RunOptions.UseFullDomainData ~= 1 && PLOT.DGMForwardSensorData == 1;
    figure(PLOT.Figure_WavePropagation_SensorData)
    for t=1:RunOptions.NumberofTimeSteps
        plot(1:1:length(DataVrblsWave.EdgeSensorInd),sqrt(DataVrblsWave.vxDataTimeSteps(DataVrblsWave.EdgeSensorInd,t).^2 + DataVrblsWave.vyDataTimeSteps(DataVrblsWave.EdgeSensorInd,t).^2),'o')
        ylim([0,600])
        title(PLOT.Figure_WavePropagation_SensorData_Title,'FontWeight','bold')
        pause(0.1)
    end
    % SelectedSensor = 52; %For 20 sensors
    SelectedSensor = 17; %For 10 sensors
    plot(dt:dt:RunOptions.NumberofTimeSteps*dt,sqrt(DataVrblsWave.vxDataTimeSteps(SelectedSensor,:).^2 + DataVrblsWave.vyDataTimeSteps(SelectedSensor,:).^2),'r')
    WavePropAtSensorTitle = sprintf('Sample Wave Propagation at Sensor %d',SelectedSensor);
    title(WavePropAtSensorTitle,'FontWeight','bold')
end