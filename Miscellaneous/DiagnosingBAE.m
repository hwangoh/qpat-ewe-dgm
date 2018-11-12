close all
clc

%Start by running QPAT_EWE_DGM2D_Driver up to keyboard after "Display Selected Options" section

%=========================================================================%
%                           Loading Objects
%=========================================================================%
addpath(genpath('C:\Users\hgoh009\Dropbox\Work\PhD\PhD_Stuff\Codes\Outputs'))
addpath(genpath('D:\Dropbox\Work\PhD\PhD_Stuff\Codes\Outputs'))

%=== Load Mesh and Parameters ===%
RunOptions.SaveFileNameMeshAndParameters = sprintf('MeshAndParameters-%sD-%sI-%dSensors-%sFinalTime',RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime);
load(RunOptions.SaveFileNameMeshAndParameters);
RunOptions.NumberofTimeSteps = ceil(RunOptions.FinalTime/dt); 
close all

%=== Load Reconstructions ===%
RunOptions.SaveFileNameReconstructions = sprintf('Reconstructions-%s-%s-%sNoise-%sD-%sI-%dSensors-%sFinalTime-%s-%s%s-%sCorr-%dAESamples-%sAE',RunOptions.SaveFileNameDomain,RunOptions.SaveFileNameStateVector,RunOptions.SaveFileNameNoiseLevel,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.SaveFileNameExpT,RunOptions.SaveFileNameCov_ENoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.SaveFileNamePriorCorr,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
load(RunOptions.SaveFileNameReconstructions);

%=== Load AE Samples SVD ===%
RunOptions.SaveFileNameAESamples = sprintf('AESamples-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
if RunOptions.UseACPrior == 1
    SaveFileNameAESamples = sprintf('AESamples-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE-ACPrior',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
end
load(RunOptions.SaveFileNameAESamples,'p0AS_FEM','T11ASamplesDataTimeSteps','T22ASamplesDataTimeSteps','T12ASamplesDataTimeSteps','vxASamplesDataTimeSteps','vyASamplesDataTimeSteps',...
                                      'p0IS_FEM','T11ISamplesDataTimeSteps','T22ISamplesDataTimeSteps','T12ISamplesDataTimeSteps','vxISamplesDataTimeSteps','vyISamplesDataTimeSteps',...
                                      'MeshANodes','MeshAElements','MeshADimensns','PrecomputedIntrplteObjectsA','pinfoA','SensorsA','xA','yA','NpA','KA',...
                                      'errorT11','errorT22','errorT12','errorvx','errorvy',...
                                      'WT11','WT22','WT12','Wvx','Wvy',...
                                      'ST11','ST22','ST12','Svx','Svy');

%=== Load AE Stats ===%
RunOptions.SaveFileNameAESampleStats = sprintf('AEStats-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
load(RunOptions.SaveFileNameAESampleStats,'SampleMeanvx')
load(RunOptions.SaveFileNameAESampleStats,'SampleCovvx')

%=========================================================================%
%                       Plotting Reconstructions
%=========================================================================%
VisualizingServerOutputs

%=========================================================================%
%                           Plotting AE Stats
%=========================================================================%
close all
SelectedSensor = 2;
SampleNumber = 47;

%=== Plotting Diagonal of Covariance vs. Sample Mean vs. Noise Model Variance ===%
CovDiagonal = reshape(diag(SampleCovvx),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
for SelectedSensor = 30:76;
    figure
    plot(1:1:RunOptions.NumberofTimeSteps,SampleMeanvx(SelectedSensor,:),'r');
    hold on
    plot(1:1:RunOptions.NumberofTimeSteps,DataVrblsWave.Evx(SelectedSensor)*ones(RunOptions.NumberofTimeSteps,1),'g')
    plot(1:1:RunOptions.NumberofTimeSteps,sqrt(CovDiagonal(SelectedSensor,:)),'b')
    plot(1:1:800,zeros(800,1),'k')
    ylim([-100,100])
    MeanCovGraphTitle = 'E[\epsilon] vs. sqrt(diag(\Gamma_{\epsilon})) vs. \sigma_E';
    SensorMeanCovGraphTitle = sprintf('Sensor %d: %s',SelectedSensor,MeanCovGraphTitle);
    title(SensorMeanCovGraphTitle,'FontWeight','bold')
    legend('AE Sample Mean','Noise Model','AE Sample Covariance')
    pause(1)
    close all
end
keyboard

%=== Plotting Sample Covariance ===%
figure
imagesc(1:1:RunOptions.NumberofTimeSteps*DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps*DataVrblsWave.NumberofSensors:-1:1,SampleCovvx);
view(2)
colorbar
axis off

%=========================================================================%
%              Plotting Forward Wave Propagation at Sensor
%=========================================================================%
%=== Plotting Sample Wave Propagation at Sensor ===%
figure
vxASamplesDataTimeSteps_SampleN = reshape(vxASamplesDataTimeSteps(:,SampleNumber),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
plot(1:1:RunOptions.NumberofTimeSteps,vxASamplesDataTimeSteps_SampleN(SelectedSensor,:),'r')

%=== Plotting Sample Wave Propagation at Sensor ===%
hold on
vxISamplesDataTimeSteps_SampleN = reshape(vxISamplesDataTimeSteps(:,SampleNumber),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
plot(1:1:RunOptions.NumberofTimeSteps,vxISamplesDataTimeSteps_SampleN(SelectedSensor,:),'b')

SampleWavePropTitle = sprintf('Sample Wave Propagation: Sample %d, Sensor %d - Accurate vs Less Accurate Model',SampleNumber,SelectedSensor);
title(SampleWavePropTitle,'FontWeight','bold')
legend('Accurate','Less Accurate')

%=== Plotting Approximation Error of Sample Wave Propagation at Sensor ===%
figure
errorvx_SamplenN = reshape(errorvx(:,SampleNumber),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
Diff = errorvx_SamplenN(SelectedSensor,:);
plot(1:1:RunOptions.NumberofTimeSteps,Diff)
SampleDiffWavePropTitle = sprintf('Approximation Error of Sample Wave Propagation: Sample %d, Sensor %d',SampleNumber,SelectedSensor);
title(SampleDiffWavePropTitle,'FontWeight','bold')

%=========================================================================%
%             Plotting All Sample Wave Propagation at Sensor
%=========================================================================%
figure
    hold on
for n=1:30%RunOptions.AEN_Samples
    %=== Accurate Model ===%
    vxASamplesDataTimeSteps_SampleN = reshape(vxASamplesDataTimeSteps(:,n),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    plot(1:1:RunOptions.NumberofTimeSteps,vxASamplesDataTimeSteps_SampleN(SelectedSensor,:),'r')
    
    %=== Less Accurate Model ===%
    vxISamplesDataTimeSteps_SampleN = reshape(vxISamplesDataTimeSteps(:,n),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    plot(1:1:RunOptions.NumberofTimeSteps,vxISamplesDataTimeSteps_SampleN(SelectedSensor,:),'b')
    
    SampleWavePropTitle = sprintf('Sample Wave Propagation: All Samples, Sensor %d - Accurate vs Less Accurate Model',SelectedSensor);
    title(SampleWavePropTitle,'FontWeight','bold')
    legend('Accurate','Less Accurate')
end

%=========================================================================%
%                Plotting Approximation Error at Sensor
%=========================================================================%
errorvx_Sample = reshape(errorvx(:,SampleNumber),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
figure
plot(1:1:RunOptions.NumberofTimeSteps,errorvx_Sample(SelectedSensor,:))

ApproxErrorPropTitle = sprintf('Approximation Error at selected sensor %d',SelectedSensor);
title(ApproxErrorPropTitle,'FontWeight','bold')

%=========================================================================%
%              Plotting All Approximation Errors at Sensor
%=========================================================================%
close all
SelectedSensor = 58;
figure
hold on
for n=1:30%RunOptions.AEN_Samples
    errorvx_SampleN = reshape(errorvx(:,n),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    plot(1:1:RunOptions.NumberofTimeSteps,errorvx_SampleN(SelectedSensor,:),'r')
    plot(1:1:800,zeros(800,1),'k')    
    AllEpsilonTitle = sprintf('Approximation Error: All Samples, Sensor %d',SelectedSensor);
    title(AllEpsilonTitle,'FontWeight','bold')
end

%=========================================================================%
%              Plotting Approximation Error Energy
%=========================================================================%
AENormExpSquared = zeros(RunOptions.NumberofTimeSteps,1);
AETraceCovTime = zeros(RunOptions.NumberofTimeSteps,1);
AEEnergy = zeros(RunOptions.NumberofTimeSteps,1);
for t=1:RunOptions.NumberofTimeSteps
    AENormExpSquared(t) = SampleMeanvx(:,t)'*SampleMeanvx(:,t);
    AETraceCovTime(t) = sum(CovDiagonal(:,t));
    AEEnergy(t) = AENormExpSquared(t) + AETraceCovTime(t);
end
NoiseNormExpSquared = zeros(RunOptions.NumberofTimeSteps,1);
NoiseTraceCovTime = zeros(RunOptions.NumberofTimeSteps,1);
NoiseEnergy = zeros(RunOptions.NumberofTimeSteps,1);
ReshapeEvx = reshape(DataVrblsWave.Evx,DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
for t=1:RunOptions.NumberofTimeSteps
    NoiseNormExpSquared(t) = 0;
    NoiseTraceCovTime(t) = sum(ReshapeEvx(:,t).^2);
    NoiseEnergy(t) = NoiseNormExpSquared(t) + NoiseTraceCovTime(t);
end

figure
hold on
plot(1:1:RunOptions.NumberofTimeSteps,AEEnergy,'r')
plot(1:1:RunOptions.NumberofTimeSteps,NoiseEnergy,'g')

%=========================================================================%
%                       Plotting Singular Values
%=========================================================================%
figure
plot(1:1:size(Svx,2),diag(Svx).^2)
title('Singular Values of Wvx','FontWeight','bold')

%=========================================================================%
%                   Simulating Forward Wave Propagation
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%
%%% Accurate Model %%%
%%%%%%%%%%%%%%%%%%%%%%
%=== Plotting Samples ===%
figure
FemPlot2D(MeshANodes,MeshAElements,p0AS_FEM(:,SampleNumber));

%=== Plotting Sample Wave Propagation ===%
PLOT.TRI_DGMMeshD = delaunay(xA,yA);
[p0AS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),p0AS_FEM(:,SampleNumber),xA,yA,NpA*KA,PrecomputedIntrplteObjectsA);
IniCond = p0AS_DGM/(RunOptions.ElasticMediumDensity*RunOptions.SpecificHeatCoeff);
[T11ASTimeSteps,T22ASTimeSteps,T12ASTimeSteps,vxASTimeSteps,vyASTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,xA,yA,NpA,KA,pinfoA,dt,PLOT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Less Accurate Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Plotting Samples ===%
figure
FemPlot2D(MeshI.Nodes,MeshI.Elements,p0IS_FEM(:,SampleNumber));

%=== Plotting Sample Wave Propagation ===%
DGMMeshI.pinfo = EWE_DGM2D_PrecomputeUpwindFluxPNonConf(RunOptions,DGMMeshI.pinfo,DGMMeshI.Norder,DGMMeshI.rho,DGMMeshI.lambda,DGMMeshI.mu,MeshI.DomainIndices,RunOptions.FluidDomainWithSolidLayerMeshI,RunOptions.SolidMeshI);  
PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshI;
[p0IS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshI.Elements,1),MeshI.Nodes*(1/RunOptions.ScalingOptical),p0IS_FEM(:,SampleNumber),DGMMeshI.x,DGMMeshI.y,DGMMeshI.Np*DGMMeshI.K,PrecomputedIntrplteObjectsI);
IniCond = p0IS_DGM/(RunOptions.ElasticMediumDensity*RunOptions.SpecificHeatCoeff);
[T11ISTimeSteps,T22ISTimeSteps,T12ISTimeSteps,vxISTimeSteps,vyISTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,DGMMeshI.x,DGMMeshI.y,DGMMeshI.Np,DGMMeshI.K,DGMMeshI.pinfo,dt,PLOT);

%=========================================================================%
%                Forward Wave Propagation on True p0
%=========================================================================%
PLOT.DGMPlotBirdsEyeView = 1;
figure
FemPlot2D(MeshD.Nodes,MeshD.Elements,full(DataVrblsWave.p0));
%=== Accurate Wave Propagation ===%
[p0TrueAS_FEM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshD.Elements,1),MeshD.Nodes*(1/RunOptions.ScalingOptical),DataVrblsWave.p0,MeshANodes(:,1),MeshANodes(:,2),size(MeshANodes,1),0);
PLOT.TRI_DGMMeshD = delaunay(xA,yA);
[p0TrueAS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),p0TrueAS_FEM,xA,yA,NpA*KA,0);
IniCond = p0TrueAS_DGM/(RunOptions.ElasticMediumDensity*RunOptions.SpecificHeatCoeff);
[T11TrueASTimeSteps,T22TrueASTimeSteps,T12TrueASTimeSteps,vxTrueASTimeSteps,vyTrueASTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,xA,yA,NpA,KA,pinfoA,dt,PLOT);
T11TrueDataASTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
T22TrueDataASTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
T12TrueDataASTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
vxTrueDataASTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
vyTrueDataASTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
%=== Interpolation to Form Sensory Data ===%
for t=1:RunOptions.NumberofTimeSteps
    for s=1:DataVrblsWave.NumberofSensors
        if RunOptions.FullqVectorData == 1;
            T11TrueDataASTimeSteps(s,t) = SensorsA{s}.l_iatsensor*T11TrueASTimeSteps(SensorsA{s}.id,t);
            T22TrueDataASTimeSteps(s,t) = SensorsA{s}.l_iatsensor*T22TrueASTimeSteps(SensorsA{s}.id,t);
            T12TrueDataASTimeSteps(s,t) = SensorsA{s}.l_iatsensor*T12TrueASTimeSteps(SensorsA{s}.id,t);
        end
        vxTrueDataASTimeSteps(s,t) = SensorsA{s}.l_iatsensor*vxTrueASTimeSteps(SensorsA{s}.id,t);
        vyDTrueDataASTimeSteps(s,t) = SensorsA{s}.l_iatsensor*vyTrueASTimeSteps(SensorsA{s}.id,t);
    end
end
clear T11TrueASTimeSteps T22TrueASTimeSteps T12TrueASTimeSteps vxTrueASTimeSteps vyTrueASTimeSteps

%=== Less Accurate Wave Propagation ===%
[p0TrueIS_FEM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshD.Elements,1),MeshD.Nodes*(1/RunOptions.ScalingOptical),DataVrblsWave.p0,MeshI.Nodes(:,1),MeshI.Nodes(:,2),size(MeshI.Nodes,1),0);
PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshI;
DGMMeshI.pinfo = EWE_DGM2D_PrecomputeUpwindFluxPNonConf(RunOptions,DGMMeshI.pinfo,DGMMeshI.Norder,DGMMeshI.rho,DGMMeshI.lambda,DGMMeshI.mu,MeshI.DomainIndices,RunOptions.FluidDomainWithSolidLayerMeshI,RunOptions.SolidMeshI);
[p0TrueIS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshI.Elements,1),MeshI.Nodes*(1/RunOptions.ScalingOptical),p0TrueIS_FEM,DGMMeshI.x,DGMMeshI.y,DGMMeshI.Np*DGMMeshI.K,PrecomputedIntrplteObjectsI);
IniCond = p0TrueIS_DGM/(RunOptions.ElasticMediumDensity*RunOptions.SpecificHeatCoeff);
[T11TrueISTimeSteps,T22TrueISTimeSteps,T12TrueISTimeSteps,vxTrueISTimeSteps,vyTrueISTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,DGMMeshI.x,DGMMeshI.y,DGMMeshI.Np,DGMMeshI.K,DGMMeshI.pinfo,dt,PLOT);
T11TrueDataISTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
T22TrueDataISTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
T12TrueDataISTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
vxTrueDataISTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
vyTrueDataISTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
%=== Interpolation to Form Sensory Data ===%
for t=1:RunOptions.NumberofTimeSteps
    for s=1:DataVrblsWave.NumberofSensors
        if RunOptions.FullqVectorData == 1;
            T11TrueDataISTimeSteps(s,t) = DataVrblsWave.SensorsI{s}.l_iatsensor*T11TrueISTimeSteps(DataVrblsWave.SensorsI{s}.id,t);
            T22TrueDataISTimeSteps(s,t) = DataVrblsWave.SensorsI{s}.l_iatsensor*T22TrueISTimeSteps(DataVrblsWave.SensorsI{s}.id,t);
            T12TrueDataiSTimeSteps(s,t) = DataVrblsWave.SensorsI{s}.l_iatsensor*T12TrueISTimeSteps(DataVrblsWave.SensorsI{s}.id,t);
        end
        vxTrueDataISTimeSteps(s,t) = DataVrblsWave.SensorsI{s}.l_iatsensor*vxTrueISTimeSteps(DataVrblsWave.SensorsI{s}.id,t);
        vyTrueDataISTimeSteps(s,t) = DataVrblsWave.SensorsI{s}.l_iatsensor*vyTrueISTimeSteps(DataVrblsWave.SensorsI{s}.id,t);
    end
end
clear T11TrueISTimeSteps T22TrueISTimeSteps T12TrueISTimeSteps vxTrueISTimeSteps vyTrueISTimeSteps

%=== Noise Model ===%
RunOptions.Cov_ENoiseLevel = 0.001;
%Max of Data
if RunOptions.FullqVectorData == 1;
    MaxT11 = max(DataVrblsWave.T11DataTimeSteps(:));
    MaxT22 = max(DataVrblsWave.T22DataTimeSteps(:));
    MaxT12 = max(DataVrblsWave.T12DataTimeSteps(:));
end
Maxvx = max(DataVrblsWave.vxDataTimeSteps(:));
Maxvy = max(DataVrblsWave.vyDataTimeSteps(:));
NoiseSTD = Maxvx*RunOptions.Cov_ENoiseLevel*ones(DataVrblsWave.NumberofSensors);
%Min minus max of Data
DataVrblsWave.MaxMinusMinT11 = zeros(1,DataVrblsWave.NumberofSensors);
DataVrblsWave.MaxMinusMinT22 = zeros(1,DataVrblsWave.NumberofSensors);
DataVrblsWave.MaxMinusMinT12 = zeros(1,DataVrblsWave.NumberofSensors);
DataVrblsWave.MaxMinusMinvx = zeros(1,DataVrblsWave.NumberofSensors);
DataVrblsWave.MaxMinusMinvy = zeros(1,DataVrblsWave.NumberofSensors);
for s = 1:DataVrblsWave.NumberofSensors
    if RunOptions.FullqVectorData == 1;
        DataVrblsWave.MaxMinusMinT11(s) = max(DataVrblsWave.T11DataTimeSteps(s,:)) - min(DataVrblsWave.T11DataTimeSteps(s,:));
        DataVrblsWave.MaxMinusMinT22(s) = max(DataVrblsWave.T22DataTimeSteps(s,:)) - min(DataVrblsWave.T22DataTimeSteps(s,:));
        DataVrblsWave.MaxMinusMinT12(s) = max(DataVrblsWave.T12DataTimeSteps(s,:)) - min(DataVrblsWave.T12DataTimeSteps(s,:));
    end
    DataVrblsWave.MaxMinusMinvx(s) = max(DataVrblsWave.vxDataTimeSteps(s,:)) - min(DataVrblsWave.vxDataTimeSteps(s,:));
    DataVrblsWave.MaxMinusMinvy(s) = max(DataVrblsWave.vyDataTimeSteps(s,:)) - min(DataVrblsWave.vyDataTimeSteps(s,:));
end
NoiseSTD = RunOptions.Cov_ENoiseLevel*DataVrblsWave.MaxMinusMinvx;

%=== Plotting ===%
for SelectedSensor = 30:DataVrblsWave.NumberofSensors
    figure
    set(gcf, 'Position', get(0, 'Screensize'));
    plot(1:1:RunOptions.NumberofTimeSteps,DataVrblsWave.vxDataTimeSteps(SelectedSensor,:),'r')
    hold on
    plot(1:1:RunOptions.NumberofTimeSteps,vxTrueDataASTimeSteps(SelectedSensor,:),'b')
    plot(1:1:RunOptions.NumberofTimeSteps,vxTrueDataISTimeSteps(SelectedSensor,:),'g')
    plot(1:1:RunOptions.NumberofTimeSteps,DataVrblsWave.vxDataTimeSteps(SelectedSensor,:)+2*sqrt(NoiseSTD(SelectedSensor)),'k--')
    plot(1:1:RunOptions.NumberofTimeSteps,DataVrblsWave.vxDataTimeSteps(SelectedSensor,:)-2*sqrt(NoiseSTD(SelectedSensor)),'k--')
    SensorDataTitle = sprintf('Sensor Data of All Forward Operators at Sensor %d',SelectedSensor);
    title(SensorDataTitle,'FontWeight','bold')
%     ylim([-400 400])
    legend('Noiseless Data','Accurate','Less Accurate')
    keyboard
    close all
end


