close all
clc

NumberofSensors = size(SensorsI,2);
NumberofTimeSteps = RunOptions.NumberofTimeSteps;
PLOT.PriorSamples = 1;
PLOT.Figure_WavePropagation_Data = figure;
PLOT.Figure_WavePropagation_Data_Title = 'Elastic Wave Propagation';
movegui(PLOT.Figure_WavePropagation_Data,'southwest');
drawnow

%=========================================================================%
%                           Accurate Sample
%=========================================================================%
%=== Generating Sample ===%
[~,~,~,p0AS_FEM] = SmoothnessPrior_AutoCorr(MeshANodes,Prior.Exp_p0,Prior.AC_Var_p0,0.001,1,PLOT);

%=== Plotting Interpolated Sample ===%
[p0AS_DGM,~] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),p0AS_FEM,xA,yA,NpA*KA,PrecomputedIntrplteObjectsA);
figure
trisurf(PLOT.TRI_DGMMeshA,xA,yA,p0AS_DGM);
if PLOT.DGMPlotBirdsEyeView == 1;
    view(2)
end
shading interp %thanks Ru!
zlim([0,800]);
caxis([0,800]);
colorbar
colormap(jet(256))
title('Interpolated Accurate Sample','FontWeight','bold')

%=========================================================================%
%                           Inaccurate Sample
%=========================================================================%
%=== Plotting Interpolated Sample ===%
figure
[p0IS_FEM, PrecomputedIntrplteObjectsAtoIFEM] = IntrplteOver2DTriangulatedMesh(size(MeshAElements,1),MeshANodes*(1/RunOptions.ScalingOptical),p0AS_FEM(:,1),MeshINodes(:,1),MeshINodes(:,2),size(MeshINodes,1),0);
[p0IS_DGM, ~] = IntrplteOver2DTriangulatedMesh(size(MeshIElements,1),MeshINodes*(1/RunOptions.ScalingOptical),p0IS_FEM,xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
trisurf(PLOT.TRI_DGMMeshI,xI,yI,p0IS_DGM);
if PLOT.DGMPlotBirdsEyeView == 1;
    view(2)
end
shading interp %thanks Ru!
zlim([0,300]);zlim([0,300]);
caxis([0,300]);
colorbar
colormap(jet(256))
title('Interpolated Inaccurate Sample','FontWeight','bold')

%=========================================================================%
%                           Forward Wave Propagation
%=========================================================================%
%=== Accurate Model ===%
PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshA;
IniCond = p0AS_DGM/(RunOptions.ElasticMediumDensity*RunOptions.SpecificHeatCoeff);
[T11ASTimeSteps,T22ASTimeSteps,T12ASTimeSteps,vxASTimeSteps,vyASTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,xA,yA,NpA,KA,pinfoA,dt,PLOT);
T11ASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
T22ASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
T12ASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
vxASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
vyASDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
%=== Interpolation to Form Sensory Data ===%
for t=1:NumberofTimeSteps
    for s=1:NumberofSensors
        vxASDataTimeSteps(s,t) = SensorsA{s}.l_iatsensor*vxASTimeSteps(SensorsA{s}.id,t);
        vyASDataTimeSteps(s,t) = SensorsA{s}.l_iatsensor*vyASTimeSteps(SensorsA{s}.id,t);
    end
end

%=== Inaccurate Model ===%
PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshI;
IniCond = p0IS_DGM/(RunOptions.ElasticMediumDensity*RunOptions.SpecificHeatCoeff);
[T11ISTimeSteps,T22ISTimeSteps,T12ISTimeSteps,vxISTimeSteps,vyISTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT);
T11ISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
T22ISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
T12ISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
vxISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
vyISDataTimeSteps = sparse(NumberofSensors,NumberofTimeSteps);
%=== Interpolation to Form Sensory Data ===%
for t=1:NumberofTimeSteps
    for s=1:NumberofSensors
        vxISDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxISTimeSteps(SensorsI{s}.id,t);
        vyISDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyISTimeSteps(SensorsI{s}.id,t);
    end
end

%=========================================================================%
%                           Plotting Data at Sensor
%=========================================================================%
close all
SelectedSensor = 57;
figure
plot(dt:dt:RunOptions.NumberofTimeSteps*dt,sqrt(vxASDataTimeSteps(SelectedSensor,:).^2 + vyASDataTimeSteps(SelectedSensor,:).^2),'r')
ylim([0,max(sqrt(full(vxASDataTimeSteps(SelectedSensor,:)).^2 + full(vyASDataTimeSteps(SelectedSensor,:)).^2))]);
% WavePropAtSensorTitle = sprintf('Accurate Sample Wave Propagation at Sensor %d',SelectedSensor);
% title(WavePropAtSensorTitle,'FontWeight','bold')

figure
plot(dt:dt:RunOptions.NumberofTimeSteps*dt,sqrt(vxISDataTimeSteps(SelectedSensor,:).^2 + vyISDataTimeSteps(SelectedSensor,:).^2),'r')
ylim([0,max(sqrt(full(vxASDataTimeSteps(SelectedSensor,:)).^2 + full(vyASDataTimeSteps(SelectedSensor,:)).^2))]);
% WavePropAtSensorTitle = sprintf('Inaccurate Sample Wave Propagation at Sensor %d',SelectedSensor);
% title(WavePropAtSensorTitle,'FontWeight','bold')

figure
plot(dt:dt:RunOptions.NumberofTimeSteps*dt,sqrt(full(vxASDataTimeSteps(SelectedSensor,:)).^2 + full(vyASDataTimeSteps(SelectedSensor,:)).^2) - sqrt(full(vxISDataTimeSteps(SelectedSensor,:)).^2 + full(vyISDataTimeSteps(SelectedSensor,:)).^2),'r')
ylim([-200,400]);
% WavePropAtSensorTitle = sprintf('Error of Sample Wave Propagation at Sensor %d',SelectedSensor);
% title(WavePropAtSensorTitle,'FontWeight','bold')