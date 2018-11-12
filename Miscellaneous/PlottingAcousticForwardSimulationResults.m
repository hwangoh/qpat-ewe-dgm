close all
clear all
clc

addpath('C:\Users\hgoh009\Dropbox\Work\PhD\PhD_Stuff\Codes\QPAT_EWE_DGM2D_Codes_Hwan')

%% %%%%%%%%%%%%%%%%%%%
%%% Save File Name %%%
%%%%%%%%%%%%%%%%%%%%%%
%=== Parameter Type ===%
RunOptions.TestAcousticParameters = 0; %Use test acoustic parameters
RunOptions.RealAcousticParameters = 1; %Use realistic acoustic parameters
%=== Acoustic Forward Problem ===%
RunOptions.FluidDomainWithSolidLayerMeshD = 1; %Use Fluid domain with solid layer representing the skull for data mesh
RunOptions.FluidDomainWithSolidLayerMeshI = 0; %Use Fluid domain with solid layer representing the skull for inverse mesh
RunOptions.FluidMeshD = 0; %Use purely fluid domain for data mesh
RunOptions.FluidMeshI = 1; %Use purely fluid domain for inverse mesh
RunOptions.SolidMeshD = 0; %Use purely solid domain for data mesh
RunOptions.SolidMeshI = 0; %Use purely solid domain for inverse mesh
RunOptions.VelocitiesData = 1; %Output of forward problem sensory data of velocities, only works when sensors are placed on the boundary
RunOptions.NumberofSensorsOnOneBoundaryEdge = 20; %Number of sensors on one boundary edge of the domain
if RunOptions.TestAcousticParameters == 1;
    RunOptions.FinalTime = 0.004; %Final Time
end
if RunOptions.RealAcousticParameters == 1;
    RunOptions.FinalTime = 0.000004; %Final Time
end

%=== Adding noise ===%
RunOptions.AddNoise = 0; %Add noise?
RunOptions.NoiseLevel = 0.05; %Scalar multiplier of noise draws
RunOptions.NoiseMax = 1; %Use max minus min of data
RunOptions.NoiseMinMax = 0; %Use max minus min of data

%=== Filenames of Outputs ===%
if RunOptions.TestAcousticParameters == 1;
    RunOptions.SaveFileNameParameterType = 'TestPrmtrs';
    TempString = num2str(RunOptions.FinalTime);
    RunOptions.SaveFileNameFinalTime = sscanf(TempString(3:end),'%s');
end
if RunOptions.RealAcousticParameters == 1;
    RunOptions.SaveFileNameParameterType = 'RealPrmtrs';
    TempString = num2str(RunOptions.FinalTime);
    RunOptions.SaveFileNameFinalTime = sscanf(TempString(1),'%s');
end
if RunOptions.FluidDomainWithSolidLayerMeshD == 1 %For saving data
    RunOptions.SaveFileNameDataDomain = 'SLDataDomain';
end
if RunOptions.FluidMeshD == 1 %For saving data
    RunOptions.SaveFileNameDataDomain = 'FDataDomain';
end
if RunOptions.AddNoise == 0;
    RunOptions.NoiseLevel = 0;
end
TempString = num2str(RunOptions.NoiseLevel);
RunOptions.SaveFileNameNoiseLevel = sscanf(TempString(3:end),'%s');
if RunOptions.NoiseMax == 1
    RunOptions.SaveFileNameNoiseType = 'Max';
end
if RunOptions.NoiseMinMax == 1
    RunOptions.SaveFileNameNoiseType = 'MinMax';
end
clear TempString

%% %%%%%%%%%%%%%%%%%%%%
%%% Loading Outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%
MeshSizes = {'002', '001', '0009', '0007', '00055', '00038', '00034', '00032', '0003'};
for ii=1:length(MeshSizes)
    TrelisMeshDElementSize = MeshSizes{ii};
    RunOptions.SaveFileNameData = sprintf('%s-%s-%sNoise%s-%sD-%dSensors-%sFinalTime',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDataDomain,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,TrelisMeshDElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime);
    load(RunOptions.SaveFileNameData,'vxDataTimeSteps','vyDataTimeSteps')
    eval([sprintf('Data%s',TrelisMeshDElementSize) '= sqrt(vxDataTimeSteps.^2 + vyDataTimeSteps.^2);']);
end

%% %%%%%%%%%%%
%%% Errors %%%
%%%%%%%%%%%%%%
OutputErrors = zeros(length(MeshSizes)-1,1);
for ii=1:length(MeshSizes)-1
    TrelisMeshDElementSize = MeshSizes{ii};
    LessAccData = eval(sprintf('Data%s',TrelisMeshDElementSize));
    VectorizedError = Data0003 - LessAccData;
    VectorizedError = VectorizedError(:);
    OutputErrors(ii) = norm(VectorizedError,2);
end
figure
Meshes = ['8', '7', '6', '5', '4', '3', '2', '1'];
plot(OutputErrors,'-o');
set(gca,'xticklabel',Meshes.')
ylim([0 60]);

%% %%%%%%%%%
%%% Data %%%
%%%%%%%%%%%%
SelectedSensor = 16;
if RunOptions.TestAcousticParameters == 1;
    dt = 1.14248568461526e-05;
    NumberofTimeSteps = 351;
end
if RunOptions.RealAcousticParameters == 1;
    dt = 5.385729166666684e-09;
    NumberofTimeSteps = 743;
end
figure
plot(dt:dt:NumberofTimeSteps*dt,Data0003(SelectedSensor,:),'b')
ylim([0 2]);


