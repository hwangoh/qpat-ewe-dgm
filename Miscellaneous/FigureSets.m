close all
clear all
clc
restoredefaultpath


% addpath(genpath('D:\Dropbox\Work\PhD\PhD_Stuff\Codes\Outputs'))
% addpath(genpath('C:\Users\hgoh009\Dropbox\Work\PhD\PhD_Stuff\Codes\Outputs'))
% addpath(genpath('C:\Users\TCB\Dropbox\Work\PhD\PhD_Stuff\Codes\Outputs'))
% addpath(genpath('C:\Users\Hwan\Downloads'))
addpath('C:\Users\hgoh009\Downloads')
% addpath(genpath('C:\Users\TCB\Downloads'))

%% %%%%%%%%%%%%%%%%%%%
%%% Save File Name %%%
%%%%%%%%%%%%%%%%%%%%%%

%=== Mesh Properties ===%
RunOptions.TrelisMeshDElementSize = '0003'; %Entry for generating Trelis data mesh, main purpose is for the file name when saving
RunOptions.TrelisMeshIElementSize = '0007'; %Entry for generating Trelis inverse mesh, main purpose is for the file name when saving
RunOptions.TrelisMeshAElementSize = '00032'; %Entry for generating Trelis accurate mesh, main purpose is for the file name when saving

%=== Acoustic Forward Problem ===%
RunOptions.TestAcousticParameters = 0; %Use test acoustic parameters
RunOptions.RealAcousticParameters = 1; %Use realistic acoustic parameters
RunOptions.NumberofSensorsOnOneBoundaryEdge = 10; %Number of sensors on one boundary edge of the domain
if RunOptions.TestAcousticParameters == 1;
    RunOptions.FinalTime = 0.008; %Final Time
end
if RunOptions.RealAcousticParameters == 1;
    RunOptions.FinalTime = 0.000004; %Final Time
end

%=== Noise Properties ===%
RunOptions.AddNoise = 1; %Add noise?
RunOptions.NoiseMinMax = 0; %Use max minus min of data at each sensor
RunOptions.NoiseMinMaxS = 0; %Use max minus min of data at each sensor
RunOptions.NoiseMax = 1; %Use max of all data data
RunOptions.NoiseMaxS = 0; %Use max of data at each sensor 

%=== Inverse Problem Properties ===%
RunOptions.UseACPrior = 0;
RunOptions.AE_EWE = 0; %Noise model accounts for approximation error due to discretization
RunOptions.AEN_Samples = 1000; 
RunOptions.EWE_LS_LogExpTransform = 1; %Use log-exp transform positivity constraint
RunOptions.EWE_LS_SigmoidTransform = 0; %Use sigmoidal transform positivity constraint
RunOptions.EWE_LS_ExponentialTransform = 0; %Use exponential transform positivity constraint

%% %%%%%%%%%%%%%%%%%%%
%%% Figure Outputs %%%
%%%%%%%%%%%%%%%%%%%%%%
RelativeErrors = zeros(3,3);
NoiseLevels = [0.05, 0.01, 0.001];
% NoiseLevels = [0.3, 0.2, 0.1];
Cases = {[0 0; 1 1],[1 1; 0 0],[1 0; 0 1]};
count = 1;
for ii=1:length(NoiseLevels)
    RunOptions.NoiseLevel = NoiseLevels(ii);
    RunOptions.Cov_ENoiseLevel = RunOptions.NoiseLevel;
    for jj=1:length(Cases)
        RunOptions.FluidDomainWithSolidLayerMeshD = Cases{jj}(1,1); %Use Fluid domain with solid layer representing the skull for data mesh
        RunOptions.FluidDomainWithSolidLayerMeshI = Cases{jj}(1,2); %Use Fluid domain with solid layer representing the skull for inverse mesh
        RunOptions.FluidMeshD = Cases{jj}(2,1); %Use purely fluid domain for data mesh
        RunOptions.FluidMeshI = Cases{jj}(2,2); %Use purely fluid domain for inverse mesh
        FilenamesofRunOptions
        SaveFileNameReconstructions = sprintf('Reconstructions-%s',RunOptions.SaveFileName)
        try
            load(SaveFileNameReconstructions,'DataVrblsWave','Prior','p0TrueIntrplte','AcousticInverseItrtnInfo','PLOT','MeshINodes','Trmntn')
            p0Recon = AcousticInverseItrtnInfo.p0Recon;
            if RunOptions.NoiseLevel == NoiseLevels(1)
                RelativeErrors(ii,jj) = AcousticInverseItrtnInfo.p0RelativeError(Trmntn,1);
                subplot(3,3,count)
                trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),full(AcousticInverseItrtnInfo.p0Recon(:,Trmntn)));
            end
            if RunOptions.NoiseLevel == NoiseLevels(2)
                RelativeErrors(ii,jj) = AcousticInverseItrtnInfo.p0RelativeError(Trmntn,1);
                subplot(3,3,count)
                trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),full(AcousticInverseItrtnInfo.p0Recon(:,Trmntn)));
            end
            if RunOptions.NoiseLevel == NoiseLevels(3)
                RelativeErrors(ii,jj) = AcousticInverseItrtnInfo.p0RelativeError(Trmntn,1);
                subplot(3,3,count)
                trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),full(AcousticInverseItrtnInfo.p0Recon(:,Trmntn)));
            end
            view(2)
            zlim(PLOT.InitialPressurezAxis);
            shading interp
            caxis([0 200])
            colormap(jet(256))
        end
        count = count + 1
    end
end













