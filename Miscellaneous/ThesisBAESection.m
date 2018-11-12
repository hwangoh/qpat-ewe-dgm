% Instructions:
% 1) Load AEStats
% 2) CovDiagonal = diag(SampleCovvx);
% 3) CovDiagonal = reshape(CovDiagonal,36,743);
% 4) save('CovDiagonal','CovDiagonal','-v7.3')
% 5) save('SampleMeanvx','SampleMeanvx','-v7.3')
% 6) Load AESamples
% 7) save('errorvx','errorvx','-v7.3')
% 8) Load DataVrblsWave.Evx

close all
SelectedSensor = 30;
% dt = 1.14248568461526e-05; %dt of 0003 Trelis mesh with test parameters
dt = 5.385729166666684e-09; %dt of 0003 Trelis mesh with real parameters

%=========================================================================%
%                           Plotting AE Stats
%=========================================================================%
figure
plot(dt:dt:RunOptions.NumberofTimeSteps*dt,SampleMeanvx(SelectedSensor,:),'r');
hold on
plot(dt:dt:RunOptions.NumberofTimeSteps*dt,DataVrblsWave.Evx(SelectedSensor)*ones(RunOptions.NumberofTimeSteps,1),'g')
plot(dt:dt:RunOptions.NumberofTimeSteps*dt,sqrt(CovDiagonal(SelectedSensor,:)),'b')
plot(dt:dt:800*dt,zeros(800,1),'k')
ylim([-200,200])
xlim([0,800*dt])

%=========================================================================%
%         Plotting Approximation Errors at Sensor for 30 Samples
%=========================================================================%
figure
hold on
for n=1:30 %RunOptions.AEN_Samples
    errorvx_SampleN = reshape(errorvx(:,n),DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    plot(dt:dt:RunOptions.NumberofTimeSteps*dt,errorvx_SampleN(SelectedSensor,:),'r')
    plot(dt:dt:800*dt,zeros(800,1),'k')  
    ylim([-800,800])
    xlim([0,800*dt])
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
plot(dt:dt:RunOptions.NumberofTimeSteps*dt,AEEnergy,'r')
plot(dt:dt:RunOptions.NumberofTimeSteps*dt,NoiseEnergy,'g')