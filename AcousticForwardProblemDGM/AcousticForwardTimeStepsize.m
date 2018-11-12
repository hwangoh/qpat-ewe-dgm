%%%%%%%%%%%%%%%%%%%%%%
%%% Time Step Size %%%
%%%%%%%%%%%%%%%%%%%%%%
%=== Timo's CFL Condition ===%
if RunOptions.TimeStepUseCFLTimo ~= 0;
    cp = max(RunOptions.AcousticMediumWaveSpeed_cp,RunOptions.ElasticMediumWaveSpeed_cp);
    cs = RunOptions.ElasticMediumWaveSpeed_cs;
    maxwavevelocity = max(max(cp,cs)); %1xK row vector with each entry representing the max wave speed in that element
    if exist('hmaxhmin','var')
        hmax = hmaxhmin(:,1);
        hmin = hmaxhmin(:,1);
    end
    if exist('MeshD','var')
        hmin = abs(min(DGMMeshD.x(1) - DGMMeshD.x(2), DGMMeshD.y(DGMMeshD.Np) - DGMMeshD.y(DGMMeshD.Np-1)));
    end
    dt = min(hmin*RunOptions.TimeStepUseCFLTimo./(DGMMeshD.Norder.^2.*maxwavevelocity(:)));
end
%=== Set Own Stepsize ===%
if RunOptions.TimeStepUseCFLTimo == 0;
    dt = RunOptions.TimeStepSizeLSERK4;
end
% dt = 1.14248568461526e-05; %dt of 0003 Trelis mesh with test parameters
% dt = 5.385729166666684e-09; %dt of 0003 Trelis mesh with real parameters

%=== Number of Time Steps ===%
RunOptions.NumberofTimeSteps = ceil(RunOptions.FinalTime/dt); %Number of time steps