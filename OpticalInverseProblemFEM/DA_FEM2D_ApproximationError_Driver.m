%=== File Names for Saving and Loading ===%
RunOptions.SaveFileNameAESampleStats = sprintf('AEStatsOpt-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%sCorr-%dAESamples-%sAE',RunOptions.SaveFileNameDomain,RunOptions.SaveFileNameStateVector,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.SaveFileNamePriorCorr,RunOptions.AE_DA_N_Samples,RunOptions.TrelisMeshAElementSize);
RunOptions.SaveFileNameAEErrorModel = sprintf('AEErrorModelOpt-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%s%s-%sCorr-%dAESamples-%sAE',RunOptions.SaveFileNameDomain,RunOptions.SaveFileNameStateVector,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.SaveFileNameCov_ENoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.SaveFileNamePriorCorr,RunOptions.AE_DA_N_Samples,RunOptions.TrelisMeshAElementSize);
RunOptions.SaveFileNameAESampleStatsAcoustic = sprintf('AEStats-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%sCorr-%dAESamples-%sAE',RunOptions.SaveFileNameDomain,RunOptions.SaveFileNameStateVector,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.SaveFileNamePriorCorr,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
RunOptions.SaveFileNameAEErrorModelAcoustic = sprintf('AEErrorModel-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%s%s-%sCorr-%dAESamples-%sAE',RunOptions.SaveFileNameDomain,RunOptions.SaveFileNameStateVector,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.SaveFileNameCov_ENoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.SaveFileNamePriorCorr,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);

%% =======================================================================%
%                    Computing AE Sample Statistics
%=========================================================================%
if RunOptions.GenerateAndSave_AESamples == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading Acoustic AE Stats %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf('Loading Acoustic AE enhanced error model statistics');
load(RunOptions.SaveFileNameAEErrorModelAcoustic,'L_ET11','L_ET22','L_ET12','L_Evx','L_Evy')
load(RunOptions.SaveFileNameAESampleStatsAcoustic,'SampleMeanT11','SampleMeanT22','SampleMeanT12','SampleMeanvx','SampleMeanvy')
T11ErrorMean = SampleMeanT11;
T22ErrorMean = SampleMeanT22;
T12ErrorMean = SampleMeanT12;
vxErrorMean = SampleMeanvx;
vyErrorMean = SampleMeanvy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Constructing Sample Statistics %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SampleMeanOpt,SampleCovOpt] = DA_FEM2D_AE_ConstructStats_SVD(RunOptions,DataVrblsOptical,MeshI,PrmtrsI,PrmtrsPrp,Prior,...
                                                              DGMMeshI,PrecomputedIntrplteObjectsI,...
                                                              L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                              T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                                              SensorsI,dt,PLOT);
end
clear L_ET11 L_ET22 L_ET12 L_Evx L_Evy T11ErrorMean T22ErrorMean T12ErrorMean vxErrorMean vyErrorMean DGMMeshI SensorsI

%=== Saving and Loading AE Sample Stats ===%
if RunOptions.Save_AEStats == 1
    printf('Saving AE sample statistics:');
    save(RunOptions.SaveFileNameAESampleStats,'SampleMeanOpt','SampleCovOpt','-v7.3')
end 
if RunOptions.Load_AEStats == 1
    printf('Loading AE sample statistics:');
    load(RunOptions.SaveFileNameAESampleStats,'SampleMeanOpt','SampleCovOpt')
end

%% =======================================================================%
%                           Computing Noise Model
%=========================================================================%
%=== Enhanced Error Model Mean ===%
ErrorMeanOpt = ErrorMeanOpt + SampleMeanOpt;

%=== Enhanced Error Model Covariance Matrix Square Root ===%
printf('Computing Cholesky decomposition of sample covariances');
if RunOptions.Save_AEErrorModel == 1
    tic
    L_E = inv(chol(diag(DataVrblsOptical.E) + SampleCovOpt));
    toc    
    %=== Saving and Loading AE Enhanced Error Model Statistics ===%
    printf('Saving AE enhanced error model statistics');
    save(RunOptions.SaveFileNameAEErrorModel,'L_E','-v7.3')
end
if RunOptions.Load_AEErrorModel == 1
    printf('Loading AE enhanced error model statistics');
    load(RunOptions.SaveFileNameAEErrorModel,'L_E')
end
clear SampleMeanOpt SampleCovOpt
