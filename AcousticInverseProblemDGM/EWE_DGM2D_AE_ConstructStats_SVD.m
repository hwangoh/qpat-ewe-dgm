function [T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy] = EWE_DGM2D_AE_ConstructStats_SVD(RunOptions,NumberofSensors,ET11,ET22,ET12,Evx,Evy,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,errorT11,errorT22,errorT12,errorvx,errorvy)

% EWE_DGM2D_AE_ConstructStats_SVD calculates the approximation error mean and approximation error 
% covariance for the enhanced error model. It is also used to compute the noise model
%
% Inputs:
%   RunOptions
%   NumberofSensors
%   ETij, Evi - Standard deviation of noise
%   TijErrorMean, viErrorMean - Noise model means
%   errorTij, errorvi - approximation error samples
%
% Outputs:
%      TijErrorMean, viErrorMean - approximation error means
%      L_ETij, L_Evi - Approximation error covariances
%
% Hwan Goh 27/07/2018, University of Auckland, New Zealand

N_Samples = RunOptions.AEN_Samples;
NumberofTimeSteps = RunOptions.NumberofTimeSteps;

%% =======================================================================%
%                     Approximation Error Statistics
%=========================================================================%
if RunOptions.Save_AEStats == 1
%=== Computing Sample Means ===%
printf('Computing sample means');
SampleMeanT11 = (1/N_Samples)*sum(errorT11,2);
SampleMeanT22 = (1/N_Samples)*sum(errorT22,2);
SampleMeanT12 = (1/N_Samples)*sum(errorT12,2);
SampleMeanvx = (1/N_Samples)*sum(errorvx,2);
SampleMeanvy = (1/N_Samples)*sum(errorvy,2);

%=== Computing Sample Covariances ===%
printf('Computing sample covariances (SVD)');
if RunOptions.FullqVectorData == 1;
    WT11 = 1/(sqrt(N_Samples - 1))*bsxfun(@minus,errorT11,SampleMeanT11);
    WT22 = 1/(sqrt(N_Samples - 1))*bsxfun(@minus,errorT22,SampleMeanT22);
    WT12 = 1/(sqrt(N_Samples - 1))*bsxfun(@minus,errorT12,SampleMeanT12);
else
    WT11 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
    WT22 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
    WT12 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
end
Wvx = 1/(sqrt(N_Samples - 1))*bsxfun(@minus,errorvx,SampleMeanvx);
Wvy = 1/(sqrt(N_Samples - 1))*bsxfun(@minus,errorvy,SampleMeanvy);

if RunOptions.FullqVectorData == 1;
    [UT11,ST11,~] = svd(WT11);
    [UT22,ST22,~] = svd(WT22);
    [UT12,ST12,~] = svd(WT12);
else
    [UT11,ST11,~] = svds(WT11);
    [UT22,ST22,~] = svds(WT22);
    [UT12,ST12,~] = svds(WT12);
end
[Uvx,Svx,~] = svd(Wvx);
[Uvy,Svy,~] = svd(Wvy);

UST11 = UT11*ST11;
UST22 = UT22*ST22;
UST12 = UT12*ST12;
USvx = Uvx*Svx;
USvy = Uvy*Svy;

SampleCovT11 = UST11*UST11';
SampleCovT22 = UST22*UST22';
SampleCovT12 = UST12*UST12';
SampleCovvx = USvx*USvx';
SampleCovvy = USvy*USvy';

%=== Reshaping sample means to be used later ===%
SampleMeanT11 = reshape(SampleMeanT11,NumberofSensors,NumberofTimeSteps);
SampleMeanT22 = reshape(SampleMeanT22,NumberofSensors,NumberofTimeSteps);
SampleMeanT12 = reshape(SampleMeanT12,NumberofSensors,NumberofTimeSteps);
SampleMeanvx = reshape(SampleMeanvx,NumberofSensors,NumberofTimeSteps);
SampleMeanvy = reshape(SampleMeanvy,NumberofSensors,NumberofTimeSteps);

%=== Saving SVD Objects ===%
if RunOptions.FullqVectorData == 1; %since ST_ij and Sv_i are sparse matrices, might as well save them sparse
    ST11 = sparse(ST11);
    ST22 = sparse(ST22);
    ST12 = sparse(ST12);
end
Svx = sparse(Svx);
Svy = sparse(Svy);

%=== Saving and Loading AE Sample Stats ===%
printf('Saving AE sample statistics:'); %NOTE old version saves SVD components in 'AESamples' MAT file!
save(RunOptions.SaveFileNameAEStats,'SampleMeanT11','SampleMeanT22','SampleMeanT12','SampleMeanvx','SampleMeanvy',...
                                          'SampleCovT11','SampleCovT22','SampleCovT12','SampleCovvx','SampleCovvy',...
                                          'WT11','WT22','WT12','Wvx','Wvy',...
                                          'UT11','UT22','UT12','Uvx','Uvy',...
                                          'ST11','ST22','ST12','Svx','Svy',...
                                          '-v7.3')                                            
end
if RunOptions.Load_AESampleMeans == 1
    printf('Loading AE sample means:');
    load(RunOptions.SaveFileNameAEStats,'SampleMeanT11','SampleMeanT22','SampleMeanT12','SampleMeanvx','SampleMeanvy')
end
if RunOptions.Load_AESampleCovs == 1
    printf('Loading AE sample covariances:');
    load(RunOptions.SaveFileNameAEStats,'SampleCovT11','SampleCovT22','SampleCovT12','SampleCovvx','SampleCovvy') 
end
                          
printf('Computation of sample means and covariances with SVD complete');

%% =======================================================================%
%                           Computing Noise Model
%=========================================================================%
%=== Enhanced Error Model Mean ===%
T11ErrorMean = T11ErrorMean + SampleMeanT11;
T22ErrorMean = T22ErrorMean + SampleMeanT22;
T12ErrorMean = T12ErrorMean + SampleMeanT12;
vxErrorMean = vxErrorMean + SampleMeanvx;
vyErrorMean = vyErrorMean + SampleMeanvy;

%=== Enhanced Error Model Covariance Matrix Square Root ===%
if RunOptions.Save_AEErrorModel == 1
    printf('Computing Cholesky decomposition of sample covariances');
    tic
    if RunOptions.FullqVectorData == 1;
        L_ET11 = inv(chol(diag(ET11.^2) + SampleCovT11));
        L_ET22 = inv(chol(diag(ET22.^2) + SampleCovT22));
        L_ET12 = inv(chol(diag(ET12.^2) + SampleCovT12));
    else
        L_ET11 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
        L_ET22 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
        L_ET12 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
    end
    L_Evx = inv(chol(diag(Evx.^2) + SampleCovvx));
    L_Evy = inv(chol(diag(Evy.^2) + SampleCovvy));
    toc
    
    %=== Saving and Loading AE Enhanced Error Model Statistics ===%
    L_Evx = sparse(L_Evx);
    L_Evy = sparse(L_Evy);
    printf('Saving AE enhanced error model covariances');
    save(RunOptions.SaveFileNameAEErrorModel,'L_ET11','L_ET22','L_ET12','L_Evx','L_Evy','-v7.3')
end
if RunOptions.Load_AEErrorModel == 1
    printf('Loading AE enhanced error model covariances');
    load(RunOptions.SaveFileNameAEErrorModel,'L_ET11','L_ET22','L_ET12','L_Evx','L_Evy')
end
