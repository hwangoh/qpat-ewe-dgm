function [T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy] = EWE_DGM2D_AE_ConstructStats_QR(RunOptions,NumberofSensors,hIS_FEM,ET11,ET22,ET12,Evx,Evy,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,errorT11,errorT22,errorT12,errorvx,errorvy)

% EWE_DGM2D_AE_ConstructStats_QR calculates the approximation error mean and approximation error 
% covariance for the enhanced error model. This is done using QR decomposition.
%
% Inputs:
%   RunOptions
%   NumberofSensors
%   hIS_FEM - Prior model sample draws, mainly used for QR construction of
%              approximation error statistics for full AE error model (not enhanced)
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
EffectiveRank = 300;

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

SampleMeanT11 = reshape(SampleMeanT11,NumberofSensors,NumberofTimeSteps);
SampleMeanT22 = reshape(SampleMeanT22,NumberofSensors,NumberofTimeSteps);
SampleMeanT12 = reshape(SampleMeanT12,NumberofSensors,NumberofTimeSteps);
SampleMeanvx = reshape(SampleMeanvx,NumberofSensors,NumberofTimeSteps);
SampleMeanvy = reshape(SampleMeanvy,NumberofSensors,NumberofTimeSteps);

SampleMeanh = (1/N_Samples)*sum(hIS_FEM,2);

%=== Computing QR Decompositions ===%
printf('Computing QR decompositions of sample errors');

errorT11 = bsxfun(@minus,errorT11,SampleMeanT11(:));
errorT22 = bsxfun(@minus,errorT22,SampleMeanT22(:));
errorT12 = bsxfun(@minus,errorT12,SampleMeanT12(:));
errorvx = bsxfun(@minus,errorvx,SampleMeanvx(:));
errorvy = bsxfun(@minus,errorvy,SampleMeanvy(:));

hIS_FEM = bsxfun(@minus,hIS_FEM,SampleMeanh(:));

tic
[Qh Rh] = qr(hIS_FEM);
[QT11 RT11] = qr(errorT11);
[QT22 RT22] = qr(errorT22);
[QT12 RT12] = qr(errorT12);
[Qvx Rvx] = qr(errorvx);
[Qvy Rvy] = qr(errorvy);
toc

Qh = Qh(:,1:EffectiveRank);
QT11 = QT11(:,1:EffectiveRank);
QT22 = QT22(:,1:EffectiveRank);
QT12 = QT12(:,1:EffectiveRank);
Qvx = Qvx(:,1:EffectiveRank);
Qvy = Qvy(:,1:EffectiveRank);

Rh = Rh(1:EffectiveRank,:);
RT11 = RT11(1:EffectiveRank,:);
RT22 = RT22(1:EffectiveRank,:);
RT12 = RT12(1:EffectiveRank,:);
Rvx = Rvx(1:EffectiveRank,:);
Rvy = Rvy(1:EffectiveRank,:);

%=== Saving and Loading AE Sample Stats ===%
printf('Saving AE sample statistics QR:');
save(RunOptions.SaveFileNameAEStats,'SampleMeanT11','SampleMeanT22','SampleMeanT12','SampleMeanvx','SampleMeanvy',...
                                          'Qh','Rh','QT11','RT11','QT22','RT22','QT12','RT12','Qvx','Rvx','Qvy','Rvy',...
                                          '-v7.3')
end
if RunOptions.Load_AESampleMeans == 1
    printf('Loading AE sample means:');
    load(RunOptions.SaveFileNameAEStats,'SampleMeanT11','SampleMeanT22','SampleMeanT12','SampleMeanvx','SampleMeanvy')
end
if RunOptions.Load_AESampleCovs == 1
    printf('Loading AE sample covariances QR:');
    load(RunOptions.SaveFileNameAEStats,'Qh','Rh','QT11','RT11','QT22','RT22','QT12','RT12','Qvx','Rvx','Qvy','Rvy')
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
L_NT11 = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,ET11.^(-1)); %This is the noise model, usually denoted as L_E but we change it here to avoid cholupdate updating the wrong thing
L_NT22 = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,ET22.^(-1));
L_NT12 = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,ET12.^(-1));
L_Nvx = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,Evx.^(-1));
L_Nvy = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,Evy.^(-1));

if RunOptions.Save_AEErrorModel == 1
    printf('Computing matrix square roots of sample covariances from QR decompositions');
    tic
    if RunOptions.FullqVectorData == 1;
        MT11 = (1/(N_Samples - 1))*(RT11*RT11');
        MT22 = (1/(N_Samples - 1))*(RT22*RT22');
        MT12 = (1/(N_Samples - 1))*(RT12*RT12');        
        L_MT11 = inv(chol(inv(MT11) + QT11'*(L_NT11'*L_NT11)*QT11));
        L_MT22 = inv(chol(inv(MT22) + QT22'*(L_NT22'*L_NT22)*QT22));
        L_MT12 = inv(chol(inv(MT12) + QT11'*(L_NT12'*L_NT12)*QT12));        
        CT11 = ((L_NT11'*L_NT11)*QT11*L_MT11);
        CT22 = ((L_NT22'*L_NT22)*QT22*L_MT22);
        CT12 = ((L_NT12'*L_NT12)*QT12*L_MT12);        
        for ii=1:size(CT11,2);
            L_ET11 = cholupdate(full(L_NT11),CT11(:,ii),'-');
            L_ET22 = cholupdate(full(L_NT22),CT22(:,ii),'-');
            L_ET12 = cholupdate(full(L_NT12),CT12(:,ii),'-');
        end
    else
        L_ET11 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
        L_ET22 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
        L_ET12 = sparse(NumberofSensors*NumberofTimeSteps,NumberofSensors*NumberofTimeSteps);
    end
    Mvx = (1/(N_Samples - 1))*(Rvx*Rvx');
    Mvy = (1/(N_Samples - 1))*(Rvy*Rvy');    
    L_Mvx = inv(chol(inv(Mvx) + Qvx'*(L_Nvx'*L_Nvx)*Qvx));
    L_Mvy = inv(chol(inv(Mvy) + Qvy'*(L_Nvy'*L_Nvy)*Qvy));    
    Cvx = ((L_Nvx'*L_Nvx)*Qvx*L_Mvx);
    Cvy = ((L_Nvy'*L_Nvy)*Qvy*L_Mvy);   
    for ii=1:size(Cvx,2);
        disp(' ')
        disp('------------------------------------------------------')
        disp('Forming Matrix Square Roots Using Cholupdate')
        disp('------------------------------------------------------')
        printf(['For the Case ' RunOptions.SaveFileName]);
        printf(['Column Number: ' num2str(ii) ' of ' num2str(size(Cvx,2))]);
        tic
        L_Evx = cholupdate(full(L_Nvx),Cvx(:,ii),'-');
        L_Evy = cholupdate(full(L_Nvy),Cvy(:,ii),'-');
        toc
    end
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
