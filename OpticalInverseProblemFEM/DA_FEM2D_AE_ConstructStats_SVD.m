function [SampleMeanOpt,SampleCovOpt]=DA_FEM2D_AE_ConstructStats_SVD(RunOptions,DataVrblsOptical,MeshI,PrmtrsI,PrmtrsPrp,Prior,DGMMeshI,PrecomputedIntrplteObjectsI,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,SensorsI,dt,PLOT)

% AEScattering calculates the modelling error mean and modelling error 
% covariance for the enhanced error model. Specifically, it is to account
% for errors arising from using an estimate for the scattering coefficient.
%
% Inputs:
%   MeshA - Accurate optical mesh properties
%   MeshI - Inverse optical mesh properties
%   Data:
% %      Gruneisen - Grunaisen parameter
%   PrmtrsI:
%      mu_s - Estimated scattering coefficient
% %   PrmtrsPrp:
%      g - mean of the cosine of the scattering angle
%   Prior:
%      InformSmooth_Corr_mu_a - 2*1 vector = [correlation_x, correlation_y]. The larger the number,
%          the more points near to a marginalisation point are correlated
%          in the x and y direction.
%      InformSmooth_Bounds_mu_a - approximately [mu_a_min, mu_a_max]
%      Exp_mu_a - initial guess of mu_a 
%      InformSmooth_Corr_mu_w - 2*1 vector = [correlation_x, correlation_y]. The larger the number,
%          the more points near to a marginalisation point are correlated
%          in the x and y direction.
%      InformSmooth_Bounds_mu_w - approximately [mu_w_min, mu_w_max]
%      Exp_mu_w - initial guess of mu_a 
%
% Outputs:
%      SampleMeanOpt - modelling error mean
%      SampleCovOpt - modelling error covariance
%
% Hwan Goh 30/03/2018, University of Auckland, New Zealand

PLOT.DGMForward = 0; %Suppress plotting of wave propagation
PLOT.PriorMargPoints = 0; %Suppress plotting of marginalisation points
PLOT.PriorSamples = 0; %Suppress plotting of prior samples
PLOT.DGMForwardSensorData = 0; %Suppress plotting of sensory data
N_Samples = RunOptions.AE_DA_N_Samples;
printf([num2str(N_Samples) ' samples to be computed for approximation error (SVD) ']);
SaveFileNameAESamples = sprintf('AESamplesOpt-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%sCorr-%dAESamples-%sAE',RunOptions.SaveFileNameDomain,RunOptions.SaveFileNameStateVector,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.SaveFileNamePriorCorr,RunOptions.AE_DA_N_Samples,RunOptions.TrelisMeshAElementSize);
pinfoI = DGMMeshI.pinfo;
NorderI = DGMMeshI.Norder;
rhoI = DGMMeshI.rho;
lambdaI = DGMMeshI.lambda;
muI = DGMMeshI.mu;
VertexNodesGlobalIndicesFEMI = DGMMeshI.VertexNodesGlobalIndicesFEM;

%=========================================================================%
%                         Generating Samples
%=========================================================================%
%% %%%%%%%%%%%%%%%%%%%
%%% Accurate Model %%%
%%%%%%%%%%%%%%%%%%%%%%
printf('Computing samples of accurate model (SVD)');

%=== Initial Pressure Samples ===%
[~,~,~,~,mu_aAS]=SmoothnessPrior_Informative(MeshI.Nodes,MeshI.Elements,2*MeshI.Dimensns,Prior.InformSmooth_Bounds_mu_a,Prior.InformSmooth_Corr_mu_a,Prior.Exp_mu_a,Prior.InformSmooth_Normalize,N_Samples,PLOT);

%=== Forward Data Samples ===%
for n=1:N_Samples;
    printf(['\nComputing Accurate Model Sample ' num2str(n) ' of ' num2str(N_Samples)]);
    printf(['For the Case ' RunOptions.SaveFileName]);
    %absorption coefficient
    PrmtrsAS.mu_a = mu_aAS(:,n);
    PrmtrsAS.mu_a_elmts = reshapeNodes(MeshI,PrmtrsAS.mu_a);
    %diffusion coefficient
    PrmtrsAS.kappa = zeros(MeshI.N_Nodes,1);
    for j=1:MeshI.N_Nodes;
        PrmtrsAS.kappa(j) = 1/(2*(PrmtrsAS.mu_a(j) + ((1-PrmtrsPrp.g)*PrmtrsI.mu_sEst)));
    end
    PrmtrsAS.kappa_elmts = reshapeNodes(MeshI,PrmtrsAS.kappa);
    %Light fluence
    [phiAS,~,~,~,~,~] = DA_FEM2D_OpticalForward(RunOptions,MeshI,PrmtrsAS);
    %Initial pressure distribution
    p0AStemp = DataVrblsOptical.Gruneisen*((PrmtrsAS.mu_a).*(phiAS)/(RunOptions.SpecificHeatCoeff));
    p0AS(:,n) = p0AStemp(:);
end
save(SaveFileNameAESamples,'mu_aAS','p0AS','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%
%%% Less Accurate Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf('Computing samples of less accurate model (SVD)');
%=== Precomputable Stuff for Acoustic Forward Problem ===%
pinfoI = EWE_DGM2D_PrecomputeUpwindFluxPNonConf(RunOptions,pinfoI,NorderI,rhoI,lambdaI,muI,MeshI.DomainIndices,RunOptions.FluidDomainWithSolidLayerMeshI,RunOptions.SolidMeshI);   
MassMatrixFEM = ConstructMassMatrixFEM2D(MeshI.Elements,MeshI.Nodes);
PLOT.TRI_DGMMeshD = delaunay(DGMMeshI.x,DGMMeshI.y); 
PLOT.TRI_DGMMeshI = PLOT.TRI_DGMMeshD;
DataVrblsWave.NumberofSensors = size(SensorsI,2);
DataVrblsWave.SensorsI = SensorsI;
%=== \int ||A(x)\nabla u(x)|| dx ===%
[Prior.Exp_p0,~,invCov,Prior.traceCov_p0,~]=SmoothnessPrior_Informative(MeshI.Nodes,MeshI.Elements,2*MeshI.Dimensns,Prior.InformSmooth_Bounds_p0,Prior.InformSmooth_Corr_p0,Prior.Exp_p0,0,PLOT);
Prior.L_pr = chol(invCov); %Regularization operator
clear invCov
if PLOT.PriorHist == 1;
    figure(PLOT.Figure_PriorHist)
    hist((chol(Prior.Gam))'*randn(MeshI.N_Nodes,1) + Prior.Exp_p0);
    title(PLOT.Figure_PriorHist_Title,'FontWeight','bold')
end
DataVrblsWave.p0Recon = RunOptions.EWE_SD_p0InitialGuess*ones(size(MeshI.Nodes,1),1);

%=== Forward Data Samples ===%
for n=1:N_Samples;
    [p0DGMIS,~] = IntrplteOver2DTriangulatedMesh(MeshI.N_Elm,MeshI.Nodes*(1/RunOptions.ScalingOptical),p0AS(:,n),DGMMeshI.x,DGMMeshI.y,DGMMeshI.Np*DGMMeshI.K,PrecomputedIntrplteObjectsI);
    p0DGMIS = reshape(p0DGMIS,DGMMeshI.Np,DGMMeshI.K)*RunOptions.ScalingOptical;
    IniCond = p0DGMIS./(rhoI*RunOptions.SpecificHeatCoeff);
    %=== Acoustic Forward Problem ===%
    [T11TimeSteps,T22TimeSteps,T12TimeSteps,vxTimeSteps,vyTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,DGMMeshI.x,DGMMeshI.y,DGMMeshI.Np,DGMMeshI.K,pinfoI,dt,PLOT);
    DataVrblsWave.T11DataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    DataVrblsWave.T22DataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    DataVrblsWave.T12DataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    DataVrblsWave.vxDataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    DataVrblsWave.vyDataTimeSteps = sparse(DataVrblsWave.NumberofSensors,RunOptions.NumberofTimeSteps);
    %Interpolation to Form Sensory Data
    for t=1:RunOptions.NumberofTimeSteps
        for s=1:DataVrblsWave.NumberofSensors
            if RunOptions.FullqVectorData == 1;
                DataVrblsWave.T11DataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T11TimeSteps(SensorsI{s}.id,t);
                DataVrblsWave.T22DataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T22TimeSteps(SensorsI{s}.id,t);
                DataVrblsWave.T12DataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*T12TimeSteps(SensorsI{s}.id,t);
            end
            DataVrblsWave.vxDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vxTimeSteps(SensorsI{s}.id,t);
            DataVrblsWave.vyDataTimeSteps(s,t) = SensorsI{s}.l_iatsensor*vyTimeSteps(SensorsI{s}.id,t);
        end
    end
    %=== Acoustic Inverse Problem ===%
    %Initial Estimate
    DataVrblsWave.p0 = p0AS(:,n);
    %Computing the Iterations
    [p0IS(:,n),~,~] = EWE_DGM2D_LineSearch_SteepestDescent(RunOptions,DataVrblsWave,@EWE_DGM2D_LSExpRK4,...
                                                           L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                           T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                                           MeshI.N_Elm,MeshI.Nodes,MeshI.Elements,MeshI.N_Elm,MeshI.Nodes,MeshI.N_Nodes,Prior,MassMatrixFEM,...
                                                           PrecomputedIntrplteObjectsI,pinfoI,DGMMeshI.x*RunOptions.ScalingOptical,DGMMeshI.y*RunOptions.ScalingOptical,...
                                                           DGMMeshI.Np,DGMMeshI.K,rhoI,VertexNodesGlobalIndicesFEMI,dt,PLOT);
    clear T11TimeSteps T22TimeSteps T12TimeSteps vxTimeSteps vyTimeSteps 
end
clear DataVrblsWave L_ET11 L_ET22 L_ET12 L_Evx L_Evy T11ErrorMean T22ErrorMean T12ErrorMean vxErrorMean vyErrorMean DGMMeshI PrecomputedIntrplteObjectsI MassMatrixFEM
save(SaveFileNameAESamples,'mu_aAS','p0AS','mu_aIS','p0IS','-v7.3')

%% =======================================================================%
%                         Approximation Error
%=========================================================================%
%=== Computing Error Vectors ===%
error = p0AS - p0IS;

%=== Computing Sample Means ===%
printf('Computing sample means');
SampleMeanOpt = (1/N_Samples)*sum(error,2);

%=== Computing Sample Covariances ===%
printf('Computing sample covariances (SVD)');
W = 1/(sqrt(N_Samples - 1))*bsxfun(@minus,error,SampleMeanOpt);
[U,S,~] = svd(W);
US = U*S;

SampleCovOpt = US*US';

save(SaveFileNameAESamples,'mu_aAS','p0AS','mu_aIS','p0IS','W','U','S','-v7.3')
                          
printf('Computation of sample means and covariances with SVD complete');
 