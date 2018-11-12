function [hRecon,AcousticInverseItrtnInfo,Trmntn,PosteriorCovariance] = QPAT_EWE_DGM2D_AcousticInverse(RunOptions,DataVrblsWave,MeshD,MeshI,Prior,DGMMeshI,dt,PrecomputedIntrplteObjectsI,PLOT)

disp(' ')
disp('----------------------')
disp('QPAT Inverse Problem')
disp('----------------------')

%%== Shortening Labels ===%
NpI = DGMMeshI.Np;
KI = DGMMeshI.K;
xI = DGMMeshI.x;
yI = DGMMeshI.y;
NorderI = DGMMeshI.Norder;
pinfoI = DGMMeshI.pinfo;
rhoI = DGMMeshI.rho;
lambdaI = DGMMeshI.lambda;
muI = DGMMeshI.mu;
SensorsI = DataVrblsWave.SensorsI;
NumberofSensors = DataVrblsWave.NumberofSensors;
NumberofTimeSteps = RunOptions.NumberofTimeSteps;

%% =======================================================================%
%                          Noise Model and Prior
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Precompute Some Upwind Flux Terms for Non-Conforming Mesh %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RunOptions.UseMyRectangularMesh == 1
    MeshI.DomainIndices = 'Not constructed'; %Domain Indices
end
pinfoI = EWE_DGM2D_PrecomputeUpwindFluxPNonConf(RunOptions,pinfoI,NorderI,rhoI,lambdaI,muI,MeshI.DomainIndices,RunOptions.FluidDomainWithSolidLayerMeshI,RunOptions.SolidMeshI);   

%%%%%%%%%%%%%%%%%%%
%%% Noise Model %%%
%%%%%%%%%%%%%%%%%%%
%=== Adding Noise - Max minus Min of Data ===% Each of these have dimension 1 by NumberofSensors*NumberofTimeSteps 
if RunOptions.NoiseMinMax == 1;
    DataVrblsWave.ET11 = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.ET22 = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.ET12 = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.Evx = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.Evy = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.ET11 = RunOptions.Cov_ENoiseLevel*(max(DataVrblsWave.T11DataTimeSteps(:))-min(DataVrblsWave.T11DataTimeSteps(:)))*DataVrblsWave.ET11;
    DataVrblsWave.ET22 = RunOptions.Cov_ENoiseLevel*(max(DataVrblsWave.T22DataTimeSteps(:))-min(DataVrblsWave.T22DataTimeSteps(:)))*DataVrblsWave.ET22;
    DataVrblsWave.ET12 = RunOptions.Cov_ENoiseLevel*(max(DataVrblsWave.T12DataTimeSteps(:))-min(DataVrblsWave.T12DataTimeSteps(:)))*DataVrblsWave.ET12;
    DataVrblsWave.Evx = RunOptions.Cov_ENoiseLevel*(max(DataVrblsWave.vxDataTimeSteps(:))-min(DataVrblsWave.vxDataTimeSteps(:)))*DataVrblsWave.Evx;
    DataVrblsWave.Evy = RunOptions.Cov_ENoiseLevel*(max(DataVrblsWave.vyDataTimeSteps(:))-min(DataVrblsWave.vyDataTimeSteps(:)))*DataVrblsWave.Evy;
end

%=== Adding Noise - Max minus Min of Data at Sensor ===% 
if RunOptions.NoiseMinMaxS == 1;
    DataVrblsWave.MaxMinusMinT11 = zeros(1,DataVrblsWave.NumberofSensors);
    DataVrblsWave.MaxMinusMinT22 = zeros(1,DataVrblsWave.NumberofSensors);
    DataVrblsWave.MaxMinusMinT12 = zeros(1,DataVrblsWave.NumberofSensors);
    DataVrblsWave.MaxMinusMinvx = zeros(1,DataVrblsWave.NumberofSensors);
    DataVrblsWave.MaxMinusMinvy = zeros(1,DataVrblsWave.NumberofSensors);
    for s = 1:DataVrblsWave.NumberofSensors
        if RunOptions.FullqVectorData == 1;
            DataVrblsWave.MaxMinusMinST11(s) = max(DataVrblsWave.T11DataTimeSteps(s,:)) - min(DataVrblsWave.T11DataTimeSteps(s,:));
            DataVrblsWave.MaxMinusMinST22(s) = max(DataVrblsWave.T22DataTimeSteps(s,:)) - min(DataVrblsWave.T22DataTimeSteps(s,:));
            DataVrblsWave.MaxMinusMinST12(s) = max(DataVrblsWave.T12DataTimeSteps(s,:)) - min(DataVrblsWave.T12DataTimeSteps(s,:));
        end
        DataVrblsWave.MaxMinusMinSvx(s) = max(DataVrblsWave.vxDataTimeSteps(s,:)) - min(DataVrblsWave.vxDataTimeSteps(s,:));
        DataVrblsWave.MaxMinusMinSvy(s) = max(DataVrblsWave.vyDataTimeSteps(s,:)) - min(DataVrblsWave.vyDataTimeSteps(s,:));
    end
    DataVrblsWave.ET11 = repmat(DataVrblsWave.MaxMinusMinST11',NumberofTimeSteps,1);
    DataVrblsWave.ET22 = repmat(DataVrblsWave.MaxMinusMinST22',NumberofTimeSteps,1);
    DataVrblsWave.ET12 = repmat(DataVrblsWave.MaxMinusMinST12',NumberofTimeSteps,1);
    DataVrblsWave.Evx = repmat(DataVrblsWave.MaxMinusMinSvx',NumberofTimeSteps,1);
    DataVrblsWave.Evy = repmat(DataVrblsWave.MaxMinusMinSvy',NumberofTimeSteps,1);
    DataVrblsWave.ET11 = RunOptions.Cov_ENoiseLevel*DataVrblsWave.ET11;
    DataVrblsWave.ET22 = RunOptions.Cov_ENoiseLevel*DataVrblsWave.ET22;
    DataVrblsWave.ET12 = RunOptions.Cov_ENoiseLevel*DataVrblsWave.ET12;
    DataVrblsWave.Evx = RunOptions.Cov_ENoiseLevel*DataVrblsWave.Evx;
    DataVrblsWave.Evy = RunOptions.Cov_ENoiseLevel*DataVrblsWave.Evy;
end

%=== Adding Noise - Max of Data ===% Each of these have dimension 1 by NumberofSensors*NumberofTimeSteps 
if RunOptions.NoiseMax == 1;
    DataVrblsWave.ET11 = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.ET22 = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.ET12 = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.Evx = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.Evy = ones(NumberofSensors*NumberofTimeSteps,1);
    DataVrblsWave.ET11 = RunOptions.Cov_ENoiseLevel*max(abs(DataVrblsWave.T11DataTimeSteps(:)))*DataVrblsWave.ET11;
    DataVrblsWave.ET22 = RunOptions.Cov_ENoiseLevel*max(abs(DataVrblsWave.T22DataTimeSteps(:)))*DataVrblsWave.ET22;
    DataVrblsWave.ET12 = RunOptions.Cov_ENoiseLevel*max(abs(DataVrblsWave.T12DataTimeSteps(:)))*DataVrblsWave.ET12;
    DataVrblsWave.Evx = RunOptions.Cov_ENoiseLevel*max(abs(DataVrblsWave.vxDataTimeSteps(:)))*DataVrblsWave.Evx;
    DataVrblsWave.Evy = RunOptions.Cov_ENoiseLevel*max(abs(DataVrblsWave.vyDataTimeSteps(:)))*DataVrblsWave.Evy;
end
%=== Adding Noise - Max of Data At Each Sensor ===%
if RunOptions.NoiseMaxS == 1;
    DataVrblsWave.MaxST11 = zeros(1,DataVrblsWave.NumberofSensors);
    DataVrblsWave.MaxST22 = zeros(1,DataVrblsWave.NumberofSensors);
    DataVrblsWave.MaxST12 = zeros(1,DataVrblsWave.NumberofSensors);
    DataVrblsWave.MaxSvx = zeros(1,DataVrblsWave.NumberofSensors);
    DataVrblsWave.MaxSvy = zeros(1,DataVrblsWave.NumberofSensors);
    for s = 1:DataVrblsWave.NumberofSensors
        if RunOptions.FullqVectorData == 1;
            DataVrblsWave.MaxST11(s) = max(abs(DataVrblsWave.T11DataTimeSteps(s,:))) ;
            DataVrblsWave.MaxST22(s) = max(abs(DataVrblsWave.T22DataTimeSteps(s,:)));
            DataVrblsWave.MaxST12(s) = max(abs(DataVrblsWave.T12DataTimeSteps(s,:)));
        end
        DataVrblsWave.MaxSvx(s) = max(abs(DataVrblsWave.vxDataTimeSteps(s,:)));
        DataVrblsWave.MaxSvy(s) = max(abs(DataVrblsWave.vyDataTimeSteps(s,:)));
    end
    DataVrblsWave.ET11 = repmat(DataVrblsWave.MaxST11',NumberofTimeSteps,1);
    DataVrblsWave.ET22 = repmat(DataVrblsWave.MaxST22',NumberofTimeSteps,1);
    DataVrblsWave.ET12 = repmat(DataVrblsWave.MaxST12',NumberofTimeSteps,1);
    DataVrblsWave.Evx = repmat(DataVrblsWave.MaxSvx',NumberofTimeSteps,1);
    DataVrblsWave.Evy = repmat(DataVrblsWave.MaxSvy',NumberofTimeSteps,1);
    DataVrblsWave.ET11 = RunOptions.Cov_ENoiseLevel*DataVrblsWave.ET11;
    DataVrblsWave.ET22 = RunOptions.Cov_ENoiseLevel*DataVrblsWave.ET22;
    DataVrblsWave.ET12 = RunOptions.Cov_ENoiseLevel*DataVrblsWave.ET12;
    DataVrblsWave.Evx = RunOptions.Cov_ENoiseLevel*DataVrblsWave.Evx;
    DataVrblsWave.Evy = RunOptions.Cov_ENoiseLevel*DataVrblsWave.Evy;
end

%=== Error Mean ===%
T11ErrorMean = sparse(NumberofSensors,NumberofTimeSteps);
T22ErrorMean = sparse(NumberofSensors,NumberofTimeSteps);
T12ErrorMean = sparse(NumberofSensors,NumberofTimeSteps);
vxErrorMean = sparse(NumberofSensors,NumberofTimeSteps);
vyErrorMean = sparse(NumberofSensors,NumberofTimeSteps);

%=== Conventional Noise Model ===% %Still require L_ETij for the source term of the adjoint wave-field
L_ET11 = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,DataVrblsWave.ET11.^(-1));
L_ET22 = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,DataVrblsWave.ET22.^(-1));
L_ET12 = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,DataVrblsWave.ET12.^(-1));
L_Evx = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,DataVrblsWave.Evx.^(-1));
L_Evy = sparse(1:1:NumberofSensors*NumberofTimeSteps,1:1:NumberofSensors*NumberofTimeSteps,DataVrblsWave.Evy.^(-1));
if RunOptions.VelocitiesData == 1 %Since max(TijDataTimeSteps) = 0, then 1/0 = Inf. Instead, make them 0
    L_ET11 = 0*speye(NumberofSensors*NumberofTimeSteps);
    L_ET22 = 0*speye(NumberofSensors*NumberofTimeSteps);
    L_ET12 = 0*speye(NumberofSensors*NumberofTimeSteps);
end

%=== Approximation Error ===%
if RunOptions.AE_EWE == 1;
    EWE_DGM2D_ApproximationError_Driver
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Smoothness Prior %%%
%%%%%%%%%%%%%%%%%%%%%%%%
%=== \int ||A(x)\nabla u(x)|| dx ===%
if RunOptions.UseInformSmoothPrior == 1
    if RunOptions.EWE_LS_ExponentialTransform == 1;
        Prior.Exp_h = log(Prior.Exp_h);
        Prior.InformSmooth_Bounds_h(2) = log(220);
        Prior.InformSmooth_Bounds_h(1) = 1e-16;
    end
    disp('Computing informative smoothness prior')
    [Prior.Exp_h,Cov_pr,invCov_pr,Prior.traceCov_h,~]=SmoothnessPrior_Informative(MeshI.Nodes,MeshI.Elements,2*MeshI.Dimensns,Prior.InformSmooth_Bounds_h,Prior.InformSmooth_Corr_h,Prior.Exp_h,Prior.InformSmooth_Normalize,0,PLOT);
    Prior.L_pr = chol(invCov_pr); %Regularization operator
    clear invCov_pr
    if PLOT.PriorHist == 1;
        figure(PLOT.Figure_PriorHist)
        hist((chol(Prior.Gam))'*randn(MeshI.N_Nodes,1) + Prior.Exp_h);
        title(PLOT.Figure_PriorHist_Title,'FontWeight','bold')
    end
end
%=== Auto Correlation Prior ===%
if RunOptions.UseACPrior == 1
    disp('Computing autocorrelation smoothness prior')
    [Prior.L_pr,Prior.traceCov_h,Cov_pr,~] = SmoothnessPrior_AutoCorr(MeshI.Nodes,Prior.Exp_h,Prior.AC_Var_h,Prior.AC_Corr_h,0,PLOT);
end
%=== Standard Tikhonov Prior ===%
% Prior.L_pr = (1/30)*eye(MeshI.N_Nodes);
% Prior.traceCov_p0 = trace(Prior.L_pr);
L_Evx = 10000*L_Evx;
L_Evy = 10000*L_Evy;
%% =======================================================================%
%                         Line Search Iterations
%=========================================================================%
%=== Initial Estimate ===%
if RunOptions.EWE_LS_LogExpTransform == 1
    RunOptions.EWE_LS_hInitialGuess = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*RunOptions.EWE_LS_hInitialGuess)+1);
end
DataVrblsWave.hRecon = RunOptions.EWE_LS_hInitialGuess*ones(size(MeshI.Nodes,1),1);
% DataVrblsWave.hRecon = full(DataVrblsWave.h); %This is used to test line search function by considering the true h as the initial guess. Warning, this is on the DGMMeshD and so interpolation of this will be required when no inverse crime is committed
% Prior.Exp_h = DataVrblsWave.h; %To be used in adjoint-state method to check whether gradient is correct

%=== Reconstructing Initial Pressure ===%
if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensity == 1
[hRecon, Hyperprior_alpha, AcousticInverseItrtnInfo, Trmntn] = EWE_DGM2D_LineSearch_h(RunOptions,DataVrblsWave,@EWE_DGM2D_LSExpRK4,...
                                                                      L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                                      T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                                                      MeshD.N_Elm,MeshD.Nodes,MeshI.Elements,MeshI.N_Elm,MeshI.Nodes,MeshI.N_Nodes,Prior,...
                                                                      PrecomputedIntrplteObjectsI,pinfoI,xI*RunOptions.ScalingOptical,yI*RunOptions.ScalingOptical,...
                                                                      NpI,KI,rhoI,dt,PLOT);
end

%% =======================================================================%
%                         Uncertainty Estimates
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computing Uncertainty Estimates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('------------------------------------------------------')
disp('Computing Uncertainty Estimates')
disp('------------------------------------------------------')

SaveFileNameReconstructions = sprintf('Reconstructions-%s',RunOptions.SaveFileName);
load(SaveFileNameReconstructions,'AcousticInverseItrtnInfo','Trmntn')
hRecon = AcousticInverseItrtnInfo.hRecon(:,Trmntn);
RunOptions.EWE_DGM2D_UseGaussNewton = 1; RunOptions.EWE_DGM2D_UseNewton = 0; %Uncertainty estimate is constructed using Gauss-Newton direction, not Newton
if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 1;
    alphaRecon = AcousticInverseItrtnInfo.alphaRecon(Trmntn);
else
    alphaRecon = 1;
end

if RunOptions.ApproximatePosteriorCovariance == 1
    PosteriorCovariance = EWE_DGM2D_ApproxUncertaintyEstimates(RunOptions,@EWE_DGM2D_LSExpRK4,SensorsI,Cov_pr,hRecon,alphaRecon,...
                                                               L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                               Prior,...
                                                               MeshI.Elements,MeshI.N_Elm,MeshI.Nodes,MeshI.N_Nodes,...
                                                               PrecomputedIntrplteObjectsI,pinfoI,xI*RunOptions.ScalingOptical,yI*RunOptions.ScalingOptical,...
                                                               NpI,KI,rhoI,dt,PLOT);
else %Construct inverse of posterior explicitly covariance column by column and then invert
    T11AdjTime0 = 'Not Required';
    T22AdjTime0 = 'Not Required';
    if RunOptions.EWE_LS_SigmoidTransform == 1 || RunOptions.EWE_LS_ExponentialTransform == 1
        error('EWE_DGM2D_AdjStateHessianAction has not been set up for these positivity constraints')
    end
    invPosteriorCovariance = zeros(MeshI.N_Nodes,MeshI.N_Nodes);
    for ii=1:MeshI.N_Nodes
        BasisVector = zeros(MeshI.N_Nodes,1);
        BasisVector(ii) = 1;
        invPosteriorCovariance(:,ii) = EWE_DGM2D_AdjStateHessianAction(RunOptions,SensorsI,@EWE_DGM2D_LSExpRK4,BasisVector,hRecon,alphaRecon,...
                                                                       L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                                       T11AdjTime0,T22AdjTime0,...
                                                                       Prior,MeshI.Elements,MeshI.N_Elm,MeshI.Nodes,MeshI.N_Nodes,...
                                                                       pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,...
                                                                       ii,'Not Required',PLOT,'transp');  
    end
    PosteriorCovariance = inv(invPosteriorCovariance);
    RunOptions.SaveFileNamePostCov = sprintf('PostCov-%s',RunOptions.SaveFileName);
    save(RunOptions.SaveFileNamePostCov,'PosteriorCovariance','-v7.3')
    keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting Uncertainty Estimates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
hTrueIntrplte = IntrplteOver2DTriangulatedMesh(MeshD.N_Elm,MeshD.Nodes,DataVrblsWave.h,MeshI.Nodes(:,1),MeshI.Nodes(:,2),MeshI.N_Nodes,0);
%=== Loading Reconstructions and Posterior Covariance ===%
RunOptions.SaveFileNameReconstructions = sprintf('Reconstructions-%s',RunOptions.SaveFileName);
load(RunOptions.SaveFileNameReconstructions,'AcousticInverseItrtnInfo','Trmntn')
hRecon = full(AcousticInverseItrtnInfo.hRecon(:,Trmntn));
if RunOptions.ApproximatePosteriorCovariance == 0
    RunOptions.SaveFileNamePostCov = sprintf('PostCov-%s',RunOptions.SaveFileName);
    load(RunOptions.SaveFileNamePostCov,'PosteriorCovariance')
end
if RunOptions.ApproximatePosteriorCovariance == 1
    RunOptions.SaveFileNamePostCov = sprintf('ApproxPostCov-%s',RunOptions.SaveFileName);
    load(RunOptions.SaveFileNamePostCov,'ApproxPosteriorCovariance')
    PosteriorCovariance = ApproxPosteriorCovariance;
    clear ApproxPosteriorCovariance
end
UEDiag = diag(PosteriorCovariance);

%=== Forming Sample Posterior Covariance of h ===%
N_Samples = 5000;
[V,D] = eig(PosteriorCovariance);
L_PostCov = sqrt(D)*V';
z_s=bsxfun(@plus,hRecon,(L_PostCov.')*randn(length(MeshI.Nodes),N_Samples));
for ii=1:N_Samples
    h_s(:,ii) = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*z_s(:,ii))+1);
end
SampleMeanh = (1/N_Samples)*sum(h_s,2);
SamplePostCovh = 1/(N_Samples - 1)*(h_s*h_s') - (N_Samples/(N_Samples - 1))*(SampleMeanh*SampleMeanh');
UEDiag = diag(SamplePostCovh);

%=== Plotting Reconstructions ===%
PLOT.TRIFEM=delaunay(MeshI.Nodes(:,1),MeshI.Nodes(:,2));
figure
trisurf(PLOT.TRIFEM,MeshI.Nodes(:,1),MeshI.Nodes(:,2),hRecon);
if PLOT.DGMPlotBirdsEyeView == 1;
    view(2)
end
if PLOT.DGMPlotUsezLim == 1;
    zlim(PLOT.AbsorbedEnergyDensityzAxis);
end
caxis([0,200])
shading interp %thanks Ru!
colormap(jet(256))
hold on
%=== Plotting Cross-Section ===%
x_axis = -0.01:0.0005:0.01;
% y_axis = 0.5*x_axis + 0.005;
y_axis = 0.003*ones(size(x_axis));
plot3([x_axis(1) ; x_axis(end)],[y_axis(1) ; y_axis(end)],300*ones(2,1),'--k')
hold off

%=== Plotting Cross-Section With Confidence Intervals ===%
Coeff = zeros(3,1);
hTrueInterpolatedValues = zeros(size(x_axis,2),1);
hReconInterpolatedValues = zeros(size(x_axis,2),1);
PostUEInterpolatedValues = zeros(size(x_axis,2),1);
for ii = 1:size(x_axis,2);
    K = pointLocation(PrecomputedIntrplteObjectsI.DT,[x_axis(ii),y_axis(ii)]); %returns the element that the gridpoint lies in.
    Node_1 = PrecomputedIntrplteObjectsI.DT(K,1);
    Node_2 = PrecomputedIntrplteObjectsI.DT(K,2);
    Node_3 = PrecomputedIntrplteObjectsI.DT(K,3);
    
    %Interpolated Values for hTrue
    z = [hTrueIntrplte(Node_1); hTrueIntrplte(Node_2); hTrueIntrplte(Node_3)];
    B = cell2mat(PrecomputedIntrplteObjectsI.PrecompCoeffMatricesInv(K));
    Coeff = B*z;
    hTrueInterpolatedValues(ii) = Coeff(1) + Coeff(2)*x_axis(ii) + Coeff(3)*y_axis(ii);
    
    %Interpolated Values for hRecon
    z = [hRecon(Node_1); hRecon(Node_2); hRecon(Node_3)];
    B = cell2mat(PrecomputedIntrplteObjectsI.PrecompCoeffMatricesInv(K));
    Coeff = B*z;
    hReconInterpolatedValues(ii) = Coeff(1) + Coeff(2)*x_axis(ii) + Coeff(3)*y_axis(ii);
    
    %Interpolated Values for Uncertainty Estimates
    z = [UEDiag(Node_1); UEDiag(Node_2); UEDiag(Node_3)];
    B = cell2mat(PrecomputedIntrplteObjectsI.PrecompCoeffMatricesInv(K));
    Coeff = B*z;
    PostUEInterpolatedValues(ii) = Coeff(1) + Coeff(2)*x_axis(ii) + Coeff(3)*y_axis(ii);
end

figure
ciplot(hReconInterpolatedValues-sqrt(PostUEInterpolatedValues),hReconInterpolatedValues+sqrt(PostUEInterpolatedValues),x_axis,'k');
hold on
ciplot(hReconInterpolatedValues-2*sqrt(PostUEInterpolatedValues),hReconInterpolatedValues+2*sqrt(PostUEInterpolatedValues),x_axis,'k');
ciplot(hReconInterpolatedValues-3*sqrt(PostUEInterpolatedValues),hReconInterpolatedValues+3*sqrt(PostUEInterpolatedValues),x_axis,'k');
plot(x_axis,hTrueInterpolatedValues,'r');
plot(x_axis,hReconInterpolatedValues,'b');
ylim([-500,500]);
keyboard




