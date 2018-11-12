function [mu_aRecon,OpticalInverseItrtnInfo,Trmntn] = QPAT_DA_FEM2D_OpticalInverse(RunOptions,DataVrblsOptical,MeshD,MeshI,PrmtrsD,PrmtrsI,PrmtrsPrp,Prior,DGMMeshI,PrecomputedIntrplteObjectsI,SensorsI,dt,PLOT)

% QPAT_DA_FEM2D_OpticalInverse is the main function file for the optical inverse problem

%% =======================================================================%
%                          Noise Model and Prior
%=========================================================================%                 
%%%%%%%%%%%%%%%%%%%
%%% Noise Model %%%
%%%%%%%%%%%%%%%%%%%
%=== Conventional Noise Model ===%
DataVrblsOptical.E = RunOptions.Cov_ENoiseLevel*max(DataVrblsOptical.p0Recon(:))*ones(MeshI.N_Nodes,1);
DataVrblsOptical.invE = DataVrblsOptical.E.^(-1);
DataVrblsOptical.L_E = sparse(1:1:MeshI.N_Nodes,1:1:MeshI.N_Nodes,DataVrblsOptical.E.^(-1));

if RunOptions.DA_GN_CoupledErrorModel == 1
    RunOptions.SaveFileNameReconstructions = sprintf('Reconstructions-%s',RunOptions.SaveFileName);
    load(RunOptions.SaveFileNameReconstructions,'AcousticInverseItrtnInfo','Trmntn')
    p0Recon = full(AcousticInverseItrtnInfo.p0Recon(:,Trmntn));
    RunOptions.SaveFileNamePostCov = sprintf('PostCov-%s',RunOptions.SaveFileName);
    load(RunOptions.SaveFileNamePostCov,'PosteriorCovariance')
    AcoustPostCov = PosteriorCovariance;
    clear PosteriorCovariance
    %=== Forming Sample Posterior Covariance of p0 ===%
    N_Samples = 20000;
    [V,D] = eig(AcoustPostCov); 
    L_PostCov = sqrt(D)*V';
    z_s=bsxfun(@plus,p0Recon,(L_PostCov.')*randn(length(MeshI.Nodes),N_Samples));
    for ii=1:N_Samples
        p0_s(:,ii) = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*z_s(:,ii))+1);
    end
    SampleMeanp0 = (1/N_Samples)*sum(p0_s,2);
    AcoustPostCov = 1/(N_Samples - 1)*(p0_s*p0_s') - (N_Samples/(N_Samples - 1))*(SampleMeanp0*SampleMeanp0');
    [V,D] = eig(inv(AcoustPostCov));
    DataVrblsOptical.L_E = sqrt(D)*V';
end

%=== Approximation Error ===%
printf('Constructing noise regularization matrix with approximation error');
if RunOptions.AE_DA == 1;
    DA_FEM2D_ApproximationError_Driver
end
                                    
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Smoothness Prior %%%
%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing smoothness prior')
[Prior.Exp_mu_a,~,invCov,Prior.traceCov_mu_a,~]=SmoothnessPrior_Informative(MeshI.Nodes,MeshI.Elements,2*MeshI.Dimensns,Prior.InformSmooth_Bounds_mu_a,Prior.InformSmooth_Corr_mu_a,Prior.Exp_mu_a,Prior.InformSmooth_Normalize,0,PLOT);
Prior.L_pr = chol(invCov); %Regularization operator

%% =======================================================================%
%                          Gauss-Newton Iterations
%=========================================================================%
%=== Initial Estimate ===%
PrmtrsI.mu_a = Prior.Exp_mu_a;

%=== Computing the Iterations ===%
[mu_aRecon,OpticalInverseItrtnInfo,Trmntn] = DA_FEM2D_LineSearch_GaussNewton(RunOptions,DataVrblsOptical,@DA_FEM2D_GNFwrd,@DA_FEM2D_GNJacobian,...
                                                                            [DataVrblsOptical.L_E*DataVrblsOptical.p0Recon(:); Prior.L_pr*Prior.Exp_mu_a],DataVrblsOptical.L_E,Prior,...
                                                                             PrmtrsD,PrmtrsI,PrmtrsPrp,MeshD,MeshI,PLOT);

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

SaveFileNameReconstructions = sprintf('Reconstructions-Optical-%s',RunOptions.SaveFileName);
load(SaveFileNameReconstructions,'OpticalInverseItrtnInfo','Trmntn')
mu_aRecon = OpticalInverseItrtnInfo.mu_aRecon(:,Trmntn);

%=== Jacobian at mu_aRecon ===%
[~,PrmtrsRecon,phiRecon,ARecon,~,~,~] = DA_FEM2D_GNFwrd(RunOptions,DataVrblsOptical,MeshI,mu_aRecon,PrmtrsI.mu_s,PrmtrsPrp,'To be computed','To be computed');
Jmu_aRecon = DA_FEM2D_GNJacobian(MeshI,PrmtrsRecon,phiRecon,ARecon);

%=== Posterior Covariance ===%
invPosteriorCovariance = Jmu_aRecon'*DataVrblsOptical.L_E'*DataVrblsOptical.L_E*Jmu_aRecon + Prior.L_pr'*Prior.L_pr;
if RunOptions.DA_GN_CoupledErrorModel == 1
    Cov_E = AcoustPostCov;
    invPosteriorCovariance = Jmu_aRecon'*inv(Cov_E)*Jmu_aRecon + Prior.L_pr'*Prior.L_pr;  
end
PosteriorCovariance = inv(invPosteriorCovariance);
clear invPosteriorCovariance
RunOptions.SaveFileNamePostCov = sprintf('PostCov-Optical-%s',RunOptions.SaveFileName);
save(RunOptions.SaveFileNamePostCov,'PosteriorCovariance','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting Uncertainty Estimates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
[mu_aTrueIntrplte] = IntrplteOver2DTriangulatedMesh(MeshD.N_Elm,MeshD.Nodes,PrmtrsD.mu_a,MeshI.Nodes(:,1),MeshI.Nodes(:,2),MeshI.N_Nodes,0);
%=== Loading Reconstructions and Posterior Covariance ===%
RunOptions.SaveFileNameReconstructions = sprintf('Reconstructions-Optical-%s',RunOptions.SaveFileName);
load(RunOptions.SaveFileNameReconstructions,'OpticalInverseItrtnInfo','Trmntn')
mu_aRecon = full(OpticalInverseItrtnInfo.mu_aRecon(:,Trmntn));
RunOptions.SaveFileNamePostCov = sprintf('PostCov-Optical-%s',RunOptions.SaveFileName);
load(RunOptions.SaveFileNamePostCov,'PosteriorCovariance')
UEDiag = diag(PosteriorCovariance);

%=== Plotting Reconstructions ===%
PLOT.TRIFEM=delaunay(MeshI.Nodes(:,1),MeshI.Nodes(:,2));
figure
trisurf(PLOT.TRIFEM,MeshI.Nodes(:,1),MeshI.Nodes(:,2),mu_aRecon*10^-3);
if PLOT.DGMPlotBirdsEyeView == 1;
    view(2)
end
zlim([-0.1,0.5])
shading interp %thanks Ru!
colormap(jet(256))
caxis([0,0.25])
% title(PLOT.Figure_Currentmu_aRecon_Title,'FontWeight','bold')
hold on
keyboard
%=== Plotting Cross-Section ===%
x_axis = -0.01:0.0005:0.01;
% y_axis = 0.5*x_axis + 0.005;
y_axis = 0.003*ones(size(x_axis));
plot3([x_axis(1) ; x_axis(end)],[y_axis(1) ; y_axis(end)],0.3*ones(2,1),'--k')
hold off

%=== Plotting Cross-Section With Confidence Intervals ===%
Coeff = zeros(3,1);
mu_aTrueInterpolatedValues = zeros(size(x_axis,2),1);
mu_aReconInterpolatedValues = zeros(size(x_axis,2),1);
PostUEInterpolatedValues = zeros(size(x_axis,2),1);

for ii = 1:size(x_axis,2);
    K = pointLocation(PrecomputedIntrplteObjectsI.DT,[x_axis(ii),y_axis(ii)]); %returns the element that the gridpoint lies in.
    Node_1 = PrecomputedIntrplteObjectsI.DT(K,1);
    Node_2 = PrecomputedIntrplteObjectsI.DT(K,2);
    Node_3 = PrecomputedIntrplteObjectsI.DT(K,3);

    %Interpolated Values for p0True
    z = [mu_aTrueIntrplte(Node_1); mu_aTrueIntrplte(Node_2); mu_aTrueIntrplte(Node_3)];
    B = cell2mat(PrecomputedIntrplteObjectsI.PrecompCoeffMatricesInv(K));
    Coeff = B*z;
    mu_aTrueInterpolatedValues(ii) = Coeff(1) + Coeff(2)*x_axis(ii) + Coeff(3)*y_axis(ii);
    
    %Interpolated Values for p0Recon
    z = [mu_aRecon(Node_1); mu_aRecon(Node_2); mu_aRecon(Node_3)];
    B = cell2mat(PrecomputedIntrplteObjectsI.PrecompCoeffMatricesInv(K));
    Coeff = B*z;
    mu_aReconInterpolatedValues(ii) = Coeff(1) + Coeff(2)*x_axis(ii) + Coeff(3)*y_axis(ii);
    
    %Interpolated Values for Uncertainty Estimates
    z = [UEDiag(Node_1); UEDiag(Node_2); UEDiag(Node_3)];
    B = cell2mat(PrecomputedIntrplteObjectsI.PrecompCoeffMatricesInv(K));
    Coeff = B*z;
    PostUEInterpolatedValues(ii) = Coeff(1) + Coeff(2)*x_axis(ii) + Coeff(3)*y_axis(ii);
end

figure
ciplot((mu_aReconInterpolatedValues-sqrt(PostUEInterpolatedValues))*10^-3,(mu_aReconInterpolatedValues+sqrt(PostUEInterpolatedValues))*10^-3,x_axis,'k');
hold on
ciplot((mu_aReconInterpolatedValues-2*sqrt(PostUEInterpolatedValues))*10^-3,(mu_aReconInterpolatedValues+2*sqrt(PostUEInterpolatedValues))*10^-3,x_axis,'k');
ciplot((mu_aReconInterpolatedValues-3*sqrt(PostUEInterpolatedValues))*10^-3,(mu_aReconInterpolatedValues+3*sqrt(PostUEInterpolatedValues))*10^-3,x_axis,'k');
plot(x_axis,mu_aTrueInterpolatedValues*10^-3,'r');
plot(x_axis,mu_aReconInterpolatedValues*10^-3,'b');
ylim([-0.1,0.3]);
keyboard



