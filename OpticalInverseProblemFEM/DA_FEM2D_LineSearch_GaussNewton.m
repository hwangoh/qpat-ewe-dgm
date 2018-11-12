function [mu_aRecon, OpticalInverseItrtnInfo, Trmntn]= DA_FEM2D_LineSearch_GaussNewton(RunOptions,DataVrblsOptical,FrwdFunction,ConstructJacQPATFrwd,RegData,L_E,Prior,PrmtrsD,PrmtrsI,PrmtrsPrp,MeshD,MeshI,PLOT)

% DA_FEM2D_LineSearch_GaussNewton computes the Gauss-Newton approximation of the absorption
% coefficient for the QPAT inverse problem using a given forward operator.
%
% Inputs:
%   RunOptions:
%      DA_GN_MaxItrtn - Max number of Iterations to be computed
%      DA_GN_ErrTol - Termination condition, difference in error between iterations
%      DA_GN_LS_fminsearch - Use fminsearch for line search
%      DA_GN_LS_PaulsFunction - Use Paul's line search function
%      DA_GN_LS_StepSize - Initial stepsize for Paul's line search function
%      DA_GN_LS_MaxStepSize - Max step size permitted for Paul's line search function
%      DA_GN_LS_MaxItrtn - Max number of iterations for Paul's line search function
%      DA_GN_LS_PostivityConstraint - Whether positivity constraints are implemented for Paul's line search function
%      DA_GN_LS_mu_aMin - Minimum value permitted mu_a for Paul's line search function
%   DataVrblsOptical: 
%      Gruneisen - Gruneisen Parameter
%   FrwdFunction - Script for the forward function of mu_a
%   ConstructJacFrwd - Script for the Jacobian of FrwdFunction
%   RegData - Regularized data vector: [L_1*Data(:); Prior.L_2*Prior.mu_abar]
%   L_E - Error model
%   Prior:
%      L_pr - Regularization operator
%      Exp_mu_a - Prior
%   PrmtrsD:
%      mu_a - the actual absorption coefficient
%   PrmtrsI:
%      mu_a - the absorption coefficient
%      kappa - the diffusion coefficient
%   PrmtrsPrp:
%      g - mean of the cosine of the scattering angle, used for computing the diffusion coefficient
%   MeshD: Data mesh properties, only used for interpolation of actual parameter values to compute the relative error
%   MeshI: 
%      N_Nodes - Number of nodes
%   PLOT - To plot or not to plot, that is the question
%
% Outputs:
%   mu_aRecon - the estimated value of the absorption coefficient
%   OpticalInverseItrtnInfo - Matrices where each column contains the Steepest-Descent approximation for one iteration. This structure contains:
%                              mu_aReconItrtn, fmu_aErrorSizeItrtn, mu_aRelativeErrorItrtn, mu_aExpErrorItrtn
%   Trmtn - The Steepest-Descent iteration number at which the iterations were
%           terminated
%
% Hwan Goh, University of Auckland, New Zealand 22/12/2015

%% =======================================================================%
%                   Setting Up Gauss-Newton Objects
%=========================================================================%
%=== Shortening Labels ===%
mu_aRecon = PrmtrsI.mu_a;
N_Nodes = MeshI.N_Nodes;
L_pr = Prior.L_pr;
traceCov_mu_a = Prior.traceCov_mu_a;
Exp_mu_a = Prior.Exp_mu_a;

%=== Setting User Defined Options ===%
DA_GN_MaxItrtn=RunOptions.DA_GN_MaxItrtn;
DA_GN_ErrTol=RunOptions.DA_GN_ErrTol;
DA_GN_LS_StepSize=RunOptions.DA_GN_LS_StepSize;

%=== Interpolated Actual mu_a Values For Computing the Relative Error ===%
if RunOptions.InverseCrime == 0;
    PrmtrsD.mu_aIntrplte = IntrplteOver2DTriangulatedMesh(MeshD.N_Elm,MeshD.Nodes,PrmtrsD.mu_a,MeshI.Nodes(:,1),MeshI.Nodes(:,2),MeshI.N_Nodes,0);
end

%=== Storage Vectors ===%
mu_aReconItrtn = sparse(N_Nodes,DA_GN_MaxItrtn); 
fmu_aErrorSizeItrtn = sparse(DA_GN_MaxItrtn,1); %||f_d - f(mu_a)||^2
mu_aRelativeErrorItrtn = sparse(DA_GN_MaxItrtn,1); %100*||mu_aTrue - mu_aRecon||/||mu_aTrue||
mu_aExpectedErrorItrtn = sparse(DA_GN_MaxItrtn,1); %||mu_aTrue - mu_aRecon||^2/(Tr(Cov_mu_a)+||E[mu_a_pr]||^2)
PrmtrsI.kappa = sparse(N_Nodes,1);

%=========================================================================%
%                         Gauss-Newton Iterations
%=========================================================================%
% if RunOptions.InvMethodDGMOneStage == 1 || RunOptions.InvMethodDGMTwoStage == 1;
%     %JacIntrplte = ConstructJacIntrplte(MeshI,varargin{:}); %Constant with respect to the argument, therefore only depends on the coordinates of the nodes for both meshes
%     load IntrplteJac-7thOrder3BndElm288
% else
%     JacIntrplte = 'Not Required';
% end

for ii=1:DA_GN_MaxItrtn;
    printf(['\nComputing Gauss-Newton iteration number ' num2str(ii) ' of ' num2str(DA_GN_MaxItrtn)]); 
    tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plotting Current Reconstructed mu_a at Each Step %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(PLOT.Figure_Currentmu_aRecon)
    PLOT.TRI_MeshI=delaunay(MeshI.Nodes(:,1),MeshI.Nodes(:,2));
    trisurf(PLOT.TRI_MeshI,MeshI.Nodes(:,1),MeshI.Nodes(:,2),mu_aRecon*10^-3);
    view(2)
    shading interp %thanks Ru!
%     colorbar
    caxis([0,0.25])
    colormap(jet(256))
    axis 'image'
%     title(PLOT.Figure_Currentmu_aRecon_Title,'FontWeight','bold')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Current f(mu_a) and Regularized f(mu_a) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ii==1;
        [fmu_a,PrmtrsI,phi,A,R,b_pattern,~] = FrwdFunction(RunOptions,DataVrblsOptical,MeshI,mu_aRecon,PrmtrsI.mu_s,PrmtrsPrp,'To be computed','To be computed');
    else
        [fmu_a,PrmtrsI,phi,A,R,b_pattern,~] = FrwdFunction(RunOptions,DataVrblsOptical,MeshI,mu_aRecon,PrmtrsI.mu_s,PrmtrsPrp,R,b_pattern);
    end
    Regfmu_a = [L_E*fmu_a;L_pr*mu_aRecon];
        
    %%%%%%%%%%%%%
    %%% Error %%%
    %%%%%%%%%%%%%
    error = RegData - Regfmu_a;
    erSize = norm(error);
    if RunOptions.InverseCrime == 0;
        mu_aError = PrmtrsD.mu_aIntrplte - mu_aRecon;
        mu_aRelError = 100*(sqrt(dot(mu_aError,mu_aError)))/(sqrt(dot(PrmtrsD.mu_aIntrplte,PrmtrsD.mu_aIntrplte)));
    else
        mu_aError = PrmtrsD.mu_a - mu_aRecon;
        mu_aRelError = 100*(sqrt(dot(mu_aError,mu_aError)))/(sqrt(dot(PrmtrsD.mu_a,PrmtrsD.mu_a)));
    end
    mu_aExpError = dot(mu_aError,mu_aError)/(traceCov_mu_a + dot(Exp_mu_a,Exp_mu_a));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Storing the Current Values %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu_aReconItrtn(:,ii) = mu_aRecon;
    fmu_aErrorSizeItrtn(ii) = erSize;
    mu_aRelativeErrorItrtn(ii) = mu_aRelError;
    mu_aExpectedErrorItrtn(ii) = mu_aExpError;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Displaying the Current Values %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Relative Error of Current Estimate = ' num2str(mu_aRelativeErrorItrtn(ii))]);
    disp(['Expected Error of Current Estimate = ' num2str(mu_aExpectedErrorItrtn(ii))]);
    disp(' ')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Termination Conditions %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ii>1 && (abs(fmu_aErrorSizeItrtn(ii)^2 - fmu_aErrorSizeItrtn(ii-1)^2)<DA_GN_ErrTol)
        printf('Iterations have reached tolerance in difference in errors');
        break
    end
    
    %% ===================================================================%
    %                       Computing New mu_aRecon
    %=====================================================================%
    %=== Jacobian ===%
    JacODTFrwd = ConstructJacQPATFrwd(MeshI,PrmtrsI,phi,A);
    RegJacForward = [L_E*JacODTFrwd;L_pr];
    Dirctn = (RegJacForward'*RegJacForward)\(RegJacForward'*error);

    %=== New mu_a Value ===%
    if RunOptions.DA_GN_LS_fminsearch == 1
        disp('Using Line Search with fminsearch to compute next estimate');
        options = optimset('Display','iter','PlotFcns',@optimplotfval);
        [DA_GN_LS_StepSize,~] = fminsearch(@(alpha) GNLSFunctional(RunOptions,DataVrblsOptical,FrwdFunction,RegData,MeshI,mu_aRecon+alpha*Dirctn,...
                                                                   PrmtrsI.mu_s,PrmtrsPrp,R,b_pattern,L_E,L_pr),0,options);
    end
    if RunOptions.DA_GN_LS_PaulsFunction == 1;
        disp('Using Line Search with Pauls Line Search Function to compute next estimate');
        [DA_GN_LS_StepSize,~] = DA_FEM2D_GN_PaulsLinesearch(RunOptions,DataVrblsOptical,FrwdFunction,RegData,error,Dirctn,MeshI,PrmtrsI,PrmtrsPrp,R,b_pattern,L_E,L_pr,ii);
    end
    mu_aRecon = mu_aRecon + DA_GN_LS_StepSize*Dirctn;
    toc
end

%=== If Gauss-Newton Iterations Were Terminated ===%
Trmntn = ii;
printf(['Terminated after ' num2str(Trmntn) ' iterations.']);

if (Trmntn==DA_GN_MaxItrtn)
    printf('Maximum iterations reached.');
end

mu_aRecon = mu_aReconItrtn(:,Trmntn);

%=========================================================================%
%                                Outputs
%=========================================================================%
OpticalInverseItrtnInfo.mu_aRecon = mu_aReconItrtn;
OpticalInverseItrtnInfo.fmu_aErrorSize = fmu_aErrorSizeItrtn;
OpticalInverseItrtnInfo.mu_aRelativeError = mu_aRelativeErrorItrtn;
OpticalInverseItrtnInfo.mu_aExpectedError = mu_aExpectedErrorItrtn;

SaveFileNameReconstructions = sprintf('Reconstructions-Optical-%s',RunOptions.SaveFileName)
mu_aIntrplte = PrmtrsD.mu_aIntrplte;
mu_aRecon = full(mu_aRecon);
save(SaveFileNameReconstructions,'DataVrblsOptical','Prior','mu_aRecon','mu_aIntrplte','OpticalInverseItrtnInfo','PLOT','MeshI','Trmntn','-v7.3')
end

%=========================================================================%
%                        Line Search Functional
%=========================================================================%
function GNLSError = GNLSFunctional(RunOptions,DataVrblsOptical,FrwdFunction,RegData,MeshI,Newmu_a,mu_s,PrmtrsPrp,R,b_pattern,L_E,L_pr)

% GNLSFunctional computes ||RegData - RegfNewmu_a|| using the new mu_a
% value where Newmu_a = mu_a + GNStepSize*Dirctn. Then GNLSFunctional can
% be considered to be a function of GNStepSize and thus used with
% fminsearch for linesearch purposes.
%
% Hwan Goh 29/12/2015, University of Auckland, New Zealand
[fmu_a,~,~,~,~,~] = FrwdFunction(RunOptions,DataVrblsOptical,MeshI,Newmu_a,mu_s,PrmtrsPrp,R,b_pattern); 
RegfNewmu_a = [L_E*fmu_a;L_pr*Newmu_a];
GNLSError = norm(RegData - RegfNewmu_a);

end

