function ApproxPosteriorCovariance = EWE_DGM2D_ApproxUncertaintyEstimates(RunOptions,FrwdFunction,SensorsI,Cov_pr,hRecon,alphaRecon,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,PrecomputedIntrplteObjectsI,pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT)

% EWE_DGM2D_ApproxUncertaintyEstimates computes the approximation posterior
% covariance for the QPAT inverse problem.
%
% Inputs:
%   RunOptions:
%   FrwdFunction - Script for the forward function of h
%   SensorsI: 1 by number of sensors array containing 3 by 1 cells 
%        - id: Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%        - xy: coordinates of the sensor
%        - l_iatsensor: 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%   Cov_pr - Prior covariance
%   hRecon - Reconstructed absorbed energy density. This is required when positivity transform is implemented, because the mapping is linearized around the MAP estimate
%   L_ET11 - Error Model for T11, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   L_ET22 - Error Model for T22, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   L_ET12 - Error Model for T12, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   L_Evx - Error Model for vx, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   L_Evy - Error Model for vy, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   Prior - Prior model statistics
%   MeshIElm: Number of elements by 3 array containing the indices of the nodes in each element
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshINodes: 2 by Number of nodes array for the coordinates of the nodes on the inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   PrecompIntrplteObjectsI - Objects that depend on the inverse Mesh nodes. May have been computed earlier and so can 
%                             be called here to avoid repeating computations. Set to 0 in function call if you
%                             want to compute new objects for a new mesh.
%   pinfoI - Invesion mesh information regarding p-refinement for p-nonconforming meshes
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI: Number of grid points in one element of the DGM inversion mesh
%   KI: Number of elements of the DGM inversion mesh
%   rhoI: Medium density of inverse mesh for computing initial condition from h
%   dt: Time step size
%   PLOT - To plot or not to plot, that is the question
%
% Outputs:
%   ApproxPosteriorCovariance - approximation to the posterior covariance
%
% Hwan Goh 3/7/2018, University of Auckland, New Zealand (Back in NZ!)

disp(' ')
disp('------------------------------------------------------')
printf(['Constructing Uncertainty Estimates ']);
disp('------------------------------------------------------')
tic
printf(['For the Case ' RunOptions.SaveFileName]);

PLOT.DGMForward = 0; %Suppress plotting of wave propagation
PLOT.AdjWaveField = 0; %Suppress plotting of adjoint wave field propagation
PLOT.DGMForwardSensorData = 0; %Suppress plotting of sensory data
NumberofSensors = size(SensorsI,2);
NumberofTimeSteps = RunOptions.NumberofTimeSteps;
Cov_prHalf = chol(Cov_pr); %Matrix square root of prior covariance
RunOptions.SaveFileNamePostCov = sprintf('ApproxPostCov-%s',RunOptions.SaveFileName);

%% =======================================================================%
%                      Halko Probabilistic Algorithm
%=========================================================================%
%=== If Positivity Constraint is Implemented ===%
if RunOptions.EWE_LS_LogExpTransform == 1
    preLogExpThRecon = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*hRecon)-1);
    pDerivhRecon = exp(Prior.LogExpTkappa*preLogExpThRecon)./(exp(Prior.LogExpTkappa*preLogExpThRecon)+1); %Derivative of positivity transform at prePosThRecon
else
    pDerivhRecon = 'Not Required';
end
if RunOptions.EWE_LS_SigmoidTransform == 1 || RunOptions.EWE_LS_ExponentialTransform == 1
    error('EWE_DGM2D_AdjStateHessianAction has not been set up for these positivity constraints')
end

%=== Stage A ===%
k = 200; %Approximate rank
W = randn(MeshIN_Nodes,2*k);
Y = zeros(MeshIN_Nodes,2*k);

for ii=1:2*k;
    disp(' ')
    disp('------------------------------------------------------')
    disp('Computing PPHoDM^3*w')
    disp('------------------------------------------------------')
    printf(['For the Case ' RunOptions.SaveFileName]);
    printf(['Current iteration number: ' num2str(ii) ' of ' num2str(2*k)]);
    TempW = W(:,ii);
    for jj=1:3 %Computing PPHoDM^3*(w)
        TempW = EWE_DGM2D_PPHoDM(RunOptions,TempW,FrwdFunction,SensorsI,Cov_prHalf,pDerivhRecon,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,MeshIN_Elm,MeshINodes,pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,PLOT);
    end
    Y(:,ii) = TempW;
end
[Q,~] = qr(Y,0); %Q is (MeshIN_Nodes x 2k). If Y is rank deficient, then Q is potentially thinner, (n x d) where d <= 2k
% clear W Y TempW

save(RunOptions.SaveFileNamePostCov,'W','Y','Q','-v7.3')

%=== Stage B ===%
PPHoDMQ = zeros(MeshIN_Nodes,2*k); 
for ii=1:2*k %Computing PPHoDM*Q
    disp(' ')
    disp('------------------------------------------------------')
    disp('Computing PPHoDM*Q')
    disp('------------------------------------------------------')
    printf(['For the Case ' RunOptions.SaveFileName]);
    printf(['Current iteration number: ' num2str(ii) ' of ' num2str(2*k)]);
    PPHoDMQ(:,ii) = EWE_DGM2D_PPHoDM(RunOptions,Q(:,ii),FrwdFunction,SensorsI,Cov_prHalf,pDerivhRecon,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,MeshIN_Elm,MeshINodes,pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,PLOT);
end
B = Q'*PPHoDMQ; %B is (2k x 2k)
[V_B,D_B] = eig(B);
% clear PPHoDMQ B

%=== Eigen Decomposition of PPHoDM ===%
V = Q*V_B; %Eigenvectors of PPHoDM
D = speye(size(D_B,1));
for ii=1:size(D_B,1)
    D(ii,ii) = D_B(ii,ii)/(D_B(ii,ii)+1);
end
% clear V_B V_D

ApproxPosteriorCovariance = Cov_pr - Cov_prHalf*V*D*V'*Cov_prHalf';

save(RunOptions.SaveFileNamePostCov,'ApproxPosteriorCovariance','W','Y','Q','B','V_B','D_B','V','D','Cov_pr','Cov_prHalf','-v7.3')
keyboard

%% =======================================================================%
%            Prior-Preconditioned Hessian of the Data Misfit
%=========================================================================%
function PPHoDM = EWE_DGM2D_PPHoDM(RunOptions,h,FrwdFunction,SensorsI,Cov_prHalf,pDerivhRecon,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,MeshIN_Elm,MeshINodes,pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,PLOT)

% EWE_DGM2D_PPHDM compute the action of the Prior-Preconditioned Hessian of
% the Data Misfit. This operation represents the action of a 
% (MeshIN_Nodes x MeshIN_Nodes) matrix.
%
% Hwan Goh 4/7/2018, University of Auckland, New Zealand (Back in NZ!)

CovHalfh = Cov_prHalf*h;
if RunOptions.EWE_LS_LogExpTransform == 1
    CovHalfh = pDerivhRecon.*CovHalfh;
end
[T11AdjTime0,T22AdjTime0,~,~,~,~,~] = EWE_DGM2D_HessianMisfit(RunOptions,CovHalfh,FrwdFunction,SensorsI,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,MeshIN_Elm,MeshINodes,pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,PLOT);
T11AdjTime0 = PrecomputedIntrplteObjectsI.InterpMatrix'*T11AdjTime0(:);
T22AdjTime0 = PrecomputedIntrplteObjectsI.InterpMatrix'*T22AdjTime0(:);
PPHoDM = Cov_prHalf'*(-pDerivhRecon.*(T11AdjTime0+T22AdjTime0));