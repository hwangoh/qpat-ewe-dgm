function [hRecon, alphaRecon, AcousticInverseItrtnInfo, Trmntn]= EWE_DGM2D_LineSearch_h(RunOptions,DataVrblsWave,FrwdFunction,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,MeshDN_Elm,MeshDNodes,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,Prior,PrecomputedIntrplteObjectsI,pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT)

% EWE_DGM2D_LineSearch_h computes the approximation of the absorbed energy
% density for the QPAT inverse problem using a given forward operator using line
% search.
%
% Inputs:
%   RunOptions:
%      InverseCrime - Whether or not we are being inverse criminals
%      EWE_LS_MaxItrtn - Max number of Iterations to be computed
%      EWE_LS_ErrTol - Termination condition, difference in error between iterations
%      EWE_LS_ExponentialTransform - Use exponential transform
%      EWE_LS_Analytical - Use analytical stepsize; elastic wave propagation is linear with respect to initial condition
%      EWE_LS_fminsearch - Use fminsearch to computer stepsize
%   DataVrblsWave:
%      T11DataTimeSteps - Time step data for T11, dimensions are NumberofSensors by NumberofTimeSteps
%      T22DataTimeSteps - Time step data for T22, dimensions are NumberofSensors by NumberofTimeSteps
%      T12DataTimeSteps - Time step data for T12, dimensions are NumberofSensors by NumberofTimeSteps
%      vxDataTimeSteps -  Time step data for vx, dimensions are NumberofSensors by NumberofTimeSteps
%      vyDataTimeSteps -  Time step data for vy, dimensions are NumberofSensors by NumberofTimeSteps
%      SensorsI: 1 by number of sensors array containing 3 by 1 cells 
%        - id: Np by 1 array containing the indices of the nodes of the element the sensor is contained in
%        - xy: coordinates of the sensor
%        - l_iatsensor: 1 by Np array representing the local basis function; [l_1(r_0,s_0),l_2(r_0,s_0),...,l_Np(r_0,s_0)] where (r_0,s_0) is such that x^k(r_0,s_0) is the coordinates of a sensor
%   FrwdFunction - Script for the forward function of h
%   L_ET11 - Error Model for T11, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   L_ET22 - Error Model for T22, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   L_ET12 - Error Model for T12, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   L_Evx - Error Model for vx, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   L_Evy - Error Model for vy, dimensions are (NumberofSensors*NumberofTimeSteps) by (NumberofSensors*NumberofTimeSteps)
%   T11ErrorMean - Mean of error model for T11, dimensions are NumberofSensors by NumberofTimeSteps
%   T22ErrorMean - Mean of error model for T22, dimensions are NumberofSensors by NumberofTimeSteps
%   T12ErrorMean - Mean of error model for T12, dimensions are NumberofSensors by NumberofTimeSteps
%   vxErrorMean - Mean of error model for vx, dimensions are NumberofSensors by NumberofTimeSteps
%   vyErrorMean - Mean of error model for vy, dimensions are NumberofSensors by NumberofTimeSteps
%   MeshDN_Elm: Number of elements in data mesh
%   MeshDNodes: 2 by Number of nodes array for the coordinates of the nodes on the data mesh
%   MeshIElm: Number of elements by 3 array containing the indices of the nodes in each element
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshINodes: 2 by Number of nodes array for the coordinates of the nodes on the inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   Prior:
%      L_pr - Regularization operator
%      Exp_h - Expected value
%      traceCov_h - Trace of covariance matrix of prior model
%   PrecompIntrplteObjectsI - Objects that depend on the inverse Mesh nodes. May have been computed earlier and so can 
%                             be called here to avoid repeating computations. Set to 0 in function call if you
%                             want to compute new objects for a new mesh.
%   pinfoI - Invesion mesh information regarding p-refinement for p-nonconforming meshes
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI: Number of grid points in one element of the DGM inversion mesh
%   KI: Number of elements of the DGM inversion mesh
%   rhoI: Medium density of inverse mesh for computing initial condition from h
%   VertexNodesGlobalIndicesFEMI: indices referring to the vertices of elements. This is for circumventing the gradient inconsistencies arising from the adjoint state method; only the gradient values at these points are considered
%   dt: Time step size
%   PLOT - To plot or not to plot, that is the question
%
% Outputs:
%   hRecon - the estimated value of the initial pressure
%   alphaRecon - the estimated hyperprior parameter
%   AcousticInverseItrtnInfo - Matrices where each column contains the Steepest-Descent approximation for one iteration. This structure contains:
%                              hReconItrtn, fhErrorSizeItrtn,hRelativeErrorItrtn, hExpErrorItrtn,DirctnItrtn
%   Trmtn - The line search iteration number at which the iterations were
%           terminated
%
% Hwan Goh, University of Auckland, New Zealand 11/09/2017
% Last Edited: 16/11/2017 - removed reliance on Globals2D and save everything into structures

%% =======================================================================%
%                   Setting Up Steepest-Descent Objects
%=========================================================================%
%=== Shortening Labels ===%
hTrue = DataVrblsWave.h;
hRecon = DataVrblsWave.hRecon;
alphaRecon = Prior.Hyperprior_alpha;
SensorsI = DataVrblsWave.SensorsI;
T11DataTimeSteps = DataVrblsWave.T11DataTimeSteps;
T22DataTimeSteps = DataVrblsWave.T22DataTimeSteps;
T12DataTimeSteps = DataVrblsWave.T12DataTimeSteps;
vxDataTimeSteps = DataVrblsWave.vxDataTimeSteps;
vyDataTimeSteps = DataVrblsWave.vyDataTimeSteps;
traceCov_h = Prior.traceCov_h;
Exp_h = Prior.Exp_h;
NumberofSensors = DataVrblsWave.NumberofSensors;
NumberofTimeSteps = RunOptions.NumberofTimeSteps;

%=== Setting User Defined Options ===%
EWE_LS_MaxItrtn=RunOptions.EWE_LS_MaxItrtn;
EWE_LS_ErrTol=RunOptions.EWE_LS_ErrTol;

%=== Interpolated Actual h Values For Computing the Relative Error ===%
if RunOptions.InverseCrime == 0;
    hTrueIntrplte = IntrplteOver2DTriangulatedMesh(MeshDN_Elm,MeshDNodes,hTrue,MeshINodes(:,1),MeshINodes(:,2),MeshIN_Nodes,0);
end

%=== Storage Vectors ===%
hReconItrtn = sparse(MeshIN_Nodes,EWE_LS_MaxItrtn); 
alphaReconItrtn = sparse(1,EWE_LS_MaxItrtn);
T11hReconItrtn = sparse(NumberofSensors*NumberofTimeSteps,EWE_LS_MaxItrtn);
T22hReconItrtn = sparse(NumberofSensors*NumberofTimeSteps,EWE_LS_MaxItrtn);
T12hReconItrtn = sparse(NumberofSensors*NumberofTimeSteps,EWE_LS_MaxItrtn);
vxhReconItrtn = sparse(NumberofSensors*NumberofTimeSteps,EWE_LS_MaxItrtn);
vyhReconItrtn = sparse(NumberofSensors*NumberofTimeSteps,EWE_LS_MaxItrtn);
fhErrorSizeItrtn = sparse(EWE_LS_MaxItrtn,1); %||f_d - f(h)||^2
hRelativeErrorItrtn = sparse(EWE_LS_MaxItrtn,1); %100*||hTrue - hRecon||/||hTrue||
hExpectedErrorItrtn = sparse(EWE_LS_MaxItrtn,1); %||hTrue - hRecon||^2/(Tr(Cov_pr)+||E[h_pr]||^2)
DirctnItrtn = sparse(MeshIN_Nodes,EWE_LS_MaxItrtn); 

%=== Setting Up Plotting ===%
PLOT.TRIFEM=delaunay(MeshINodes(:,1),MeshINodes(:,2));
PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshI; %This is for plotting forward elastic wave propagation on inversion mesh. FwrdFunction is called which is EWE_DGM2D_LSExpRK4. However, that function calls PLOT.TRI_DGMMeshD when plotting, so here we replace it with DGMMeshI but re-use the same name

%=== Iterations ===%
for ii=1:EWE_LS_MaxItrtn;
    disp(' ')
    disp('------------------------------------------------------')
    printf(['Computing line search iteration number ' num2str(ii) ' of ' num2str(EWE_LS_MaxItrtn)]); 
    disp('------------------------------------------------------')
    tic
    printf(['For the Case ' RunOptions.SaveFileName]);
    %% ===================================================================%
    %                   Updating Using Current hRecon
    %=====================================================================%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plotting Current hRecon at Each Step %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(PLOT.Figure_CurrenthRecon)
    trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),hRecon);
    if PLOT.DGMPlotBirdsEyeView == 1;
        view(2)
    end
    if PLOT.DGMPlotUsezLim == 1;
        zlim(PLOT.AbsorbedEnergyDensityzAxis);
    end
    shading interp %thanks Ru!
    caxis([0 200])
    colormap(jet(256))
    title(PLOT.Figure_CurrenthRecon_Title,'FontWeight','bold')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Current f(hRecon) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    hReconDGM = IntrplteOver2DTriangulatedMesh(MeshIN_Elm,MeshINodes,hRecon,xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
    IniCond = RunOptions.LinearHeatExpansionCoeff*reshape(hReconDGM,NpI,KI)./(rhoI*RunOptions.SpecificHeatCoeff);
    [T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps] = FrwdFunction(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT);
  
    %=== Errors ===%
    [T11InterpolatedTimeSteps,T22InterpolatedTimeSteps,T12InterpolatedTimeSteps,vxInterpolatedTimeSteps,vyInterpolatedTimeSteps,Error_qh,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps] =...
        EWE_DGM2D_ObjectiveFunctionalErrorTerm(RunOptions,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                               T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,...
                                               T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps,...
                                               T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                               SensorsI);
    if RunOptions.InverseCrime ~= 1;
        hError = hTrueIntrplte - hRecon;
        hRelError = 100*(sqrt(dot(hError,hError)))/(sqrt(dot(hTrueIntrplte,hTrueIntrplte)));
    else
        hError = hTrue(:) - hRecon;
        hRelError = 100*(sqrt(dot(hError,hError)))/(sqrt(dot(hTrue(:),hTrue(:))));
    end
    hExpError = dot(hError,hError)/(traceCov_h + dot(Exp_h,Exp_h));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Storing the Current Values %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hReconItrtn(:,ii) = hRecon;
    alphaReconItrtn(ii) = alphaRecon;
    T11hReconItrtn(:,ii) = T11InterpolatedTimeSteps(:);
    T22hReconItrtn(:,ii) = T22InterpolatedTimeSteps(:);
    T12hReconItrtn(:,ii) = T12InterpolatedTimeSteps(:);
    vxhReconItrtn(:,ii) = vxInterpolatedTimeSteps(:);
    vyhReconItrtn(:,ii) = vyInterpolatedTimeSteps(:);
    fhErrorSizeItrtn(ii) = Error_qh;
    hRelativeErrorItrtn(ii) = hRelError;
    hExpectedErrorItrtn(ii) = hExpError;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Displaying the Current Values %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Relative Error of Current Estimate = ' num2str(hRelativeErrorItrtn(ii))]);
    disp(['Expected Error of Current Estimate = ' num2str(hExpectedErrorItrtn(ii))]);
    disp(' ')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Termination Conditions %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ii>1 && (abs(hRelativeErrorItrtn(ii)^2 - hRelativeErrorItrtn(ii-1)^2)<EWE_LS_ErrTol)
        printf('Iterations have reached tolerance in difference in errors');
        break
    end
    if ii>1 && ii == EWE_LS_MaxItrtn
        printf('Max number of iterations reached');
        break
    end

    %% ===================================================================%
    %                       Computing New hRecon
    %=====================================================================%
    %=== Adjoint-State Method ===%
    [Gradient,Gradient_alpha,T11AdjTime0,T22AdjTime0] = EWE_DGM2D_AdjointStateMethod_Gradient_h(RunOptions,SensorsI,alphaRecon,...
                                                        L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                        L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps,...
                                                        hRecon,Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,...
                                                        pinfoI,xI,yI,NpI,KI,PrecomputedIntrplteObjectsI.InterpMatrix,dt,PLOT,ii);
    if RunOptions.EWE_DGM2D_UseSteepestDescent == 1;
        Dirctn = -Gradient/norm(Gradient,2);
        if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 1
            Dirctn_palpha = -[Gradient;Gradient_alpha]/norm([Gradient;Gradient_alpha],2);
            Dirctn = Dirctn_palpha(1:end-1);
            Dirctn_alpha = Dirctn_palpha(end);
        end
    end
    if RunOptions.EWE_DGM2D_UseNewton == 1 || RunOptions.EWE_DGM2D_UseGaussNewton == 1;
        if ii == 1;
            Gradient_InitialGuess = Gradient; %To be used as tolerance for conjugate gradient method
        end
        Dirctn = EWE_DGM2D_AdjointStateMethod_Newton_h(RunOptions,SensorsI,FrwdFunction,hRecon,alphaRecon,Gradient,Gradient_alpha,Gradient_InitialGuess,...
                                                             L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                             T11AdjTime0,T22AdjTime0,...
                                                             Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,...
                                                             pinfoI,PrecomputedIntrplteObjectsI,xI,yI,NpI,KI,rhoI,dt,ii,hRelError,PLOT);
        Dirctn = Dirctn/norm(Dirctn,2);
        if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 1
            Dirctn_alpha = Dirctn(end);
            Dirctn = Dirctn(1:end-1);
        end
    end
    DirctnItrtn(:,ii) = Dirctn;
    
    %=== Plotting Search Direction ===%
    if PLOT.AdjStSearchDirection == 1;
        figure(PLOT.Figure_AdjStSearchDirection)
        trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),-Dirctn);
        if PLOT.DGMPlotBirdsEyeView == 1;
            view(2)
        end
        shading interp %thanks Ru!
        colormap(jet(256))
        title(PLOT.Figure_AdjStSearchDirection_Title,'FontWeight','bold')
    end

    %=== Computing Step Size ===%
    if RunOptions.EWE_LS_Analytical == 1
        disp('Analytically computing step size to compute next estimate');
        if RunOptions.EWE_LS_ExponentialTransform == 1 || RunOptions.EWE_LS_SigmoidTransform == 1
            error('Acoustic forward problem is no longer linear with transforms; cannot use analytical stepsize')
        end
        if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 1
            error('Inverse problem is no longer linear when using a hyperprior')
        end
        StepSize = EWE_DGM2D_AnalyticalStepSize(RunOptions,FrwdFunction,Dirctn,hRecon,...
                                                  L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                  L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps,...
                                                  Prior,MeshIN_Elm,MeshINodes,PrecomputedIntrplteObjectsI,SensorsI,...
                                                  pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT);
    end
    if RunOptions.EWE_LS_LogExpTransform == 1 %Computing pre-transformed parameter
        hRecon = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*hRecon)-1);
    end
    if RunOptions.EWE_LS_SigmoidTransform == 1 %Computing pre-transformed parameter
        hRecon = (atanh((hRecon - Prior.SigTalpha)/Prior.SigTbeta))/Prior.SigTkappa + Prior.SigTdelta;
    end
    if RunOptions.EWE_LS_ExponentialTransform == 1 %Computing pre-transformed parameter
        hRecon = log(hRecon);
    end
    if RunOptions.EWE_LS_fminsearch == 1
        disp('Using Line Search with fminsearch to compute next estimate');
        options = optimset('Display','iter','PlotFcns',@optimplotfval);
        if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 0
            [StepSize,~] = fminsearch(@(step) LSFunctional(RunOptions,FrwdFunction,hRecon + step*Dirctn,'NoHyperprior',...
                                                           L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                           T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,...
                                                           T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                                           Prior,MeshIN_Elm,MeshINodes,PrecomputedIntrplteObjectsI,SensorsI,...
                                                           pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT,ii,hRelError),0,options);
        end
        if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 1
            [StepSize,~] = fminsearch(@(step) LSFunctional(RunOptions,FrwdFunction,hRecon + step*Dirctn,alphaRecon + step*Dirctn_alpha,...
                                                           L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                           T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,...
                                                           T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                                           Prior,MeshIN_Elm,MeshINodes,PrecomputedIntrplteObjectsI,SensorsI,...
                                                           pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT,ii,hRelError),0,options);
        end
    end
    printf(['Step Size: ' num2str(StepSize)])
    
    %=== New hRecon Value ===%
    hRecon = hRecon + StepSize*Dirctn;
    if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 1    
        alphaRecon = alphaRecon + StepSize*Dirctn_alpha; 
    end
    if RunOptions.EWE_LS_LogExpTransform == 1
        hRecon = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*hRecon)+1);
    end
    if RunOptions.EWE_LS_SigmoidTransform == 1
        hRecon = Prior.SigTalpha + Prior.SigTbeta*tanh(Prior.SigTkappa*(hRecon - Prior.SigTdelta));
    end
    if RunOptions.EWE_LS_ExponentialTransform == 1
        hRecon = exp(hRecon);
    end
    toc
end

%=== If Steepest-Descent Iterations Were Terminated ===%
Trmntn = ii;
printf(['Terminated after ' num2str(Trmntn) ' iterations.']);

if (Trmntn==EWE_LS_MaxItrtn)
    printf('Maximum iterations reached.');
end

hRecon = hReconItrtn(:,Trmntn);
alphaRecon = alphaReconItrtn(Trmntn);

%% =======================================================================%
%                                Outputs
%=========================================================================%
AcousticInverseItrtnInfo.hRecon = hReconItrtn;
AcousticInverseItrtnInfo.alphaRecon = alphaReconItrtn;
AcousticInverseItrtnInfo.T11hRecon = T11hReconItrtn;
AcousticInverseItrtnInfo.T22hRecon = T22hReconItrtn;
AcousticInverseItrtnInfo.T12hRecon = T12hReconItrtn;
AcousticInverseItrtnInfo.vxhRecon = vxhReconItrtn;
AcousticInverseItrtnInfo.vyhRecon = vyhReconItrtn;
AcousticInverseItrtnInfo.fhErrorSize = fhErrorSizeItrtn;
AcousticInverseItrtnInfo.hRelativeError = hRelativeErrorItrtn;
AcousticInverseItrtnInfo.hExpectedError = hExpectedErrorItrtn;
AcousticInverseItrtnInfo.Dirctn = DirctnItrtn;

SaveFileNameReconstructions = sprintf('Reconstructions-%s',RunOptions.SaveFileName)
save(SaveFileNameReconstructions,'DataVrblsWave','Prior','hTrueIntrplte','AcousticInverseItrtnInfo','PLOT','MeshINodes','Trmntn','-v7.3')
keyboard
end

%% =======================================================================%
%                        Line Search Functional
%=========================================================================%
function LS_Error = LSFunctional(RunOptions,FrwdFunction,NewhRecon,NewalphaRecon,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,Prior,MeshIN_Elm,MeshINodes,PrecomputedIntrplteObjectsI,SensorsI,pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT,CurrentIteration,hRelError)

% LSFunctional computes ||Data - Newf(h)|| using the new h
% value where Newh = h + StepSize*Dirctn. Then LSFunctional can
% be considered to be a function of StepSize and thus used with
% fminsearch for linesearch purposes.
%
% Hwan Goh 18/9/2017, University of Auckland, New Zealand

PLOT.DGMForward = 0; %Suppress plotting of wave propagation
PLOT.DGMForwardSensorData = 0; %Suppress plotting of sensory data

disp(' ')
disp('------------------------------------------------------')
disp('Using fminsearch ')
disp('------------------------------------------------------')
printf(['For the Case ' RunOptions.SaveFileName]);
printf(['Current line search iteration number: ' num2str(CurrentIteration)]);
printf(['Current relative error: ' num2str(hRelError)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computing and Interpolating f(hRecon + alpha*Dirctn) %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Computing f(hRecon + alpha*Dirctn) ===%
if RunOptions.EWE_LS_LogExpTransform == 1
    NewhRecon = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*NewhRecon)+1);
end
if RunOptions.EWE_LS_SigmoidTransform == 1
    NewhRecon = Prior.SigTalpha + Prior.SigTbeta*tanh(Prior.SigTkappa*(NewhRecon - Prior.SigTdelta));
end
if RunOptions.EWE_LS_ExponentialTransform == 1
    NewhRecon = exp(NewhRecon);
end
NewhReconDGM = IntrplteOver2DTriangulatedMesh(MeshIN_Elm,MeshINodes,NewhRecon,xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
IniCond = reshape(NewhReconDGM,NpI,KI)./(rhoI*RunOptions.SpecificHeatCoeff);
PLOT.DGMForward = 0;
[T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps] = FrwdFunction(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT); 
[~,~,~,~,~,~,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps] = EWE_DGM2D_ObjectiveFunctionalErrorTerm(RunOptions,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,SensorsI);
if RunOptions.EWE_LS_RemoveTimesSteps ~= 0;
    L_ET11ErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
    L_ET22ErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
    L_ET12ErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
    L_EvxErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
    L_EvyErrorTimeSteps(:,1:RunOptions.EWE_LS_RemoveTimesSteps) = 0;
end
L_EError = [L_ET11ErrorTimeSteps;L_ET22ErrorTimeSteps;L_ET12ErrorTimeSteps;L_EvxErrorTimeSteps;L_EvyErrorTimeSteps];

%%%%%%%%%%%%%%%%%%%%%%
%%% Relative Error %%% 
%%%%%%%%%%%%%%%%%%%%%%
if RunOptions.EWE_LS_LogExpTransform == 1 %Computing pre-transformed parameter
    NewhRecon = (1/Prior.LogExpTkappa)*log(exp(Prior.LogExpTkappa*NewhRecon)-1);
end
if RunOptions.EWE_LS_SigmoidTransform == 1 %Computing pre-transformed parameter
    NewhRecon = (atanh((NewhRecon - Prior.SigTalpha)/Prior.SigTbeta))/Prior.SigTkappa + Prior.SigTdelta; %inverse of rior.SigTalpha + Prior.SigTbeta*tanh(Prior.SigTkappa*(hRecon - Prior.SigTalpha))
end
if RunOptions.EWE_LS_ExponentialTransform == 1 %Computing pre-transformed parameter
    NewhRecon = log(NewhRecon);
end
if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 0
    LS_Error = (1/2)*norm(L_EError,2) + (1/2)*norm(Prior.L_pr*(Prior.Exp_h - NewhRecon),2);
end
if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior == 1
    LS_Error = (1/2)*norm(L_EError,2) + ((NewalphaRecon^2)/2)*norm(Prior.L_pr*(Prior.Exp_h - NewhRecon),2) + Prior.Hyperprior_Std_alpha^(-2)*(Prior.Hyperprior_Exp_alpha - NewalphaRecon)^2 - size(MeshINodes,1)*log(abs(NewalphaRecon));
    printf(['Current alpha value: ' num2str(NewalphaRecon)]);
end
printf(['Current functional error: ' num2str(LS_Error)]);
end


