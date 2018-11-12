function [cpRecon, csRecon, AcousticInverseItrtnInfo, Trmntn]= EWE_DGM2D_LineSearch_SteepestDescent_c(RunOptions,DataVrblsWave,FrwdFunction,p0DGM,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,MeshDN_Elm,MeshDNodes,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,Prior,MassMatrixFEM,PrecomputedIntrplteObjectsI,pinfoI,xI,yI,NpI,KI,rhoI,VertexNodesGlobalIndicesFEMI,dt,PLOT)

% EWE_DGM2D_LineSearch_SteepestDescent_c computes the approximation of the
% wave speeds for the QPAT inverse problem using a given forward operator using line search with the
% steepest-descent direction.
%
% Inputs:
%   RunOptions:
%      InverseCrime - Whether or not we are being inverse criminals
%      EWE_SD_MaxItrtn - Max number of Iterations to be computed
%      EWE_SD_ErrTol - Termination condition, difference in error between iterations
%      EWE_SD_ExponentialTransform - Use exponential transform
%      EWE_SD_LS_Analytical - Use analytical stepsize; elastic wave propagation is linear with respect to initial condition
%      EWE_SD_LS_fminsearch - Use fminsearch to computer stepsize
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
%   FrwdFunction - Script for the forward function of p0
%   MeshDN_Elm: Number of elements in data mesh
%   MeshDNodes: 2 by Number of nodes array for the coordinates of the nodes on the data mesh
%   MeshIElm: Number of elements by 3 array containing the indices of the nodes in each element
%   MeshIN_Elm: Number of elements in inversion FEM mesh
%   MeshINodes: 2 by Number of nodes array for the coordinates of the nodes on the inversion FEM mesh
%   MeshIN_Nodes: Number of nodes on the inversion FEM mesh
%   Prior:
%      L_pr - Regularization operator
%      Exp_p0 - Expected value
%      traceCov_p0 - Trace of covariance matrix of prior model
%   MassMatrixFEM - This is used in computing gradient in adjoint state method. This does not vary and so we precomputed it for use here
%   PrecompIntrplteObjectsI - Objects that depend on the inverse Mesh nodes. May have been computed earlier and so can 
%                             be called here to avoid repeating computations. Set to 0 in function call if you
%                             want to compute new objects for a new mesh.
%   pinfoI - Invesion mesh information regarding p-refinement for p-nonconforming meshes
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI: Number of grid points in one element of the DGM inversion mesh
%   KI: Number of elements of the DGM inversion mesh
%   rhoI: Medium density of inverse mesh for computing initial condition from p0
%   VertexNodesGlobalIndicesFEMI: indices referring to the vertices of elements. This is for circumventing the gradient inconsistencies arising from the adjoint state method; only the gradient values at these points are considered
%   dt: Time step size
%   PLOT - To plot or not to plot, that is the question
%
% Outputs:
%   cpRecon - the estimated value of the wave speed cp
%   csRecon - the estimated value of the wave speed cs
%   AcousticInverseItrtnInfo - Matrices where each column contains the Steepest-Descent approximation for one iteration.
%   Trmtn - The Steepest-Descent iteration number at which the iterations were
%           terminated
%
% Hwan Goh, University of Auckland, New Zealand 23/04/2018

%% =======================================================================%
%                   Setting Up Steepest-Descent Objects
%=========================================================================%
%=== Shortening Labels ===%
cpTrue = pinfoI(2).cpM;
csTrue = pinfoI(2).csM;
cpRecon = DataVrblsWave.cpRecon;
csRecon = DataVrblsWave.csRecon;
SensorsI = DataVrblsWave.SensorsI;
T11DataTimeSteps = DataVrblsWave.T11DataTimeSteps;
T22DataTimeSteps = DataVrblsWave.T22DataTimeSteps;
T12DataTimeSteps = DataVrblsWave.T12DataTimeSteps;
vxDataTimeSteps = DataVrblsWave.vxDataTimeSteps;
vyDataTimeSteps = DataVrblsWave.vyDataTimeSteps;
Exp_cp = Prior.Exp_cp;
Exp_cs = Prior.Exp_cs;
NumberofSensors = DataVrblsWave.NumberofSensors;
NumberofTimeSteps = RunOptions.NumberofTimeSteps;

%=== Setting User Defined Options ===%
EWE_SD_MaxItrtn=RunOptions.EWE_SD_MaxItrtn;
EWE_SD_ErrTol=RunOptions.EWE_SD_ErrTol;

%=== Storage Vectors ===%
cpReconItrtn = sparse(NpI*KI,EWE_SD_MaxItrtn); 
c_sReconItrtn = sparse(NpI*KI,EWE_SD_MaxItrtn); 
fcErrorSizeItrtn = sparse(EWE_SD_MaxItrtn,1); %||f_d - f(c)||^2
cpRelativeErrorItrtn = sparse(EWE_SD_MaxItrtn,1); %100*||cpTrue - cpRecon||/||cpTrue||
cpsRelativeErrorItrtn = sparse(EWE_SD_MaxItrtn,1); %100*||csTrue - csRecon||/||csTrue||
DirctnItrtn = sparse(NpI*KI,EWE_SD_MaxItrtn); 

%=== Setting Up Plotting ===%
PLOT.TRIFEM=delaunay(MeshINodes(:,1),MeshINodes(:,2));
PLOT.TRI_DGMMeshD = PLOT.TRI_DGMMeshI; %This is for plotting forward elastic wave propagation on inversion mesh. FwrdFunction is called which is EWE_DGM2D_LSExpRK4. However, that function calls PLOT.TRI_DGMMeshD when plotting, so here we replace it with DGMMeshI but re-use the same name

%=== Iterations ===%
for ii=1:EWE_SD_MaxItrtn;
    disp(' ')
    disp('------------------------------------------------------')
    printf(['Computing Steepest-Descent iteration number ' num2str(ii) ' of ' num2str(EWE_SD_MaxItrtn)]); 
    disp('------------------------------------------------------')
    tic
    printf(['For the Case ' RunOptions.SaveFileName]);
    %% ===================================================================%
    %            Updating Using Current cpRecon and csRecon
    %=====================================================================%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plotting Current cpRecon and csRecon at Each Step %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(PLOT.Figure_CurrentcpRecon)
    trisurf(PLOT.TRI_DGMMeshI,xI,yI,cpRecon);
    if PLOT.DGMPlotBirdsEyeView == 1;
        view(2)
    end
    shading interp %thanks Ru!
    colorbar
    colormap(jet(256))
    title(PLOT.Figure_CurrentcpRecon_Title,'FontWeight','bold')
    
    figure(PLOT.Figure_CurrentcsRecon)
    trisurf(PLOT.TRI_DGMMeshI,xI,yI,csRecon);
    if PLOT.DGMPlotBirdsEyeView == 1;
        view(2)
    end
    shading interp %thanks Ru!
    colorbar
    colormap(jet(256))
    title(PLOT.Figure_CurrentcsRecon_Title,'FontWeight','bold')
    
    %%%%%%%%%%%%%%%%%%%%
    %%% Current f(c) %%%
    %%%%%%%%%%%%%%%%%%%%
    IniCond = reshape(p0DGM,NpI,KI)./(rhoI*RunOptions.SpecificHeatCoeff);
    pinfoI(2).cpM = cpRecon;
    pinfoI(2).cpM = csRecon;
    if RunOptions.TimeLSERK4 == 1;
        [T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps] = FrwdFunction(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT);
    end    
    %=== Errors ===%
    [T11InterpolatedTimeSteps,T22InterpolatedTimeSteps,T12InterpolatedTimeSteps,vxInterpolatedTimeSteps,vyInterpolatedTimeSteps,Error_qc,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps] =...
        EWE_DGM2D_ObjectiveFunctionalErrorTerm(RunOptions,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                               T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,...
                                               T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps,...
                                               T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                               SensorsI);

    cpError = cpTrue(:) - cpRecon;
    csError = csTrue(:) - csRecon;
    cpRelError = 100*(sqrt(dot(cpError,cpError)))/(sqrt(dot(cpTrue(:),cpTrue(:))));
    csRelError = 100*(sqrt(dot(csError,csError)))/(sqrt(dot(csTrue(:),csTrue(:))));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Storing the Current Values %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cpReconItrtn(:,ii) = cpRecon;
    csReconItrtn(:,ii) = csRecon;
    T11cReconItrtn(:,ii) = T11InterpolatedTimeSteps(:);
    T22cReconItrtn(:,ii) = T22InterpolatedTimeSteps(:);
    T12cReconItrtn(:,ii) = T12InterpolatedTimeSteps(:);
    vxcReconItrtn(:,ii) = vxInterpolatedTimeSteps(:);
    vycReconItrtn(:,ii) = vyInterpolatedTimeSteps(:);
    fcErrorSizeItrtn(ii) = Error_qc;
    cpRelativeErrorItrtn(ii) = cpRelError;
    csRelativeErrorItrtn(ii) = csRelError;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Displaying the Current Values %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Relative Error of Current Estimate of cp = ' num2str(cpRelativeErrorItrtn(ii))]);
    disp(['Relative Error of Current Estimate of cs = ' num2str(csRelativeErrorItrtn(ii))]);
    disp(' ')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Termination Conditions %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ii>1 && (abs(cpRelativeErrorItrtn(ii)^2 - cpRelativeErrorItrtn(ii-1)^2)<EWE_SD_ErrTol)
        printf('Iterations have reached tolerance in difference in errors');
        break
    end
    if ii>1 && (abs(csRelativeErrorItrtn(ii)^2 - csRelativeErrorItrtn(ii-1)^2)<EWE_SD_ErrTol)
        printf('Iterations have reached tolerance in difference in errors');
        break
    end
%     if ii>1 && SDStepSize<0
%         printf('Stuck going back and forth, step size is negative');
%         break
%     end

    %% ===================================================================%
    %                 Computing New cpRecon and csRecon
    %=====================================================================%
    %=== Adjoint-State Method ===%
    gradient = EWE_DGM2D_AdjointStateMethod_c(RunOptions,SensorsI,...
                                              L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                              L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps,...
                                              cpRecon,csRecon,Prior,MeshIElements,MeshIN_Elm,MeshINodes,MeshIN_Nodes,pinfoI,xI,yI,NpI,KI,VertexNodesGlobalIndicesFEMI,dt,PLOT);
    Dirctn = -gradient;
    DirctnItrtn(:,ii) = Dirctn;
    
    %=== Plotting Search Direction ===%
    if PLOT.AdjStSearchDirection == 1;
        figure(PLOT.Figure_AdjStSearchDirection)
        trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),-gradient);
        if PLOT.DGMPlotBirdsEyeView == 1;
            view(2)
        end
        shading interp %thanks Ru!
        colormap(jet(256))
        title(PLOT.Figure_AdjStSearchDirection_Title,'FontWeight','bold')
    end

    %=== Computing Step Size ===%
    if RunOptions.EWE_SD_LS_Analytical == 1
        error('Acoustic forward problem is not linear with respect to wave speeds, cannot use analytical stepsize')
    end
    if RunOptions.EWE_SD_LS_fminsearch == 1
        disp('Using Line Search with fminsearch to compute next estimate');
        options = optimset('Display','iter','PlotFcns',@optimplotfval);
        if RunOptions.EWE_SD_ExponentialTransform == 1
            cpRecon = log(cpRecon);
            csRecon = log(csRecon);
        end
        [cSDStepSize,~] = fminsearch(@(alpha) SDLSFunctional_c(RunOptions,FrwdFunction,cpRecon + alpha(1)*cpDirctn,cpRecon + alpha(2)*cpDirctn,p0DGM,...
                                                               L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                               T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,...
                                                               T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                                               Prior,MeshIN_Elm,MeshINodes,PrecomputedIntrplteObjectsI,SensorsI,...
                                                               pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT),0,options);
                                                        
    end
    if RunOptions.EWE_SD_LS_PaulsFunction == 1;
        error('Pauls line search currently not compatible with reconstructing wave speeds')
    end
    printf(['Step Size for cp: ' num2str(cSDStepSize(1))])
    printf(['Step Size for cs: ' num2str(cSDStepSize(2))])
    
    %=== New cpRecon and csRecon Values ===%
    cpRecon = cpRecon + cSDStepSize(1)*cpDirctn;
    csRecon = csRecon + cSDStepSize(2)*csDirctn;
    if RunOptions.EWE_SD_ExponentialTransform == 1
        cpRecon = exp(cpRecon);
        csRecon = exp(csRecon);
    end
    toc
end

%=== If Steepest-Descent Iterations Were Terminated ===%
Trmntn = ii;
printf(['Terminated after ' num2str(Trmntn) ' iterations.']);

if (Trmntn==EWE_SD_MaxItrtn)
    printf('Maximum iterations reached.');
end

cpRecon = cpReconItrtn(:,Trmntn);
csRecon = csReconItrtn(:,Trmntn);

%% =======================================================================%
%                                Outputs
%=========================================================================%
AcousticInverseItrtnInfo.cpRecon = cpReconItrtn;
AcousticInverseItrtnInfo.csRecon = csReconItrtn;
AcousticInverseItrtnInfo.T11cRecon = T11cReconItrtn;
AcousticInverseItrtnInfo.T22cRecon = T22cReconItrtn;
AcousticInverseItrtnInfo.T12cRecon = T12cReconItrtn;
AcousticInverseItrtnInfo.vxcRecon = vxcReconItrtn;
AcousticInverseItrtnInfo.vycRecon = vycReconItrtn;
AcousticInverseItrtnInfo.fp0ErrorSize = fcErrorSizeItrtn;
AcousticInverseItrtnInfo.cpRelativeError = cpRelativeErrorItrtn;
AcousticInverseItrtnInfo.Dirctn = DirctnItrtn;

SaveFileNameReconstructions = sprintf('Reconstructions-%s',RunOptions.SaveFileName)
save(SaveFileNameReconstructions,'DataVrblsWave','Prior','cpTrue','csTrue','AcousticInverseItrtnInfo','PLOT','MeshINodes','Trmntn','-v7.3')
keyboard
end

%% =======================================================================%
%                        Line Search Functional
%=========================================================================%
function SD_LS_RelError = SDLSFunctional_c(RunOptions,FrwdFunction,NewcpRecon,NewcsRecon,p0DGM,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,Prior,MeshIN_Elm,MeshINodes,PrecomputedIntrplteObjectsI,SensorsI,pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT)

% SDLSFunctional computes ||Data - Newf(p0)|| using the new p0
% value where Newp0 = p0 + SDStepSize*Dirctn. Then SDLSFunctional can
% be considered to be a function of SDStepSize and thus used with
% fminsearch for linesearch purposes.
%
% Hwan Goh 18/9/2017, University of Auckland, New Zealand
disp(' ')
disp('------------------------------------------------------')
disp('Using fminsearch ')
disp('------------------------------------------------------')
printf(['For the Case ' RunOptions.SaveFileName]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computing and Interpolating f(cRecon + alpha*Dirctn) %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Computing f(cRecon + alpha*Dirctn) ===%
if RunOptions.EWE_SD_ExponentialTransform == 1
    p0DGM = exp(p0DGM);
end
Newp0ReconDGM = IntrplteOver2DTriangulatedMesh(MeshIN_Elm,MeshINodes,p0DGM,xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
IniCond = reshape(Newp0ReconDGM,NpI,KI)./(rhoI*RunOptions.SpecificHeatCoeff);
pinfoI(2).cpM = NewcpRecon;
pinfoI(2).csM = NewcsRecon;
PLOT.DGMForward = 0;
[T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps] = FrwdFunction(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT); 
[~,~,~,~,~,~,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps] = EWE_DGM2D_ObjectiveFunctionalErrorTerm(RunOptions,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,SensorsI);
L_EError = [L_ET11ErrorTimeSteps;L_ET22ErrorTimeSteps;L_ET12ErrorTimeSteps;L_EvxErrorTimeSteps;L_EvyErrorTimeSteps];

%%%%%%%%%%%%%%%%%%%%%%
%%% Relative Error %%% 
%%%%%%%%%%%%%%%%%%%%%%
SD_LS_Error = (1/2)*norm(L_EError,2) + (1/2)*norm(Prior.L_pr*(p0DGM - Prior.Exp_p0),2);
SD_LS_RelError = 100*(SD_LS_Error)/(norm([T11DataTimeSteps(:);T22DataTimeSteps(:);T12DataTimeSteps(:);vxDataTimeSteps(:);vyDataTimeSteps(:)],2));
end

