function [StepSize]=EWE_DGM2D_LS_PaulsLinesearch(RunOptions,FrwdFunction,p0Recon,Dirctn,CurrentObjFunctError,CurrLSItertn,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,Prior,MeshIN_Elm,MeshINodes,PrecomputedIntrplteObjectsI,SensorsI,pinfoI,xI,yI,NpI,KI,rhoI,dt,PLOT)

% EWE_DGM2D_LS_PaulsLinesearch finds the optimal step length to take from xHat in dir assuming
% that we have an approximately quadratic objective functional.
%
%Inputs:
%   RunOptions: a structure with the various user defineable options; Defaults: 8,1e-4,2,1,1e-6,3
%              EWE_LS_MaxItrtn - max number of iterations for linesearch
%              tol - tolerance for linesearch, terminate if
%                           abs(GNalpha_i-GNalpha_{i-1})/GNalpha_{i-1}<tol
%              EWE_LS_MaxStepSize - max step size 
%              EWE_LS_PostivityConstraint - Whether positivity constraints are implemented
%              EWE_LS_mu_aMin - minimum value for mu_a
%              pnt ~ number of data points calculated for linesearch
%                           decreasing speeds up but decreases accuracy
%                           increasing slows down but increases accuracy
%   DataVrblsOptical - T: Grunaisen parameter
%   RegData - target regulated data vector
%   Dirctn - search direction
%   CurrError - error = ||RegData - Regfmu_a|| of current Gauss-Newton iteration
%   Mesh - computational mesh
%   PrmtrsI - inputs from GNIterations
%   PrmtrsPrp:
%      g - mean of the cosine of the scattering angle, used for computing the diffusion coefficient
%   R - boundary matrix
%   b_pattern - N_Patterns by N_Nodes Matrix where each row is a lightsource vector with deactivated light sources
%   L_1 - Cholesky factor of a weighting operator
%   L_2 - Regularization operator
%   CurrGNItertn - CurrentGauss-Newton iteration number
%   Varagin:
%      For DGM:
%        FrwdMatrixI - Matrix that represents the linear forward acoustic solver on the inverse mesh
%        PrecompIntrplteObjects - Objects that depend on the Mesh nodes. May have been computed earlier and so can 
%                            be called here to avoid repeating computations. Set to 0 in function call if you
%                            want to compute new objects for a new mesh.
%
%Outputs:
%   GNStepSize - approximate optimal step size from xHat in dir
%   dir - direction vector, may have been modified if constraints present
%
% Hwan Goh 17/02/2018, University of Auckland, New Zealand
% adapted from P. J. Hadwin 11/06/2011, University of Auckland, New Zealand
% TO CONSIDER - having partially positive i.e. elements 4 and 7 have to be
% positive but rest can be neg. Idea - use indicator for positive elements
%             - Fix up debugging!!!! i.e. when GNalpha is Nan or negative

%=========================================================================%
%                             Initial Setup
%=========================================================================%
%=== Setting User Defined Options ===%
EWE_LS_MaxItrtn = RunOptions.EWE_SD_LS_MaxItrtn;
EWE_LS_MaxStepSize = RunOptions.EWE_SD_LS_MaxStepSize; 
EWE_LS_PostivityConstraint = RunOptions.EWE_SD_LS_PostivityConstraint;
EWE_LS_p0Min = RunOptions.EWE_SD_LS_p0Min;
PLOT.DGMForward = 0; %Suppress plotting of wave propagation

%=== Initial Structures ===%
StepSizes=[0; RunOptions.EWE_LS_StepSize; zeros(EWE_LS_MaxItrtn,1)];
ErrorVals=[norm(CurrentObjFunctError);zeros(EWE_LS_MaxItrtn,1)];
ErrorValsmin=ErrorVals(1); %First ymin is just error
minItrtnNum=1; %Iteration number of current ymin
HasMinChanged=[0,0]; %If a new minimum has been found, then MinChanged(2) = MinChanged(2) + 1
StepSize=StepSizes(minItrtnNum); %Corresponding stepsize that generates current ymin
i=2;

%=========================================================================%
%                         Line Search Iterations
%=========================================================================%
for i=2:EWE_LS_MaxItrtn
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculating the New Error with New Step Size %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p0Recon=p0Recon+StepSizes(i)*Dirctn;
    
    %=== Positivity Constraint ===%
%     if EWE_SD_LS_PostivityConstraint
%         D=find(p0Recon<0);
%         if ~isempty(D)
%             printf(['LSPositivity Constraint applied during ' num2str(i) 'th LS iteration, ' num2str(CurrSDItertn) 'th GN iteration']);
%         end
%         p0Recon(p0Recon<0)=EWE_SD_LS_p0Min;
%     end
    
    %== Applying Transformations ===%
    if RunOptions.EWE_LS_ExponentialTransform == 1
        p0Recon = exp(p0Recon);
    end
    if RunOptions.EWE_LS_SigmoidTransform == 1
        p0Recon = Prior.SigTalpha + Prior.SigTbeta*tanh(Prior.SigTkappa*(p0Recon - Prior.SigTalpha));
    end
    
    %=== Forward Wave Propagation ===%
    p0ReconDGM = IntrplteOver2DTriangulatedMesh(MeshIN_Elm,MeshINodes,p0Recon,xI,yI,NpI*KI,PrecomputedIntrplteObjectsI);
    IniCond = reshape(p0ReconDGM,NpI,KI)./(rhoI*RunOptions.SpecificHeatCoeff);
    [T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps] = FrwdFunction(RunOptions,IniCond,xI,yI,NpI,KI,pinfoI,dt,PLOT);
    [~,~,~,~,~,~,L_ET11ErrorTimeSteps,L_ET22ErrorTimeSteps,L_ET12ErrorTimeSteps,L_EvxErrorTimeSteps,L_EvyErrorTimeSteps] = EWE_DGM2D_ObjectiveFunctionalErrorTerm(RunOptions,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy,...
                                                                                                                                                        T11DataTimeSteps,T22DataTimeSteps,T12DataTimeSteps,vxDataTimeSteps,vyDataTimeSteps,...
                                                                                                                                                        T11ReconTimeSteps,T22ReconTimeSteps,T12ReconTimeSteps,vxReconTimeSteps,vyReconTimeSteps,...
                                                                                                                                                        T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                                                                                                                                        SensorsI);                                                                                                                                                        
    L_EError = [L_ET11ErrorTimeSteps;L_ET22ErrorTimeSteps;L_ET12ErrorTimeSteps;L_EvxErrorTimeSteps;L_EvyErrorTimeSteps];
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Relative Error %%% 
    %%%%%%%%%%%%%%%%%%%%%%
    ErrorVals(i) = (1/2)*norm(L_EError,2) + (1/2)*norm(Prior.L_pr*(p0Recon - Prior.Exp_p0),2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Checking if New Mininum Value has been Found %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [NewMin, I]=min(ErrorVals(ErrorVals>0));
    if I == minItrtnNum % if min hasn't changed
        HasMinChanged(2)=HasMinChanged(2)+1;
    else %if min changed
        ErrorValsmin=NewMin; 
        minItrtnNum=I; 
        HasMinChanged(2)=0;
        StepSize=StepSizes(minItrtnNum);
    end
    
    %=== Find Closest Step Size to Current Step Size ===%
    g=abs(StepSize-StepSizes(1:i));
    II=find(g==min(g(g~=0)), 1, 'last'); %g~=0 is a boolean vector checking each entry of g to see if the entry equals 0. Therefore g(g~=0) outputs all nonzero values of g. 
                                         %Then min(g(g~=0)) finds the smallest nonzero value of g. Thus g==min(g(g~=0)) is a boolean vector checking each entry of g to
                                         %see the entry equals its smallest value. Finally, find(g==min(g(g~=0)), 1, 'last') will output the entry number of where the entry of g is the
                                         %smallest value of g. If there are repeats then it will output the furthest entry number. 
    %=== Next Point ===%
    if (StepSize<StepSizes(II))  %If closest step size is larger than current step size
        if (HasMinChanged==0) %If minimum has not changed
            StepSizes(i+1)=StepSize+0.5*abs(diff([StepSizes(II) StepSize]));
        else %If minimum has changed
            StepSizes(i+1)=StepSize+(-1)^HasMinChanged(2)*0.5*abs(diff([StepSizes(II) StepSize]));
        end
    else %If closest step size is smaller than current step size
        if (HasMinChanged==0) %If minimum has not changed
            StepSizes(i+1)=StepSize+abs(diff([StepSizes(II) StepSize]));
        else %If minimum has not changed
            StepSizes(i+1)=StepSize+0.5*abs(diff([StepSizes(II) StepSize]));
        end
    end
    if (StepSizes(i+1)==0)
        StepSizes(i+1)=StepSize/2;
    end
    StepSizes=abs(StepSizes); jj=0;
    while ~isempty(find(StepSizes(1:i)==StepSizes(i+1), 1))
        StepSizes(i+1)=(StepSizes(i+1))*(1+0.2*(-1)^jj);
        jj=jj+1;
        if jj>10
            printf('Linesearch NOTE: Not Unique point, something wrong?')
            break
        end
    end
    StepSizes=abs(StepSizes);
    
%     %plot points
%     figure(500),clf,plot(ALPHA(1:ii),Y(1:ii),'o'),hold on 
%     plot(GNalpha,ymin,'or'),hold off,
%     axis([0,max([ALPHA;2]),min(Y(Y>0)),max(Y)]),drawnow
end
if StepSize==0
    StepSize=0.5;
end
StepSize=min(StepSize,EWE_LS_MaxStepSize);

