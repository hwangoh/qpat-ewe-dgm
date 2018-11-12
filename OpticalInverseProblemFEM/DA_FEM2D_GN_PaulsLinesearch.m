function [GNStepSize, Dirctn]=DA_FEM2D_GN_PaulsLinesearch(RunOptions,DataVrblsOptical,FrwdFunction,RegData,CurrError,Dirctn,MeshI,PrmtrsI,PrmtrsPrp,R,b_pattern,L_1,L_2,CurrGNItertn)

%DA_FEM2D_GN_PaulsLinesearch finds the optimal step length to take from xHat in dir assuming
% that we have an approximately quadratic objective functional.
%
%Inputs:
%   RunOptions: a structure with the various user defineable options; Defaults: 8,1e-4,2,1,1e-6,3
%              DA_GN_LS_MaxItrtn - max number of iterations for linesearch
%              tol - tolerance for linesearch, terminate if
%                           abs(GNalpha_i-GNalpha_{i-1})/GNalpha_{i-1}<tol
%              DA_GN_LS_MaxStepSize - max step size 
%              DA_GN_LS_PostivityConstraint - Whether positivity constraints are implemented
%              DA_GN_LS_mu_aMin - minimum value for mu_a
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
% Hwan Goh 28/12/2015, University of Auckland, New Zealand
% adapted from P. J. Hadwin 11/06/2011, University of Auckland, New Zealand
% TO CONSIDER - having partially positive i.e. elements 4 and 7 have to be
% positive but rest can be neg. Idea - use indicator for positive elements
%             - Fix up debugging!!!! i.e. when GNalpha is Nan or negative

%=========================================================================%
%                             Initial Setup
%=========================================================================%
%=== Setting User Defined Options ===%
DA_GN_LS_MaxItrtn = RunOptions.DA_GN_LS_MaxItrtn;
DA_GN_LS_MaxStepSize = RunOptions.DA_GN_LS_MaxStepSize; 
DA_GN_LS_PostivityConstraint = RunOptions.DA_GN_LS_PostivityConstraint;
DA_GN_LS_mu_aMin = RunOptions.DA_GN_LS_mu_aMin;

%=== Initial Structures ===%
StepSizes=[0; 1; zeros(DA_GN_LS_MaxItrtn,1)];
ErrorVals=[norm(CurrError);zeros(DA_GN_LS_MaxItrtn,1)];
ErrorValsmin=ErrorVals(1); %First ymin is just error
minItrtnNum=1; %Iteration number of current ymin
HasMinChanged=[0,0]; %If a new minimum has been found, then MinChanged(2) = MinChanged(2) + 1
GNStepSize=StepSizes(minItrtnNum); %Corresponding stepsize that generates current ymin
i=2;

%=========================================================================%
%                         Line Search Iterations
%=========================================================================%
for i=2:DA_GN_LS_MaxItrtn
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculating the New Error with New Step Size %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PrmtrsI.mu_a=PrmtrsI.mu_a+StepSizes(i)*Dirctn;
    
    if DA_GN_LS_PostivityConstraint
        D=find(PrmtrsI.mu_a<0);
        if ~isempty(D)
            printf(['LSPositivity Constraint applied during ' num2str(i) 'th LS iteration, ' num2str(CurrGNItertn) 'th GN iteration']);
        end
        PrmtrsI.mu_a(PrmtrsI.mu_a<0)=DA_GN_LS_mu_aMin;
    end
    
    [fmu_a,~,~,~,~,~] = FrwdFunction(RunOptions,DataVrblsOptical,MeshI,PrmtrsI.mu_a,PrmtrsI.mu_s,PrmtrsPrp,R,b_pattern);
    Regfmu_a = [L_1*fmu_a;L_2*PrmtrsI.mu_a];
        
    Error = RegData - Regfmu_a;
    ErrorVals(i) = norm(Error);
    
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
        GNStepSize=StepSizes(minItrtnNum);
    end
    
    %=== Find Closest Step Size to Current Step Size ===%
    g=abs(GNStepSize-StepSizes(1:i));
    II=find(g==min(g(g~=0)), 1, 'last'); %g~=0 is a boolean vector checking each entry of g to see if the entry equals 0. Therefore g(g~=0) outputs all nonzero values of g. 
                                         %Then min(g(g~=0)) finds the smallest nonzero value of g. Thus g==min(g(g~=0)) is a boolean vector checking each entry of g to
                                         %see the entry equals its smallest value. Finally, find(g==min(g(g~=0)), 1, 'last') will output the entry number of where the entry of g is the
                                         %smallest value of g. If there are repeats then it will output the furthest entry number. 
    %=== Next Point ===%
    if (GNStepSize<StepSizes(II))  %If closest step size is larger than current step size
        if (HasMinChanged==0) %If minimum has not changed
            StepSizes(i+1)=GNStepSize+0.5*abs(diff([StepSizes(II) GNStepSize]));
        else %If minimum has changed
            StepSizes(i+1)=GNStepSize+(-1)^HasMinChanged(2)*0.5*abs(diff([StepSizes(II) GNStepSize]));
        end
    else %If closest step size is smaller than current step size
        if (HasMinChanged==0) %If minimum has not changed
            StepSizes(i+1)=GNStepSize+abs(diff([StepSizes(II) GNStepSize]));
        else %If minimum has not changed
            StepSizes(i+1)=GNStepSize+0.5*abs(diff([StepSizes(II) GNStepSize]));
        end
    end
    if (StepSizes(i+1)==0)
        StepSizes(i+1)=GNStepSize/2;
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
if GNStepSize==0
    GNStepSize=0.5;
end
GNStepSize=min(GNStepSize,DA_GN_LS_MaxStepSize);

