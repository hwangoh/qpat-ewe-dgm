function [T11AdjTime0,T22AdjTime0,T12AdjTime0,vxAdjTime0,vyAdjTime0,T11AdjTimeStep,T22AdjTimeStep,T12AdjTimeStep,vxAdjTimeStep,vyAdjTimeStep] = EWE_DGM2D_LSExpRK4Reverse(RunOptions,T11DerivObjFunctTimeSteps,T22DerivObjFunctTimeSteps,T12DerivObjFunctTimeSteps,vxDerivObjFunctTimeSteps,vyDerivObjFunctTimeSteps,pinfoI,xI,yI,NpI,KI,dt,PLOT)

% EWE_DGM2D_LSExpRK4Reverse computes the DGM acoustic forward problem in reverse time stepping
% using the low-storage five-stage explicit fourth order Runge-Kutta method.
%
% Inputs:
%   RunOptions:
%      NumberofTimeSteps - Number of time steps to be computed
%      BC - Boundary conditions
%   T11DerivObjFunctTimeSteps - T11 term of derivative of objective functional
%   T22DerivObjFunctTimeSteps - T22 term of derivative of objective functional
%   T12DerivObjFunctTimeSteps - T12 term of derivative of objective functional
%   vxDerivObjFunctTimeSteps -  vx term of derivative of objective functional
%   vyDerivObjFunctTimeSteps -  vy term of derivative of objective functional
%   pinfoI - Information regarding p-refinement for p-nonconforming meshes on DGM inversion mesh
%   xI: x coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   yI: y coordinates of the DGM inversion mesh nodes, for interpolation onto optical inversion mesh and plotting
%   NpI - Number of nodes per element in inverse DGM grid
%   KI - Number of elements in inverse DGM grid
%   dt - Size of time step, need negative of this to do reverse time stepping
%   PLOT: To plot or not to plot, that is the question
%
% Outputs:
%   T11AdjTime0 - Adjoint variable of T11 at initial time
%   T22AdjTime0 - Adjoint variable of T22 at initial time
%   T12AdjTime0 - Adjoint variable of T12 at initial time
%   vxAdjTime0 - Adjoint variable of vx at initial time
%   vxAdjTime0 - Adjoint variable of vx at initial time
%   T11AdjTimeStep - Adjoint variable of T11
%   T22AdjTimeStep - Adjoint variable of T22
%   T12AdjTimeStep - Adjoint variable of T12
%   vxAdjTimeStep - Adjoint variable of vx
%   vyAdjTimeStep - Adjoint variable of vy
%
% Hwan Goh, University of Auckland, New Zealand 12/09/2017
% Last Edited: 22/11/2017 - set up for when sensor is not on node of inverse mesh

tic
disp('Computing Adjoint Wave Field')

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0]; 

%% =======================================================================%
%                           Initial Structures
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Structures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Time Step Storage ===%
T11AdjTimeStep = zeros(NpI*KI,RunOptions.NumberofTimeSteps);
T22AdjTimeStep = zeros(NpI*KI,RunOptions.NumberofTimeSteps);
T12AdjTimeStep = zeros(NpI*KI,RunOptions.NumberofTimeSteps);
vxAdjTimeStep = zeros(NpI*KI,RunOptions.NumberofTimeSteps);
vyAdjTimeStep = zeros(NpI*KI,RunOptions.NumberofTimeSteps);

%=== Runge-Kutta residual storage ===% 
res_T11 = zeros(NpI*KI,1); 
res_T22 = zeros(NpI*KI,1); 
res_T12 = zeros(NpI*KI,1); 
res_vx = zeros(NpI*KI,1);
res_vy = zeros(NpI*KI,1);

%% =======================================================================%
%                             Time Stepping
%=========================================================================%
%=== Outer time step loop ===%
for tstep=1:RunOptions.NumberofTimeSteps-1
    %=== Adjoint Variables at Current Time ===%
    T11AdjCurrent = T11AdjTimeStep(:,RunOptions.NumberofTimeSteps - (tstep-1));
    T22AdjCurrent = T22AdjTimeStep(:,RunOptions.NumberofTimeSteps - (tstep-1));
    T12AdjCurrent = T12AdjTimeStep(:,RunOptions.NumberofTimeSteps - (tstep-1));
    vxAdjCurrent = vxAdjTimeStep(:,RunOptions.NumberofTimeSteps - (tstep-1));
    vyAdjCurrent = vyAdjTimeStep(:,RunOptions.NumberofTimeSteps - (tstep-1));
    
    %=== Derivative of Objective Functional ===%
    T11DerivObjFunctCurrent = reshape(T11DerivObjFunctTimeSteps(:,RunOptions.NumberofTimeSteps - (tstep-1)),NpI,KI);
    T22DerivObjFunctCurrent = reshape(T22DerivObjFunctTimeSteps(:,RunOptions.NumberofTimeSteps - (tstep-1)),NpI,KI);
    T12DerivObjFunctCurrent = reshape(T12DerivObjFunctTimeSteps(:,RunOptions.NumberofTimeSteps - (tstep-1)),NpI,KI);
    vxDerivObjFunctCurrent = reshape(vxDerivObjFunctTimeSteps(:,RunOptions.NumberofTimeSteps - (tstep-1)),NpI,KI);
    vyDerivObjFunctCurrent = reshape(vyDerivObjFunctTimeSteps(:,RunOptions.NumberofTimeSteps - (tstep-1)),NpI,KI);
    
    %=== Compute RHS ===%  
    for INTRK = 1:5            
        %Spatially discretized adjoint wave field
        [rhs_T11,rhs_T22,rhs_T12,rhs_vx,rhs_vy] = EWE_DGM2D_PNonConf_StrainAdjoint(RunOptions,pinfoI,...
                                                                                   T11AdjCurrent,T22AdjCurrent,T12AdjCurrent,vxAdjCurrent,vyAdjCurrent,...
                                                                                   T11DerivObjFunctCurrent,T22DerivObjFunctCurrent,T12DerivObjFunctCurrent,vxDerivObjFunctCurrent,vyDerivObjFunctCurrent,...
                                                                                   NpI,KI);                                                                        
        %=== Initiate and increment Runge-Kutta residuals ===%
        res_T11 = rk4a(INTRK)*res_T11 + -dt*rhs_T11;
        res_T22 = rk4a(INTRK)*res_T22 + -dt*rhs_T22;
        res_T12 = rk4a(INTRK)*res_T12 + -dt*rhs_T12;
        res_vx  = rk4a(INTRK)*res_vx  + -dt*rhs_vx;
        res_vy  = rk4a(INTRK)*res_vy  + -dt*rhs_vy;
       
        %=== Update Fields ===%
        T11AdjCurrent = T11AdjCurrent + rk4b(INTRK)*res_T11;
        T12AdjCurrent = T12AdjCurrent + rk4b(INTRK)*res_T12;
        T22AdjCurrent = T22AdjCurrent + rk4b(INTRK)*res_T22;
        vxAdjCurrent  = vxAdjCurrent  + rk4b(INTRK)*res_vx;
        vyAdjCurrent  = vyAdjCurrent  + rk4b(INTRK)*res_vy;
    end
    
    %=== New values for next iteration ===%
    pDGMVeloTemp = sqrt(vxAdjCurrent.^2 + vyAdjCurrent.^2); %current pDGM derivative in time 
    T11AdjTimeStep(:,RunOptions.NumberofTimeSteps - tstep) = T11AdjCurrent;
    T22AdjTimeStep(:,RunOptions.NumberofTimeSteps - tstep) = T22AdjCurrent;
    T12AdjTimeStep(:,RunOptions.NumberofTimeSteps - tstep) = T12AdjCurrent;
    vxAdjTimeStep(:,RunOptions.NumberofTimeSteps - tstep) = vxAdjCurrent;
    vyAdjTimeStep(:,RunOptions.NumberofTimeSteps - tstep) = vyAdjCurrent;
    
  %=== Plotting ===%
  if PLOT.AdjWaveField == 1;
      figure(PLOT.Figure_AdjWaveField)
      trisurf(PLOT.TRI_DGMMeshI,xI,yI,pDGMVeloTemp);
      shading interp
      colormap(jet(256))
      if PLOT.DGMPlotBirdsEyeView == 1;
          view(2)
      end
      if PLOT.DGMPlotUsezLim == 1;
          zlim([PLOT.DGMAdjointPlotzAxis]);
      end
%       if PLOT.HoldColourAxis == 1;
%           caxis(PLOT.ColourAxis);
%       end
      title(PLOT.Figure_AdjWaveField_Title,'FontWeight','bold')
      pause(.001)
  end
  %=== Increment time ===%
  %time = time+dt; %DGMWave2DElstc does not depend on time  
  %toc
end
toc
%% =======================================================================%
%                              Output Data
%=========================================================================%
T11AdjTime0 = reshape(T11AdjTimeStep(:,1),NpI,KI);
T22AdjTime0 = reshape(T22AdjTimeStep(:,1),NpI,KI);
T12AdjTime0 = reshape(T12AdjTimeStep(:,1),NpI,KI);
vxAdjTime0 = reshape(vxAdjTimeStep(:,1),NpI,KI);
vyAdjTime0 = reshape(vyAdjTimeStep(:,1),NpI,KI);
