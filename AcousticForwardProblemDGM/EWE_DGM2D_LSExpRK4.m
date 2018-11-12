function [T11TimeSteps,T22TimeSteps,T12TimeSteps,vxTimeSteps,vyTimeSteps] = EWE_DGM2D_LSExpRK4(RunOptions,IniCond,x,y,Np,K,pinfo,dt,PLOT)

% EWE_DGM2D_LSExpRK4 computes the DGM acoustic forward problem to generate 
% sensory data using the low-storage five-stage explicit fourth order Runge-Kutta method.
%
% Inputs:
%   RunOptions:
%      NumberofTimeSteps - Number of time steps to be computed
%      BC - Boundary conditions
%   IniCond - Initial condition for QPAT wave propagation
%   x: x-coordinate of the nodes of the DGM mesh [m]
%   y: y-coordinate of the nodes of the DGM mesh [m]
%   Np - Number of nodes per element in DGM grid
%   K - Number of elements
%   dt - Size of time step
%   PLOT: To plot or not to plot, that is the question
%
% Outputs:
%   T11TimeSteps - Time step values for T11
%   T22TimeSteps - Time step values for T22
%   T12TimeSteps - Time step values for T12
%   vxTimeSteps -  Time step values for vx
%   vyTimeSteps -  Time step values for vy
%
% Hwan Goh, University of Auckland, New Zealand 11/6/2016 (Hyv‰‰ huomenta from Kuopio!)
% Last Edited: 24/11/2017 - Separated out sensory data definition from output

disp('Computing Forward Wave Field')
tic

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
%                          Initial Structures
%=========================================================================%
%=== Time Step Storage ===%
T11TimeSteps = zeros(Np*K,RunOptions.NumberofTimeSteps);
T22TimeSteps = zeros(Np*K,RunOptions.NumberofTimeSteps);
T12TimeSteps = zeros(Np*K,RunOptions.NumberofTimeSteps);
vxTimeSteps = zeros(Np*K,RunOptions.NumberofTimeSteps);
vyTimeSteps = zeros(Np*K,RunOptions.NumberofTimeSteps);

%=== Initial Time Step Values ===%
%QPAT Initial Condition
T11TimeSteps(:,1) = IniCond(:);
T22TimeSteps(:,1) = IniCond(:);
T12TimeSteps(:,1) = zeros(Np*K,1);
vxTimeSteps(:,1) = zeros(Np*K,1);
vyTimeSteps(:,1) = zeros(Np*K,1);
    
%=== Runge-Kutta residual storage ===%
res_T11 = zeros(Np*K,1);
res_T22 = zeros(Np*K,1);
res_T12 = zeros(Np*K,1);
res_vx = zeros(Np*K,1);
res_vy = zeros(Np*K,1);

%% =======================================================================%
%                             Time Stepping
%=========================================================================%
%=== Set Up for Saving Movie ===%
% mov = VideoWriter('WaveElasticLayerDomain.avi');
% mov.Quality = 100; % From [0, 100]
% mov.FrameRate = 20;
% set(gca, 'nextplot','replacechildren', 'Visible','on');
% open(mov)

%=== Outer time step loop ===%
for tstep=2:RunOptions.NumberofTimeSteps
    %printf(['\nComputing time step number ' num2str(tstep) ' of ' num2str(RunOptions.NumberofTimeSteps)]);
    %tic
    T11 = T11TimeSteps(:,tstep-1);
    T22 = T22TimeSteps(:,tstep-1);
    T12 = T12TimeSteps(:,tstep-1);
    vx = vxTimeSteps(:,tstep-1);
    vy = vyTimeSteps(:,tstep-1);
    for INTRK = 1:5
        %timelocal = time + rk4c(INTRK)*dt; %DGMWave2DLSHS does not depend on time
        %=== Compute RHS ===%
        if RunOptions.StrainVelocityForm == 1;
            [rhs_T11,rhs_T22,rhs_T12,rhs_vx,rhs_vy] = EWE_DGM2D_PNonConf_Strain(RunOptions,pinfo,Np,K,T11,T22,T12,vx,vy);
        end
        if RunOptions.StressVelocityForm == 1;
            [rhs_T11,rhs_T22,rhs_T12,rhs_vx,rhs_vy] = EWE_DGM2D_PNonConf_Stress(RunOptions,pinfo,Np,K,T11,T22,T12,vx,vy);
        end
        
        %=== Initiate and increment Runge-Kutta residuals ===%
        res_T11 = rk4a(INTRK)*res_T11 + dt*rhs_T11;
        res_T22 = rk4a(INTRK)*res_T22 + dt*rhs_T22;
        res_T12 = rk4a(INTRK)*res_T12 + dt*rhs_T12;
        res_vx  = rk4a(INTRK)*res_vx  + dt*rhs_vx;
        res_vy  = rk4a(INTRK)*res_vy  + dt*rhs_vy;
        
        %=== Update Fields ===%
        T11 = T11 + rk4b(INTRK)*res_T11;
        T12 = T12 + rk4b(INTRK)*res_T12;
        T22 = T22 + rk4b(INTRK)*res_T22;
        vx  = vx  + rk4b(INTRK)*res_vx;
        vy  = vy  + rk4b(INTRK)*res_vy;
    end   
    
    %=== New values for next iteration ===%
    pDGMVeloTemp = sqrt(vx.^2 + vy.^2); %current pDGM derivative in time
    T11TimeSteps(:,tstep) = T11;
    T22TimeSteps(:,tstep) = T22;
    T12TimeSteps(:,tstep) = T12;
    vxTimeSteps(:,tstep) = vx;
    vyTimeSteps(:,tstep) = vy;
    
    %=== Plotting ===%
    if PLOT.DGMForward == 1;
        figure(PLOT.Figure_WavePropagation_Data)
        trisurf(PLOT.TRI_DGMMeshD,x,y,pDGMVeloTemp);
        shading interp
        colorbar
        colormap(jet(256))
        if PLOT.DGMPlotBirdsEyeView == 1;
            view(2)
        end
        if PLOT.DGMPlotUsezLim == 1;
            zlim([PLOT.DGMPlotzAxis]);
        end
        if PLOT.HoldColourAxis == 1;
            caxis(PLOT.ColourAxis);
        end
        title(PLOT.Figure_WavePropagation_Data_Title,'FontWeight','bold')
%         writeVideo(mov,getframe(gcf));
    end
    if PLOT.DGMForwardQuiver == 1;
        figure(PLOT.Figure_WavePropagationQuiver_Data)
        quiver(x(:),y(:),vx,vy);
        title(PLOT.Figure_WavePropagationQuiver_Data_Title,'FontWeight','bold')
    end
end
toc

