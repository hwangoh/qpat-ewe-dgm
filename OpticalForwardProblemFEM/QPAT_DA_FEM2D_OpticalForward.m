% QPAT_DA_FEM2D_OpticalForward solves the optical forward problem
%
% Hwan Goh, University of Auckland, New Zealand - 6/7/2015

%% %%%%%%%%%%%%%%%%%%
%%% Light Fluence %%%
%%%%%%%%%%%%%%%%%%%%%                    
%Simulated phi
[DataVrblsOptical.phi_temp,~,~,DataVrblsOptical.R,DataVrblsOptical.b,~] = DA_FEM2D_OpticalForward(RunOptions,MeshD,PrmtrsD);
DataVrblsOptical.phi = DataVrblsOptical.phi_temp; %this is done to avoid the "Illegal right hand side in assignment. Too many elements" error
%Plotting
if PLOT.phi==1;
    figure(PLOT.Figure_phi)
    trisurf(PLOT.TRI_MeshD,MeshD.Nodes(:,1),MeshD.Nodes(:,2),full(DataVrblsOptical.phi));
    view(2)
    shading interp %thanks Ru!
    colorbar
    colormap(jet(256))
    axis 'image'
    title(PLOT.Figure_phi_Title,'FontWeight','bold')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Absorbed Energy Density %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataVrblsOptical.h = (PrmtrsD.mu_a).*(DataVrblsOptical.phi); %Absorbed energy density

%=== Plotting Absorbed Energy Density ===%
if PLOT.AbsorbedEnergyDensity==1;
    figure(PLOT.Figure_AbsorbedEnergyDensity)
    trisurf(PLOT.TRI_MeshD,MeshD.Nodes(:,1),MeshD.Nodes(:,2),full(DataVrblsOptical.h));
    view(2)
    shading interp %thanks Ru!
    colorbar
    colormap(jet(256))
    axis 'image'
    title(PLOT.Figure_AbsorbedEnergyDensity_Title,'FontWeight','bold')
end