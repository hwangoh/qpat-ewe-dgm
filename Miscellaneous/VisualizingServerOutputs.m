close all
clc

PLOT.DGMPlotBirdsEyeView = 1;

PLOT.Figure_AbsorbedEnergyDensity = figure;
PLOT.Figure_AbsorbedEnergyDensity_Title = 'Absorbed Energy Density';
movegui(PLOT.Figure_AbsorbedEnergyDensity,'northwest');
drawnow

PLOT.Figure_CurrenthRecon = figure;
PLOT.Figure_CurrenthRecon_Title = 'Reconstructed Absorbed Energy Density';
movegui(PLOT.Figure_CurrenthRecon,'northeast');
drawnow

PLOT.Figure_AdjStSearchDirection = figure;
PLOT.Figure_AdjStSearchDirection_Title = 'Search Direction';
movegui(PLOT.Figure_AdjStSearchDirection,'south');
drawnow

%=== Plotting Absorbed Energy Density ===%
figure(PLOT.Figure_AbsorbedEnergyDensity)
trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),hTrueIntrplte);
if PLOT.DGMPlotBirdsEyeView == 1;
    view(2)
end
shading interp %thanks Ru!
if PLOT.DGMPlotUsezLim == 1;
    zlim(PLOT.AbsorbedEnergyDensityzAxis);
end
colorbar
caxis([0 200])
colormap(jet(256))
title(PLOT.Figure_AbsorbedEnergyDensity_Title,'FontWeight','bold')

for ii = 1:100
    %=== Plotting Search Direction ===%
    figure(PLOT.Figure_AdjStSearchDirection)
    trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),full(AcousticInverseItrtnInfo.Dirctn(:,ii)));
    if PLOT.DGMPlotBirdsEyeView == 1;
        view(2)
    end
    shading interp %thanks Ru!
    colormap(jet(256))
    title(PLOT.Figure_AdjStSearchDirection_Title,'FontWeight','bold')
    
    %=== Plotting Reconstruction ===%
    figure(PLOT.Figure_CurrenthRecon)
    trisurf(PLOT.TRIFEM,MeshINodes(:,1),MeshINodes(:,2),full(AcousticInverseItrtnInfo.hRecon(:,ii)));
    if PLOT.DGMPlotBirdsEyeView == 1;
        view(2)
    end
    if PLOT.DGMPlotUsezLim == 1;
        zlim(PLOT.AbsorbedEnergyDensityzAxis);
    end
    shading interp %thanks Ru!
%     colorbar
    caxis([0 200])
    colormap(jet(256))
%     title(PLOT.Figure_CurrenthRecon_Title,'FontWeight','bold')
    ii
    keyboard
end