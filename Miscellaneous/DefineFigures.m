if PLOT.MeshD == 1
    PLOT.Figure_MeshD = figure;
    PLOT.Figure_MeshD_Title = 'Data Mesh';
    drawnow
end

if PLOT.MeshI == 1
    PLOT.Figure_MeshI = figure;
    PLOT.Figure_MeshI_Title = 'Inversion Mesh';
    drawnow
end

if PLOT.mu_aD == 1
    PLOT.Figure_mu_aD = figure;
    PLOT.Figure_mu_aD_Title = 'Absorption Coefficient';
    drawnow
end

if PLOT.mu_sD == 1
    PLOT.Figure_mu_sD = figure;
    PLOT.Figure_mu_sD_Title = 'Scattering Coefficient';
    drawnow
end

if PLOT.mu_sI == 1
    PLOT.Figure_mu_sI = figure;
    PLOT.Figure_mu_sI_Title = 'Scattering Coefficient';
    drawnow
end

if PLOT.phi==1;
    PLOT.Figure_phi = figure;
    PLOT.Figure_phi_Title = 'Light Fluence';
    drawnow
end

if PLOT.AbsorbedEnergyDensity==1;
    PLOT.AbsorbedEnergyDensity = figure;
    PLOT.Figure_AbsorbedEnergyDensity_Title = 'Absorbed Energy Density';
    movegui(PLOT.Figure_AbsorbedEnergyDensity,'north');
    drawnow
end

if PLOT.WaveDGMForwardMesh==1;
    PLOT.Figure_WaveDGMForwardMesh = figure;
    PLOT.Figure_WaveDGMForwardMesh_Title = 'Forward Mesh';
    movegui(PLOT.Figure_WaveDGMForwardMesh,'northwest');
    drawnow
end

if PLOT.WaveDGMInverseMesh==1;
    PLOT.Figure_WaveDGMInverseMesh = figure;
    PLOT.Figure_WaveDGMInverseMesh_Title = 'Inversion Mesh';
    movegui(PLOT.Figure_WaveDGMInverseMesh,'northwest');
    drawnow
end

if PLOT.IntrplteAbsorbedEnergyDensity==1;
    PLOT.Figure_AbsorbedEnergyDensityInterpolated = figure;
    PLOT.Figure_AbsorbedEnergyDensityInterpolated_Title = 'Absorbed Energy Density Interpolated';
    movegui(PLOT.Figure_AbsorbedEnergyDensityInterpolated,'north');
    drawnow
end

if PLOT.DGMForward==1
    PLOT.Figure_WavePropagation_Data = figure;
    PLOT.Figure_WavePropagation_Data_Title = 'Elastic Wave Propagation';
    movegui(PLOT.Figure_WavePropagation_Data,'southwest');
    drawnow
end

if PLOT.DGMForwardQuiver == 1;
    PLOT.Figure_WavePropagationQuiver_Data = figure;
    PLOT.Figure_WavePropagation_Data_Title = 'Elastic Wave Propagation';
    drawnow
end

if RunOptions.UseFullDomainData ~= 1 && PLOT.DGMForwardSensorData == 1;
    PLOT.Figure_WavePropagation_SensorData = figure;
    PLOT.Figure_WavePropagation_SensorData_Title = 'Sensor Data';
    drawnow
end

if PLOT.PriorMargPoints == 1;
    PLOT.Figure_PriorMargPoints = figure;
    PLOT.Figure_Prior_Title = 'Marginalisation Points';
    drawnow
end

if PLOT.PriorHist == 1;
    PLOT.Figure_PriorHist = figure;
    PLOT.Figure_PriorHist_Title = 'Prior Histogram';
    drawnow
end

if RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensity == 1
if PLOT.hRecon == 1;
    PLOT.Figure_CurrenthRecon = figure;
    PLOT.Figure_CurrenthRecon_Title = 'Reconstructed Absorbed Energy Density';
    movegui(PLOT.Figure_CurrenthRecon,'northeast');
    drawnow
end
end

if PLOT.AdjWaveField == 1;
    PLOT.Figure_AdjWaveField = figure;
    PLOT.Figure_AdjWaveField_Title = 'Adjoint Wave Field Propagation';
    movegui(PLOT.Figure_AdjWaveField,'south');
    drawnow
end

if PLOT.AdjStSearchDirection == 1;
    PLOT.Figure_AdjStSearchDirection = figure;
    PLOT.Figure_AdjStSearchDirection_Title = 'Search Direction';
    movegui(PLOT.Figure_AdjStSearchDirection,'southeast');
    drawnow
end

if PLOT.mu_aRecon == 1;
    PLOT.Figure_Currentmu_aRecon = figure;
    PLOT.Figure_Currentmu_aRecon_Title = 'Reconstructed Absorption Coefficient';
    movegui(PLOT.Figure_Currentmu_aRecon,'northeast');
    drawnow
end