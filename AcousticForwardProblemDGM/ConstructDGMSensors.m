%% =======================================================================%
%                Place Sensors and Find Corresponding Nodes
%=========================================================================%
%=== Sensor Setups ===%
vec = linspace(-0.01, 0.01, RunOptions.NumberofSensorsOnOneBoundaryEdge)';
DataVrblsWave.SensorCoords = [vec  0.01*ones(RunOptions.NumberofSensorsOnOneBoundaryEdge, 1);
                              vec -0.01*ones(RunOptions.NumberofSensorsOnOneBoundaryEdge, 1);
                              0.01*ones(RunOptions.NumberofSensorsOnOneBoundaryEdge-2, 1) vec(2:end-1);
                             -0.01*ones(RunOptions.NumberofSensorsOnOneBoundaryEdge-2, 1) vec(2:end-1)];
% DataVrblsWave.SensorCoords = [DGMMeshD.x(DGMMeshD.pinfo(2).vmapB),DGMMeshD.y(DGMMeshD.pinfo(2).vmapB)];
if RunOptions.UseFullDomainData == 1 %Using full domain data
    DataVrblsWave.SensorCoords = [DGMMeshD.x',DGMMeshD.y'];
end
DataVrblsWave.NumberofSensors = size(DataVrblsWave.SensorCoords,1);

%=== Generate Sensors ===%
DataVrblsWave.SensorsD = EWE_DGM2D_ConstructSensorArray(RunOptions,DataVrblsWave.SensorCoords,DGMMeshD.VX,DGMMeshD.VY,DGMMeshD.EToV,MeshD.Bnd_ElmInd,MeshD.Bnd_NodeInd,DGMMeshD.pinfo,DGMMeshD.Norder);
DataVrblsWave.SensorsI = EWE_DGM2D_ConstructSensorArray(RunOptions,DataVrblsWave.SensorCoords,DGMMeshI.VX,DGMMeshI.VY,DGMMeshI.EToV,MeshI.Bnd_ElmInd,MeshI.Bnd_NodeInd,DGMMeshI.pinfo,DGMMeshI.Norder);
