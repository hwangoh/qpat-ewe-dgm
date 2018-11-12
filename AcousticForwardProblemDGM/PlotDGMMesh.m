%% =======================================================================%
%                      Plotting DGM Forward Mesh
%=========================================================================%
if PLOT.WaveDGMForwardMesh==1;
    figure(PLOT.Figure_WaveDGMForwardMesh)
    triplot(MeshD.Elements,MeshD.Nodes(:,1)*(1/RunOptions.ScalingOptical),MeshD.Nodes(:,2)*(1/RunOptions.ScalingOptical)); %In metres
    hold on
    title(PLOT.Figure_WaveDGMForwardMesh_Title,'FontWeight','bold')
    axis 'image'
    
    %=== All Nodes ===%
    plot(DGMMeshD.x,DGMMeshD.y,'go');
    
    %=== Sensor Locations ===%
    plot(DataVrblsWave.SensorCoords(:,1),DataVrblsWave.SensorCoords(:,2),'ro');
       
%     %=== Boundary Nodes ===%
%     plot(DGMMeshD.x(DGMMeshD.pinfo(2).vmapB),DGMMeshD.y(DGMMeshD.pinfo(2).vmapB),'ro');
%     
%     %=== Top Boundary Nodes ===%
%         plot(DGMMeshD.x(vmapTop(i)),DGMMeshD.y(vmapTop(i)),'ro');
%         
%     %=== Left Boundary Nodes ===%
%         plot(DGMMeshD.x(vmapLeft(i)),DGMMeshD.y(vmapLeft(i)),'ro');
%  
%     %=== Right Boundary Nodes ===%
%         plot(DGMMeshD.x(vmapRight(i)),DGMMeshD.y(vmapRight(i)),'ro');
% 
%     %=== Bttm Boundary Nodes ===%
%         plot(DGMMeshD.x(vmapBttm(i)),DGMMeshD.y(vmapBttm(i)),'ro');
%
%     %=== All DGM element boundary grid points ===%
%     for i=1:DGMMeshD.Np*DGMMeshD.K;
%         plot(DGMMeshD.x(i),DGMMeshD.y(i),'ro');
%         NodesInfo(i,1) = DGMMeshD.vmapM(i);
%         NodesInfo(i,2:3) = [DGMMeshD.x(DGMMeshD.vmapM(i)),DGMMeshD.y(DGMMeshD.vmapM(i))];
%         hold on
%         pause(0.5)
%     end
%     
%     %=== All DGM grid points ===%
%     %DGM grid points are called element wise. Use this plotting code to see the order of grid points.
%     for i=1:size(DGMMeshD.x(:));
%         i
%         [DGMMeshD.x(i),DGMMeshD.y(i)]
%         plot(DGMMeshD.x(i),DGMMeshD.y(i),'ro');
%         hold on
%         keyboard
%     end
%     
%     %=== Element Vertices ===%
%     %DGM grid points are called element wise. Use this plotting code to see the order of grid points.
%     for i=1:size(VertexNodesGlobalIndices(:));
%         plot(DGMMeshD.x(DGMMeshD.VertexNodesGlobalIndicesFEM(i)),DGMMeshD.y(DGMMeshD.VertexNodesGlobalIndicesFEM(i)),'ro');
%         hold on
%         pause(0.1)
%     end
%     keyboard
end
if RunOptions.InverseCrime == 0;
if PLOT.WaveDGMInverseMesh == 1;
    figure(PLOT.Figure_WaveDGMInverseMesh)
    triplot(MeshI.Elements,MeshI.Nodes(:,1)*(1/RunOptions.ScalingOptical),MeshI.Nodes(:,2)*(1/RunOptions.ScalingOptical)); %In metres
    hold on
    axis 'image'
    title(PLOT.Figure_WaveDGMInverseMesh_Title,'FontWeight','bold')
    %=== All Nodes ===%
    plot(DGMMeshI.x,DGMMeshI.y,'go');
    %=== Sensor Locations ===%
    plot(DataVrblsWave.SensorCoords(:,1),DataVrblsWave.SensorCoords(:,2),'ro');
%     %=== Sensor Normal Vectors ===%
%     if PLOT.WaveDGMSensorNormals==1
%         for ii=1:size(DataVrblsWave.SensorsI,2)
%             ii
%             hold on
%             quiver(DataVrblsWave.SensorsI{ii}.xy(1),DataVrblsWave.SensorsI{ii}.xy(2),DataVrblsWave.SensorsI{ii}.NormalVector(1),DataVrblsWave.SensorsI{ii}.NormalVector(2),0.0005,'color','r');
%         end
%     end
end
end