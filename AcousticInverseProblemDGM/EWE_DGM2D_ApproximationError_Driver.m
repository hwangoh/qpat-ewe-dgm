printf('Constructing noise regularization matrix with approximation error');
%=== File Names for Saving and Loading ===%
RunOptions.SaveFileNameAESamples = sprintf('AESamples-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
RunOptions.SaveFileNameAESamplesAccurate = sprintf('AESamplesAccurate-%s-%s-%sD-%dSensors-%sFinalTime-%dAESamples-%sAE',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDataDomain,RunOptions.TrelisMeshDElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
if RunOptions.UseACPrior == 1
    RunOptions.SaveFileNameAESamples = sprintf('AESamples-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE-ACPrior',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
    RunOptions.SaveFileNameAESamplesAccurate = sprintf('AESamplesAccurate-%s-%s-%sD-%dSensors-%sFinalTime-%dAESamples-%sAE-ACPrior',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDataDomain,RunOptions.TrelisMeshDElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
end
if RunOptions.FluidDomainWithSolidLayerMeshD == 1
    if RunOptions.AE_VarySolidLayer == 1
        RunOptions.SaveFileNameAESamples = sprintf('AESamples-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE-vSL',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
        RunOptions.SaveFileNameAESamplesAccurate = sprintf('AESamplesAccurate-%s-%s-%sD-%dSensors-%sFinalTime-%dAESamples-%sAE-vSL',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDataDomain,RunOptions.TrelisMeshDElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
    end
end
if RunOptions.AE_ConstructStats_SVD == 1;
    RunOptions.SaveFileNameAEStats = sprintf('AEStats-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
    RunOptions.SaveFileNameAEErrorModel = sprintf('AEErrorModel-%s-%s-%sD-%sI-%s%s-%dSensors-%sFinalTime-%dAESamples-%sAE',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
    if RunOptions.UseACPrior == 1
        RunOptions.SaveFileNameAEStats = sprintf('AEStats-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE-ACPrior',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
        RunOptions.SaveFileNameAEErrorModel = sprintf('AEErrorModel-%s-%s-%sD-%sI-%s%s-%dSensors-%sFinalTime-%dAESamples-%sAE-ACPrior',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
    end
    if RunOptions.FluidDomainWithSolidLayerMeshD == 1
        if RunOptions.AE_VarySolidLayer == 1
            RunOptions.SaveFileNameAEStats = sprintf('AEStats-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAE-vSL',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
            RunOptions.SaveFileNameAEErrorModel = sprintf('AEErrorModel-%s-%s-%sD-%sI-%s%s-%dSensors-%sFinalTime-%dAESamples-%sAE-vSL',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
        end
    end
end
if RunOptions.AE_ConstructStats_QR == 1;
    RunOptions.SaveFileNameAEStats = sprintf('AEStats-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAEQR',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
    RunOptions.SaveFileNameAEErrorModel = sprintf('AEErrorModel-%s-%s-%sD-%sI-%s%s-%dSensors-%sFinalTime-%dAESamples-%sAEQR',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
    if RunOptions.UseACPrior == 1
        RunOptions.SaveFileNameAEStats = sprintf('AEStats-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAEQR-ACPrior',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);        
        RunOptions.SaveFileNameAEErrorModel = sprintf('AEErrorModel-%s-%s-%sD-%sI-%s%s-%dSensors-%sFinalTime-%dAESamples-%sAEQR-ACPrior',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
    end
    if RunOptions.FluidDomainWithSolidLayerMeshD == 1
        if RunOptions.AE_VarySolidLayer == 1
            RunOptions.SaveFileNameAEStats = sprintf('AEStats-%s-%s-%sD-%sI-%dSensors-%sFinalTime-%dAESamples-%sAEQR-vSL',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
            RunOptions.SaveFileNameAEErrorModel = sprintf('AEErrorModel-%s-%s-%sD-%sI-%s%s-%dSensors-%sFinalTime-%dAESamples-%sAEQR-vSL',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
            if RunOptions.FluidDomainWithSolidLayerMeshI == 1;
                error('Have not coded varying solid layer draws for SLvsSL case yet')
            end
        end
    end
end

%% =======================================================================%
%                    Computing AE Sample Statistics
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Constructing Accurate Mesh %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RunOptions.GenerateMesh == 1
    MeshA = QPAT_TrelisMesh2DGenerator(RunOptions,'Trelis_QPATMesh_Accurate.inp',3,1);
else
    LoadFileNameMeshA = sprintf('Mesh-%s',RunOptions.TrelisMeshAElementSize);
    load(LoadFileNameMeshA);
    MeshA = Mesh;
    clear Mesh
end
[DGMMeshA,PrecomputedIntrplteObjectsA] = EWE_DGM2D_Setup(RunOptions,MeshA,RunOptions.FluidMeshD,RunOptions.SolidMeshD,RunOptions.FluidDomainWithSolidLayerMeshD);
SensorsA = EWE_DGM2D_ConstructSensorArray(RunOptions,DataVrblsWave.SensorCoords,DGMMeshA.VX,DGMMeshA.VY,DGMMeshA.EToV,MeshA.Bnd_ElmInd,MeshA.Bnd_NodeInd,DGMMeshA.pinfo,DGMMeshA.Norder);
DGMMeshA.pinfo = EWE_DGM2D_PrecomputeUpwindFluxPNonConf(RunOptions,DGMMeshA.pinfo,DGMMeshA.Norder,DGMMeshA.rho,DGMMeshA.lambda,DGMMeshA.mu,MeshA.DomainIndices,RunOptions.FluidDomainWithSolidLayerMeshD,RunOptions.SolidMeshD);
PLOT.TRI_DGMMeshA=delaunay(DGMMeshA.x,DGMMeshA.y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Constructing Samples %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,hIS_FEM,errorT11,errorT22,errorT12,errorvx,errorvy] = EWE_DGM2D_AE_ConstructSamples(RunOptions,MeshA.Nodes,MeshA.Elements,MeshA.Dimensns,MeshA.DomainIndices,MeshI.Nodes,MeshI.Elements,...
                                                                                               DGMMeshA.x,DGMMeshA.y,DGMMeshA.Np,DGMMeshA.K,PrecomputedIntrplteObjectsA,DGMMeshA.pinfo,DGMMeshA.Norder,DGMMeshA.lambda,DGMMeshA.rho,...
                                                                                               xI,yI,NpI,KI,rhoI,PrecomputedIntrplteObjectsI,pinfoI,...
                                                                                               SensorsA,SensorsI,dt,Prior,PLOT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Constructing Sample Statistics and Noise Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RunOptions.AE_ConstructStats_SVD == 1;
    [T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy] = EWE_DGM2D_AE_ConstructStats_SVD(RunOptions,size(SensorsI,2),...
                                                                                                        DataVrblsWave.ET11,DataVrblsWave.ET22,DataVrblsWave.ET12,DataVrblsWave.Evx,DataVrblsWave.Evy,...
                                                                                                        T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,...
                                                                                                        errorT11,errorT22,errorT12,errorvx,errorvy);
end
if RunOptions.AE_ConstructStats_QR == 1;
    [T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,L_ET11,L_ET22,L_ET12,L_Evx,L_Evy] = EWE_DGM2D_AE_ConstructStats_QR(RunOptions,size(SensorsI,2),hIS_FEM,...
                                                                                                                                       DataVrblsWave.ET11,DataVrblsWave.ET22,DataVrblsWave.ET12,DataVrblsWave.Evx,DataVrblsWave.Evy,...
                                                                                                                                       T11ErrorMean,T22ErrorMean,T12ErrorMean,vxErrorMean,vyErrorMean,... 
                                                                                                                                       errorT11,errorT22,errorT12,errorvx,errorvy);    
end
                                                                                                       
clear DGMMeshA PrecomputedIntrplteObjectsA SensorsA hAS_FEM hIS_FEM errorT11 errorT22 errorT12 errorvx errorvy