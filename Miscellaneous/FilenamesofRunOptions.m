%=== Filenames of Outputs ===%
if RunOptions.FluidDomainWithSolidLayerMeshD == 1 && RunOptions.FluidMeshI == 1
    RunOptions.SaveFileNameDomain = 'SLvsF';
end
if RunOptions.FluidDomainWithSolidLayerMeshD == 1 && RunOptions.FluidDomainWithSolidLayerMeshI == 1
    RunOptions.SaveFileNameDomain = 'SLvsSL';
end
if RunOptions.FluidMeshD == 1 && RunOptions.FluidMeshI == 1
    RunOptions.SaveFileNameDomain = 'FvsF';
end
if RunOptions.FluidDomainWithSolidLayerMeshD == 1 %For saving data
    RunOptions.SaveFileNameDataDomain = 'SLDataDomain';
end
if RunOptions.FluidMeshD == 1 %For saving data
    RunOptions.SaveFileNameDataDomain = 'FDataDomain';
end
if RunOptions.EWE_LS_LogExpTransform == 1;
    RunOptions.SaveFileNamePosConstraint = 'LogExpT';
end
if RunOptions.EWE_LS_SigmoidTransform == 1;
    RunOptions.SaveFileNamePosConstraint = 'SigT';
end
if RunOptions.EWE_LS_ExponentialTransform == 1;
    RunOptions.SaveFileNamePosConstraint = 'ExpT';
end
if RunOptions.EWE_LS_LogExpTransform == 0 && RunOptions.EWE_LS_ExponentialTransform == 0 && RunOptions.EWE_LS_SigmoidTransform == 0;
    RunOptions.SaveFileNamePosConstraint = 'NoPosCon';
end

%=== AE Properties ===%
if RunOptions.AE_EWE == 0;
    RunOptions.AEN_Samples = 0;
    RunOptions.TrelisMeshAElementSize = 'No';
end

%=== Noisy Data Properties ===%
if RunOptions.AddNoise == 0;
    RunOptions.NoiseLevel = 0;
end
TempString = num2str(RunOptions.NoiseLevel);
RunOptions.SaveFileNameNoiseLevel = sscanf(TempString(3:end),'%s');
if RunOptions.NoiseMinMax == 1
    RunOptions.SaveFileNameNoiseType = 'MinMax';
end
if RunOptions.NoiseMinMaxS == 1
    RunOptions.SaveFileNameNoiseType = 'MinMaxS';
end
if RunOptions.NoiseMax == 1
    RunOptions.SaveFileNameNoiseType = 'Max';
end
if RunOptions.NoiseMaxS == 1
    RunOptions.SaveFileNameNoiseType = 'MaxS';
end

%=== Parameter Type and Final Time ===%
if RunOptions.TestAcousticParameters == 1;
    RunOptions.SaveFileNameParameterType = 'TestPrmtrs';
    TempString = num2str(RunOptions.FinalTime);
    RunOptions.SaveFileNameFinalTime = sscanf(TempString(3:end),'%s');
end
if RunOptions.RealAcousticParameters == 1;
    RunOptions.SaveFileNameParameterType = 'RealPrmtrs';
    TempString = num2str(RunOptions.FinalTime);
    RunOptions.SaveFileNameFinalTime = sscanf(TempString(1),'%s');
    clear TempString
end

%=== Save File Name ===%
RunOptions.SaveFileName = sprintf('%s-%s-%sD-%sI-%s%s-%dSensors-%sFinalTime-%s-%dAESamples-%sAE',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.SaveFileNamePosConstraint,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
if RunOptions.UseACPrior == 1
    RunOptions.SaveFileName = sprintf('%s-%s-%sD-%sI-%s%s-%dSensors-%sFinalTime-%s-%dAESamples-%sAE-ACPrior',RunOptions.SaveFileNameParameterType,RunOptions.SaveFileNameDomain,RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.SaveFileNameNoiseLevel,RunOptions.SaveFileNameNoiseType,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime,RunOptions.SaveFileNamePosConstraint,RunOptions.AEN_Samples,RunOptions.TrelisMeshAElementSize);
end

