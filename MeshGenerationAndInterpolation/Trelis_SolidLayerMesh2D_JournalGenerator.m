restoredefaultpath
addpath('QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation')

%=== For Thesis ===%
% 0.00028 - Elements: 10938, Nodes: 5612
% 0.00034 - Elements: 7427,  Nodes: 3831
% 0.00045 - Elements: 3862,  Nodes: 2020
% 0.00055 - Elements: 2802,  Nodes: 1474
% 0.0007  - Elements: 1530,  Nodes: 824
% 0.00082 - Elements: 1200,  Nodes: 649

% 0.0003ASYM - Elements: 9612, Nodes: 4941
% 0.00034ASYM - Elements: 7562, Nodes: 3900
% 0.00034ASYM - Elements: 7076, Nodes: 3653
% 0.00036ASYM - Elements: 6718, Nodes: 3472
% 0.00037ASYM - Elements: 6310, Nodes: 3264
% 0.00038ASYM - Elements: 5964, Nodes: 3089
% 0.0004ASYM - Elements: 5466, Nodes: 2834
% 0.00045ASYM - Elements: 4096, Nodes: 2137
% 0.00055ASYM - Elements: 2764, Nodes: 1455
% 0.0007ASYM - Elements: 1758, Nodes: 938
% 0.0009ASYM - Elements: 1058, Nodes: 574

%=== For Debugging ===%
% 0.001 -  Elements: 760,  Nodes: 421
% 0.002 -  Elements: 284,  Nodes: 145
% 0.001ASYM -  Elements: 856,  Nodes: 469
% 0.002ASYM -  Elements: 226,  Nodes: 134

close all
clear all
clc

%=== Choose Mesh ===%
GenerateDataMesh = 1;
GenerateInverseMesh = 0;
GenerateAccurateMesh = 0;

%=== Size of Elements In Each Domain ===%
if GenerateDataMesh == 1
    OuterSoftTissue = 0.0005;
    SkullLayer = 0.0005;
    InnerSoftTissue = 0.0005;
end
if GenerateInverseMesh == 1
    OuterSoftTissue = 0.00045;
    SkullLayer = 0.00045;
    InnerSoftTissue = 0.00045;
end
if GenerateAccurateMesh == 1
    OuterSoftTissue = 0.00045;
    SkullLayer = 0.00045;
    InnerSoftTissue = 0.00045;
end

%=== Generate Journal File ===%
% Trelis_SolidLayerMesh2D_Rectangular(GenerateDataMesh,GenerateInverseMesh,GenerateAccurateMesh,OuterSoftTissue,SkullLayer,InnerSoftTissue);
Trelis_SolidLayerMesh2D_Rectangular_Asymmetric(GenerateDataMesh,GenerateInverseMesh,GenerateAccurateMesh,OuterSoftTissue,SkullLayer,InnerSoftTissue);
% Trelis_SolidLayerMesh2D_Circular(GenerateDataMesh,GenerateInverseMesh,GenerateAccurateMesh,OuterSoftTissue,SkullLayer,InnerSoftTissue);
