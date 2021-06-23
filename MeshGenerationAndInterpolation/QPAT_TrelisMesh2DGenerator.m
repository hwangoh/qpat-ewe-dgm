function [Mesh] = QPAT_TrelisMesh2DGenerator(RunOptions,TrelisMeshString,NumberofDomains,NumberofBoundaryConditions)

% QPAT_TrelisMesh2DGenerator defines the global mesh parameters
%
% Hwan Goh, University of Auckland, New Zealand - 25/01/2018

%% =======================================================================%
%                    Reading Trelis Generated Mesh
%=========================================================================%
visualization.initialization = false; %Don't want to use Timo's visualization
[Mesh.Elements, Mesh.Nodes(:,1), Mesh.Nodes(:,2), Mesh.N_Elm, Mesh.N_Nodes, Mesh.BCType, Mesh.Bnd_NodeInd, Mesh.DomainIndices] = ReadGrid_trelis(TrelisMeshString,NumberofDomains,NumberofBoundaryConditions,visualization);
clc
Mesh.Dimensns = (1/2)*[0.02 0.02]; %Height and width of mesh, may need to manually input this because it is not outputted in the Trelis file
clear visualization

%% =======================================================================%
%              Generating Boundary Elements and Light Sources
%=========================================================================%
%=== Boundary Elements ===%
Mesh.Bnd_ElmInd = find(sum(Mesh.BCType,2)>=1);
%=== Light Sources ===%
Mesh.Lght_ElmtInd = Mesh.Bnd_ElmInd;
Mesh.Lght_Nelm = size(Mesh.Lght_ElmtInd,1);
Mesh.N_Lght = 1; %Light sources are the whole boundary. This was required for the rectangular mesh to select which edge acts as a light source; so MeshD.N_Lght \in {1,2,3,4}.
Mesh.I_s = RunOptions.OptI_s;