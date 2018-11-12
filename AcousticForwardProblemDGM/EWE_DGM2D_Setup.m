function [DGMMesh,PrecomputedIntrplteObjects] = EWE_DGM2D_Setup(RunOptions,Mesh,FluidMesh,SolidMesh,FluidDomainWithSolidLayerMesh)

% EWEDGM2DSetup generates the required structures for the
% discontinuous Galerkin method which are independent of the illumination
% pattern
%
% Hwan Goh, University of Auckland, New Zealand - 20/7/2015

%% =======================================================================%
%                       Setting Up the DGM Grid
%=========================================================================%
Nfaces = 3; %Number of faces, required for constructing BCType in the script CorrectBCTable

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Construct DGM Mesh %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
[DGMMesh,PrecomputedIntrplteObjects] = ConstructDGMMesh(RunOptions,Mesh); %No need to output objects since they are saved in Globals2D

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Domain Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Fluid Domain ===%
if FluidMesh == 1 || FluidDomainWithSolidLayerMesh == 1
    DGMMesh.rho = RunOptions.AcousticMediumDensity*ones(DGMMesh.Np,DGMMesh.K);  %Acoustic medium density for edges of all elements [kg/m^3]
    DGMMesh.mu = zeros(DGMMesh.Np,DGMMesh.K); %Lame's second parameter for edges of all elements [kg/m^3]
end

%=== Solid Domain ===%
if SolidMesh == 1
    DGMMesh.rho = RunOptions.ElasticMediumDensity*ones(DGMMesh.Np,DGMMesh.K);  %Elastic medium density for edges of all elements [kg/m^3]
    DGMMesh.mu = RunOptions.ElasticLamemu*ones(DGMMesh.Np,DGMMesh.K);
end
DGMMesh.lambda = RunOptions.ElasticLamelambda*ones(DGMMesh.Np,DGMMesh.K);  %Lame's first parameter for edges of all elements [kg/m^3]

%=== Solid Layer Domain ===%
if RunOptions.UseTrelisMesh == 1 && FluidDomainWithSolidLayerMesh == 1
    for ii=1:Mesh.N_Elm
        if Mesh.DomainIndices(ii) == 1 || Mesh.DomainIndices(ii) == 3
            DGMMesh.mu(:,ii) = 0;
            DGMMesh.rho(:,ii) = RunOptions.AcousticMediumDensity;
        else
            DGMMesh.mu(:,ii) = RunOptions.ElasticLamemu;
            DGMMesh.rho(:,ii) = RunOptions.ElasticMediumDensity;
        end
    end
end
