% QPAT_RectangularMesh2DGenerator defines the global mesh parameters
%
% Hwan Goh, University of Auckland, New Zealand - 6/7/2015

%% =======================================================================%
%                               Data Mesh
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data Mesh Domain Properties %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MeshD.dx = (0.005/MeshD.N_BndElmx); %spacing between grid points in the x direction [metres] Default = 0.005/7
MeshD.dy = (0.005/MeshD.N_BndElmy); %spacing between grid points in the y direction [metres] Default = 0.005/7
MeshD.ExtndBndx = RunOptions.BndExtndx; %total extension of boundary in +/-x direction [metres] default=0.02
MeshD.ExtndBndy = RunOptions.BndExtndy; %total extension of boundary in +/-y direction [metres] default=0.02
if (RunOptions.BndExtndx == 0) && (RunOptions.BndExtndy == 0)
    MeshD.Nx = ((2*MeshD.Dimensns(1)+MeshD.ExtndBndx)/(RunOptions.ScalingOptical*MeshD.dx))+1; %number of grid points in the x direction
    MeshD.Ny = ((2*MeshD.Dimensns(2)+MeshD.ExtndBndy)/(RunOptions.ScalingOptical*MeshD.dy))+1; %number of grid points in the y direction
else
    MeshD.Nx = (2*MeshD.Dimensns(1)+MeshD.ExtndBndx)/(RunOptions.ScalingOptical*MeshD.dx); %number of grid points in the x direction
    MeshD.Ny = (2*MeshD.Dimensns(2)+MeshD.ExtndBndy)/(RunOptions.ScalingOptical*MeshD.dy); %number of grid points in the y direction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generating the Mesh %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MeshD]=TriangulateRect(RunOptions,MeshD);

%% =======================================================================%
%                             Inverse Mesh
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inverse Mesh Domain Properties %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MeshI.dx = (0.005/MeshI.N_BndElmx); %spacing between grid points in the x direction [metres] Default = 0.005/7
MeshI.dy = (0.005/MeshI.N_BndElmy); %spacing between grid points in the y direction [metres] Default = 0.005/7
MeshI.ExtndBndx = RunOptions.BndExtndx; %total extension of boundary in +/-x direction [metres] default=20
MeshI.ExtndBndy = RunOptions.BndExtndy; %total extension of boundary in +/-y direction [metres] default=20
if (RunOptions.BndExtndx == 0) && (RunOptions.BndExtndy == 0)
    MeshI.Nx = ((2*MeshI.Dimensns(1)+MeshI.ExtndBndx)/(RunOptions.ScalingOptical*MeshI.dx))+1; %number of grid points in the x direction
    MeshI.Ny = ((2*MeshI.Dimensns(2)+MeshI.ExtndBndy)/(RunOptions.ScalingOptical*MeshI.dy))+1; %number of grid points in the y direction
else
    MeshI.Nx = (2*MeshI.Dimensns(1)+MeshI.ExtndBndx)/(RunOptions.ScalingOptical*MeshI.dx); %number of grid points in the x direction
    MeshI.Ny = (2*MeshI.Dimensns(2)+MeshI.ExtndBndy)/(RunOptions.ScalingOptical*MeshI.dy); %number of grid points in the y direction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generating the Mesh %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MeshI]=TriangulateRect(RunOptions,MeshI);
