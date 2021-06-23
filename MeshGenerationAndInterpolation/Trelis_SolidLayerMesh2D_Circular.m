function Trelis_SolidLayerMesh2D_Circular(GenerateDataMesh,GenerateInverseMesh,GenerateAccurateMesh,OuterSoftTissue,SkullLayer,InnerSoftTissue)
% this file generates input file to Trelis/Cubit software

% open file
if GenerateDataMesh == 1
    fidout = fopen('../MeshGenerationAndInterpolation/Trelis_QPATMesh_Data.jou', 'w');
end

if GenerateInverseMesh == 1
    fidout = fopen('../MeshGenerationAndInterpolation/Trelis_QPATMesh_Inverse.jou', 'w');
end
if GenerateAccurateMesh == 1
    fidout = fopen('../MeshGenerationAndInterpolation/Trelis_QPATMesh_Accurate.jou', 'w');
end

%% =======================================================================%
%                          Generating Domain
%=========================================================================%
%=== Initialize Cubit ===%
fprintf(fidout,'%s\n','reset');

%=== Create The Main Cylinder ===%
fprintf(fidout,'%s\n',sprintf('create surface circle radius 0.01 zplane'));

%% =======================================================================%
%            Generating Soft Tissue Layers and Skull Layers
%=========================================================================%
%=== Creating Layers ===%
fprintf(fidout,'%s\n',sprintf('webcut body 1 with cylinder radius 0.0098 axis z center 0 0 0'));
fprintf(fidout,'%s\n',sprintf('webcut body 2 with cylinder radius 0.009 axis z center 0 0 0'));

%==== Merge Split Subdomains ====% (note that interface remains)
fprintf(fidout,'%s\n','merge vol all');

%=== Meshing Surfaces ===%
% first generate grid to surface id 2 (here we set the mesh size)
fprintf(fidout,'%s\n',sprintf('surface 2 size %f', OuterSoftTissue));
% this line defines the meshing scheme
fprintf(fidout,'%s\n',sprintf('surface 2 scheme trimesh'));
% here the actual grid is generated
fprintf(fidout,'%s\n',sprintf('mesh surface 2'));

% idea same as for the surface 2 above
fprintf(fidout,'%s\n',sprintf('surface 4 size %f', SkullLayer));
fprintf(fidout,'%s\n',sprintf('surface 4 scheme trimesh'));
fprintf(fidout,'%s\n',sprintf('mesh surface 4'));

% idea same as for the surface 2 above
fprintf(fidout,'%s\n',sprintf('surface 5 size %f', InnerSoftTissue));
fprintf(fidout,'%s\n',sprintf('surface 5 scheme trimesh'));
fprintf(fidout,'%s\n',sprintf('mesh surface 5'));

%=== Element Sets ===%
fprintf(fidout,'%s\n',sprintf('Block 1 surface 2'));
fprintf(fidout,'%s\n',sprintf('Block 2 surface 4'));
fprintf(fidout,'%s\n',sprintf('Block 3 surface 5'));

%=== Boundaries ===%
%Outer Boundary
fprintf(fidout,'%s\n',sprintf('SideSet 1 Curve 1'));

%% =======================================================================%
%                           Output Mesh
%=========================================================================%
% export to abaqus format (there might be a better way to do this - should
% be checked in some stage)
if GenerateDataMesh == 1
    % fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/Hwan/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Data.inp" dimension 2 block all overwrite'));
    fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/hgoh009/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Data.inp" dimension 2 block all overwrite'));
end
if GenerateInverseMesh == 1
    % fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/Hwan/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Inverse.inp" dimension 2 block all overwrite'));
    fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/hgoh009/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Inverse.inp" dimension 2 block all overwrite'));
end
if GenerateAccurateMesh == 1
%     fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/Hwan/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Accurate.inp" dimension 2 block all overwrite'));
    fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/hgoh009/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Accurate.inp" dimension 2 block all overwrite'));
end

% close the file
fclose(fidout);

% call the software (note that this works only in linux and perhaps in ma
% but surely not in windows machines)
fprintf(' Generating mesh ...\n');
tic
    % call trelis in cubit folder 
%    cd cubit    
%    system('/opt/Trelis-15.2/bin/trelis -nographics -batch cubit.jou > meshgenerator.txt');
%    cd ..
toc
% please check meshgenerator.txt for possible errors (there should not be
% any)


return



