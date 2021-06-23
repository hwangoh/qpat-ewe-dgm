function Trelis_SolidLayerMesh2D_Rectangular(GenerateDataMesh,GenerateInverseMesh,GenerateAccurateMesh,OuterSoftTissue,SkullLayer,InnerSoftTissue)
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

%=== Create The Main Rectangle ===%
fprintf(fidout,'%s\n',sprintf('create surface rectangle width 0.02 height 0.02 zplane'));


%% =======================================================================%
%            Generating Soft Tissue Layers and Skull Layers
%=========================================================================%
%=== Creating Skin Layer ===%
% fprintf(fidout,'%s\n',sprintf('create brick width 0.01929 depth 0.01929 height 0.01929'));
fprintf(fidout,'%s\n',sprintf('create brick width 0.01800 depth 0.01800 height 0.01929'));
fprintf(fidout,'%s\n',sprintf('create sheet extended from surface 4 5 6 7'));
fprintf(fidout,'%s\n',sprintf('delete volume 2'));
fprintf(fidout,'%s\n',sprintf('webcut body 1 with sheet body 3'));
fprintf(fidout,'%s\n',sprintf('delete volume 3'));

%=== Creating Skull Layer ===%
% fprintf(fidout,'%s\n',sprintf('create brick width 0.0175 depth 0.0175 height 0.018'));
fprintf(fidout,'%s\n',sprintf('create brick width 0.015 depth 0.015 height 0.018'));
fprintf(fidout,'%s\n',sprintf('create sheet extended from surface 16 17 18 19'));
fprintf(fidout,'%s\n',sprintf('delete volume 5'));
fprintf(fidout,'%s\n',sprintf('webcut body 1 with sheet body 6'));
fprintf(fidout,'%s\n',sprintf('delete volume 6'));

%==== Merge Split Subdomains ====% (note that interface remains)
fprintf(fidout,'%s\n','merge vol all');

%=== Meshing Surfaces ===%
% first generate grid to surface id 2 (here we set the mesh size)
fprintf(fidout,'%s\n',sprintf('surface 13 size %f', OuterSoftTissue));
% this line defines the meshing scheme
fprintf(fidout,'%s\n',sprintf('surface 13 scheme trimesh'));
% here the actual grid is generated
fprintf(fidout,'%s\n',sprintf('mesh surface 13'));

% idea same as for the surface 2 above
fprintf(fidout,'%s\n',sprintf('surface 25 size %f', SkullLayer));
fprintf(fidout,'%s\n',sprintf('surface 25 scheme trimesh'));
fprintf(fidout,'%s\n',sprintf('mesh surface 25'));

% idea same as for the surface 2 above
fprintf(fidout,'%s\n',sprintf('surface 24 size %f', InnerSoftTissue));
fprintf(fidout,'%s\n',sprintf('surface 24 scheme trimesh'));
fprintf(fidout,'%s\n',sprintf('mesh surface 24'));

%=== Element Sets ===%
fprintf(fidout,'%s\n',sprintf('Block 1 surface 13'));
fprintf(fidout,'%s\n',sprintf('Block 2 surface 25'));
fprintf(fidout,'%s\n',sprintf('Block 3 surface 24'));

%=== Boundaries ===%
fprintf(fidout,'%s\n',sprintf('SideSet 1 Curve 1 2 3 4'));

%=== Boundary Nodes ===%
fprintf(fidout,'%s\n',sprintf('NodeSet 1 Curve 1 2 3 4'));

%% =======================================================================%
%                           Output Mesh
%=========================================================================%
% export to abaqus format (there might be a better way to do this - should
% be checked in some stage)
if GenerateDataMesh == 1
%     fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/Hwan/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Data.inp" dimension 2 block all overwrite'));
    fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/hgoh009/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Data.inp" dimension 2 block all overwrite'));
end
if GenerateInverseMesh == 1
%     fprintf(fidout,'%s\n',sprintf('export Abaqus "C:/Users/Hwan/Dropbox/Work/PhD/PhD_Stuff/Codes/QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation/Trelis_QPATMesh_Inverse.inp" dimension 2 block all overwrite'));
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



