function [rhs_e11,rhs_e22,rhs_e12,rhs_vx,rhs_vy]=EWE_DGM2D_PNonConf_Strain(RunOptions,pinfo,Np,K,e11,e22,e12,vx,vy)

% EWE_DGM2D_PNonConf_Strain.m is the main script file for the discontinuous Galerkin method
% implementation of the acoustic forward problem using the elastic wave
% equation in strain-velocity form. This code is enabled for p-nonconforming meshes and is written in a
% style more similar to Timo Lahivaara's code and also
% utilizes codes from the textbook "Nodal Discontinous Galerkin Methods, Algorithms, Analysis and
% Applications" by Jan Hesthaven and Tim Warburton, 2007
%
% Inputs: 
%   RunOptions:
%          BC - Set boundary conditions
%   DataVrbls:
%          rho_e - Mass density per unit volume
%          mu_e - First Lame parameter, also known as the shear modulus. Setting the shear modulus to zero yields the acoustic fluid wave equation
%          lambda_e = Second Lame parameter
%   pinfo: structural data with respect to the differing element orders
%   e11: Strain Tensor entry
%   e12: Strain Tensor entry
%   e22: Strain Tensor entry
%   vx: x component of the velocity of the displacement
%   vy: y component of the velocity of the displacement
%          
%  Outputs:
%   rhs_e11: Velocity of strain tensor component e11
%   rhs_e22: Velocity of strain tensor component e22
%   rhs_e12: Velocity of strain tensor component e12
%   rhs_vx: Velocity of the x component of the velocity of the displacement
%   rhs_vy: Velocity of the y component of the velocity of the displacement
%
% Notes: The textbook's version of BuildPNonCon2D creates indexing which,
%        for only some elements, differs from Timo's version of BuildPNonCon2D 
%        which I have called BuildPNonCon2DTimo.
%
% Hwan Goh, 2/01/2017 (Happy New Year!), University of Auckland, New Zealand
% Last Edited: 16/11/2017 - Removing reliance on Globals2D

Nfaces = 3;

%% %%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Structures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== rhs Storage ===% %Needs to be revised! Np and K defined before Np=pinf.Np and K=pinf.K is only fine for a uniform mesh
rhs_e11 = zeros(Np*K,1);
rhs_e22 = zeros(Np*K,1);
rhs_e12 = zeros(Np*K,1);
rhs_vx   = zeros(Np*K,1);
rhs_vy  = zeros(Np*K,1);

for N1=1:length(pinfo)
    pinf = pinfo(N1);
    Np = pinf.Np;
    K = pinf.K;
    Nfp = pinf.Nfp;
    mapO = pinf.mapO;
    mapD = pinf.mapD;
    mapN = pinf.mapN;
    if (K>0)
        ids = pinf.ids;
        Fmask = pinf.Fmask;
        e11N = e11(ids);
        e22N = e22(ids);
        e12N = e12(ids);
        vxN = vx(ids);
        vyN = vy(ids);
        %=== Interior and Exterior Terms ===%
        e11M = e11N(Fmask(:),:); 
        e22M = e22N(Fmask(:),:);
        e12M = e12N(Fmask(:),:);
        vxM   = vxN(Fmask(:),:);
        vyM   = vyN(Fmask(:),:);
        e11P = zeros(Nfaces*Nfp,K);
        e22P = zeros(Nfaces*Nfp,K);
        e12P = zeros(Nfaces*Nfp,K);
        vxP = zeros(Nfaces*Nfp,K);
        vyP = zeros(Nfaces*Nfp,K);
        for N2 = 1:length(pinfo)
            if ~isempty(pinf.fmapM{N2})
                interpMAT = pinf.interpP{N2};
                fmapM = pinf.fmapM{N2};
                vmapP = pinf.vmapP{N2};
                e11P(fmapM) = interpMAT*e11(vmapP);
                e12P(fmapM) = interpMAT*e12(vmapP);
                e22P(fmapM) = interpMAT*e22(vmapP);
                vxP(fmapM)   = interpMAT*vx(vmapP);
                vyP(fmapM)   = interpMAT*vy(vmapP);
            end
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Set Boundary Conditions %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bcvece11 = ones(Nfaces*Nfp,K);
        bcvece22 = ones(Nfaces*Nfp,K);
        bcvece12 = ones(Nfaces*Nfp,K);
        bcvecvx = ones(Nfaces*Nfp,K);
        bcvecvy = ones(Nfaces*Nfp,K);

        %Out/Absorbing
        bcvece11(mapO) = 0;
        bcvece22(mapO) = 0;
        bcvece12(mapO) = 0;
        bcvecvx(mapO) = 0;
        bcvecvy(mapO) = 0;             
        %Dirichlet
        bcvece11(mapD) = 1;
        bcvece22(mapD) = 1;
        bcvece12(mapD) = 1;
        bcvecvx(mapD) = -1;
        bcvecvy(mapD) = -1;               
        %Neumann
        bcvece11(mapN) = -1;
        bcvece22(mapN) = -1;
        bcvece12(mapN) = 1;
        bcvecvx(mapN) = 1;
        bcvecvy(mapN) = 1;
               
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Some Precomputed Values %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % _T_ means .* and _P_ means + and _D_ means ./
        %=== Sums of Strain Values ===%
        e11M_P_e22M = e11M+e22M;
        bcvece11_T_e11P_P_bcvece22_T_e22P = bcvece11.*e11P+bcvece22.*e22P;
        %=== Sums of Wave Speeds and Medium Densities ===%
        rhoP_T_csM = pinf.rhoP.*pinf.csM;
        %=== [[v_e]] ===%
        ndv  = pinf.nx.*vxM+pinf.ny.*vyM-(pinf.nx.*bcvecvx.*vxP+pinf.ny.*bcvecvy.*vyP);
        %=== [v_e] ===%
        dvx = vxM-bcvecvx.*vxP;
        dvy = vyM-bcvecvy.*vyP;
        %=== [[S]] ===%
        sn1 = 2*pinf.muM.*(e11M.*pinf.nx+pinf.ny.*e12M)+pinf.lambdaM.*e11M_P_e22M.*pinf.nx-...
             (2*pinf.muP.*(bcvece11.*e11P.*pinf.nx+pinf.ny.*bcvece12.*e12P)+pinf.lambdaP.*bcvece11_T_e11P_P_bcvece22_T_e22P.*pinf.nx); %First row of [[S]]
        sn2 = 2*pinf.muM.*(e12M.*pinf.nx+pinf.ny.*e22M)+ pinf.lambdaM.*e11M_P_e22M.*pinf.ny-...
             (2*pinf.muP.*(bcvece12.*e12P.*pinf.nx+pinf.ny.*bcvece22.*e22P)+pinf.lambdaP.*bcvece11_T_e11P_P_bcvece22_T_e22P.*pinf.ny); %Second row of [[S]]
        %=== n'[[S]] ===%
        ndsn = 2*pinf.muM.*(pinf.nxnx.*e11M + 2*pinf.nxny.*e12M + pinf.nyny.*e22M)+pinf.lambdaM.*e11M_P_e22M-...
              (2*pinf.muP.*(pinf.nxnx.*bcvece11.*e11P + 2*pinf.nxny.*bcvece12.*e12P + pinf.nyny.*bcvece22.*e22P) + pinf.lambdaP.*bcvece11_T_e11P_P_bcvece22_T_e22P);
        %=== Some other values for the flux ===%
        ndsn_T_nx = ndsn.*pinf.nx;
        ndv_T_nx = ndv.*pinf.nx;
        ndsn_T_ny = ndsn.*pinf.ny;
        ndv_T_ny = ndv.*pinf.ny;
        %=== B1cpM ===%
        B_1cpM = pinf.cpM.*(pinf.d11.*ndsn +  pinf.d13.*ndv);
       
        %% %%%%%%%%%
        %%% Flux %%%
        %%%%%%%%%%%%
        flux = zeros(Nfaces*Nfp, K, 5);
        flux(:,:,1) = B_1cpM.*pinf.r1(:,:,1) + pinf.B_2csM.*(pinf.nx.*(sn1-ndsn_T_nx)                                           + rhoP_T_csM.*pinf.nx.*(dvx-ndv_T_nx));
        flux(:,:,2) = B_1cpM.*pinf.r1(:,:,2) + pinf.B_2csM.*(pinf.ny.*(sn2-ndsn_T_ny)                                           + rhoP_T_csM.*pinf.ny.*(dvy-ndv_T_ny));
        flux(:,:,3) = B_1cpM.*pinf.r1(:,:,3) + pinf.B_2csM.*(pinf.nx/2.*(sn2-ndsn_T_ny) + pinf.ny/2.*(sn1-ndsn_T_nx)            + rhoP_T_csM.*(pinf.nx/2.*(dvy-ndv_T_ny) + pinf.ny/2.*(dvx-ndv_T_nx)));
        flux(:,:,4) = B_1cpM.*pinf.r1(:,:,4) + pinf.B_2csM.*(pinf.csM.*(sn1-ndsn_T_nx)                                          + rhoP_T_csM.*pinf.csM.*(dvx-ndv_T_nx));
        flux(:,:,5) = B_1cpM.*pinf.r1(:,:,5) + pinf.B_2csM.*(pinf.csM.*(sn2-ndsn_T_ny)                                          + rhoP_T_csM.*pinf.csM.*(dvy-ndv_T_ny));        

        %% %%%%%%%%%%%%%%%%%%%%%
        %%% Derivative Terms %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate local derivatives of fields
        de11dr = pinf.Dr*e11N; de11ds = pinf.Ds*e11N;
        de12dr = pinf.Dr*e12N; de12ds = pinf.Ds*e12N;
        de22dr = pinf.Dr*e22N; de22ds = pinf.Ds*e22N;
        dvxdr  = pinf.Dr*vxN;   dvxds = pinf.Ds*vxN;
        dvydr  = pinf.Dr*vyN;   dvyds = pinf.Ds*vyN;
        
        % Compute physical derivatives of fields
        de11dy = pinf.ry.*de11dr + pinf.sy.*de11ds;
        de11dx = pinf.rx.*de11dr + pinf.sx.*de11ds;
        de12dy = pinf.ry.*de12dr + pinf.sy.*de12ds;
        de12dx = pinf.rx.*de12dr + pinf.sx.*de12ds;
        de22dy = pinf.ry.*de22dr + pinf.sy.*de22ds;
        de22dx = pinf.rx.*de22dr + pinf.sx.*de22ds;
        dvxdy  = pinf.ry.*dvxdr  + pinf.sy.*dvxds;
        dvxdx  = pinf.rx.*dvxdr  + pinf.sx.*dvxds;
        dvydy  = pinf.ry.*dvydr  + pinf.sy.*dvyds;
        dvydx  = pinf.rx.*dvydr  + pinf.sx.*dvyds;
        
        %% %%%%%%%%%%%%%%%%%%
        %%% Compute qvelo %%%
        %%%%%%%%%%%%%%%%%%%%%
        %Defining matrix A and B entries
        A41_D_rho = pinf.A41_D_rho;
        A42_D_rho = pinf.A42_D_rho;
        B43_D_rho = pinf.B43_D_rho;
        A53_D_rho = B43_D_rho;
        B51_D_rho = A42_D_rho;
        B52_D_rho = A41_D_rho;
        
        rhs_e11(ids) = dvxdx - pinf.LIFT*(pinf.Fscale.*flux(:,:,1));
        rhs_e22(ids) = dvydy - pinf.LIFT*(pinf.Fscale.*flux(:,:,2));
        rhs_e12(ids) = (1/2)*dvxdy + (1/2)*dvydx - pinf.LIFT*(pinf.Fscale.*flux(:,:,3));
        rhs_vx(ids)  = A41_D_rho.*de11dx + A42_D_rho.*de22dx + B43_D_rho.*de12dy - pinf.LIFT*(pinf.Fscale.*flux(:,:,4));
        rhs_vy(ids)  = B51_D_rho.*de11dy + B52_D_rho.*de22dy + A53_D_rho.*de12dx - pinf.LIFT*(pinf.Fscale.*flux(:,:,5));
    end
end