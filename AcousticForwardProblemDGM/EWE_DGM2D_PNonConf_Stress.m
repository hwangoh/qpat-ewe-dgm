function [rhs_s11,rhs_s22,rhs_s12,rhs_vx,rhs_vy]=EWE_DGM2D_PNonConf_Stress(RunOptions,pinfo,Np,K,s11,s22,s12,vx,vy)

% EWE_DGM2D_PNonConf_Stress.m is the main script file for the discontinuous Galerkin method
% implementation of the acoustic forward problem using the elastic wave
% equation in stress-velocity form. This code is enabled for p-nonconforming meshes and is written in a
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
%   s11: Stress Tensor entry
%   s12: Stress Tensor entry
%   s22: Stress Tensor entry
%   vx: x component of the velocity of the displacement
%   vy: y component of the velocity of the displacement
%          
%  Outputs:
%   rhs_s11: Velocity of stress tensor component e11
%   rhs_s22: Velocity of stress tensor component e22
%   rhs_s12: Velocity of stress tensor component e12
%   rhs_vx: Velocity of the x component of the velocity of the displacement
%   rhs_vy: Velocity of the y component of the velocity of the displacement
%
% Notes: The textbook's version of BuildPNonCon2D creates indexing which,
%        for only some elements, differs from Timo's version of BuildPNonCon2D 
%        which I have called BuildPNonCon2DTimo.
%
% Hwan Goh, 27/07/2017, University of Auckland, New Zealand

Nfaces = 3;

%% %%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Structures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== rhs Storage ===% %Needs to be revised! K defined before K=pinf.K is only fine for a uniform mesh
rhs_s11 = zeros(Np*K,1);
rhs_s22 = zeros(Np*K,1);
rhs_s12 = zeros(Np*K,1);
rhs_vx   = zeros(Np*K,1);
rhs_vy  = zeros(Np*K,1);

for N1=1:length(pinfo)
    pinf = pinfo(N1);
    Np = pinf.Np;
    K = pinf.K;
    Nfp = pinf.Nfp;
    if (K>0)
        ids = pinf.ids;
        Fmask = pinf.Fmask;
        s11N = s11(ids);
        s22N = s22(ids);
        s12N = s12(ids);
        vxN = vx(ids);
        vyN = vy(ids);
        %=== Interior and Exterior Terms ===%
        s11M = s11N(Fmask(:),:); 
        s22M = s22N(Fmask(:),:);
        s12M = s12N(Fmask(:),:);
        vxM   = vxN(Fmask(:),:);
        vyM   = vyN(Fmask(:),:);
        s11P = zeros(Nfaces*Nfp,K);
        s22P = zeros(Nfaces*Nfp,K);
        s12P = zeros(Nfaces*Nfp,K);
        vxP = zeros(Nfaces*Nfp,K);
        vyP = zeros(Nfaces*Nfp,K);
        for N2 = 1:length(pinfo)
            if ~isempty(pinf.fmapM{N2})
                interpMAT = pinf.interpP{N2};
                fmapM = pinf.fmapM{N2};
                vmapP = pinf.vmapP{N2};
                s11P(fmapM) = interpMAT*s11(vmapP);
                s12P(fmapM) = interpMAT*s12(vmapP);
                s22P(fmapM) = interpMAT*s22(vmapP);
                vxP(fmapM)   = interpMAT*vx(vmapP);
                vyP(fmapM)   = interpMAT*vy(vmapP);
            end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Set Boundary Conditions %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bcvecs11 = ones(Nfaces*Nfp,K);
        bcvecs22 = ones(Nfaces*Nfp,K);
        bcvecs12 = ones(Nfaces*Nfp,K);
        bcvecvx = ones(Nfaces*Nfp,K);
        bcvecvy = ones(Nfaces*Nfp,K);
        
        %Out/Absorbing
        bcvecs11(mapO) = 0;
        bcvecs22(mapO) = 0;
        bcvecs12(mapO) = 0;
        bcvecvx(mapO) = 0;
        bcvecvy(mapO) = 0;       
        %Dirichlet
        bcvecs11(mapD) = 1;
        bcvecs22(mapD) = 1;
        bcvecs12(mapD) = 1;
        bcvecvx(mapD) = -1;
        bcvecvy(mapD) = -1;        
        %Neumann
        bcvecs11(mapN) = -1;
        bcvecs22(mapN) = -1;
        bcvecs12(mapN) = 1;
        bcvecvx(mapN) = 1;
        bcvecvy(mapN) = 1;
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Some Precomputed Values %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % _T_ means .* and _P_ means + and _D_ means ./      
        %=== Sums of Wave Speeds and Medium Densities ===%
        csP_T_csM = pinf.csP.*pinf.csM;
        muP_T_csM = pinf.muP.*pinf.csM;
        %=== [[v_e]] ===%
        ndv  = pinf.nx.*vxM+pinf.ny.*vyM-(pinf.nx.*bcvecvx.*vxP+pinf.ny.*bcvecvy.*vyP);
        %=== [v_e] ===%
        dvx = vxM-bcvecvx.*vxP;
        dvy = vyM-bcvecvy.*vyP;
        %=== [[S]] ===%
        sn1 = s11M.*pinf.nx + s12M.*pinf.ny + bcvecs11.*s11P.*-pinf.nx + bcvecs12.*s12P.*-pinf.ny;
        sn2 = s12M.*pinf.nx + s22M.*pinf.ny + bcvecs12.*s12P.*-pinf.nx + bcvecs22.*s22P.*-pinf.ny;
        %=== n'[[S]] ===%
        ndsn = pinf.nx.*sn1 + pinf.ny.*sn2;
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
        if RunOptions.UseAcousticViaElasticDGMFormulation == 1;
            flux(:,:,1) = B_1cpM.*pinf.r1(:,:,1);
            flux(:,:,2) = B_1cpM.*pinf.r1(:,:,2);
            flux(:,:,3) = B_1cpM.*pinf.r1(:,:,3);
            flux(:,:,4) = B_1cpM.*pinf.r1(:,:,4);
            flux(:,:,5) = B_1cpM.*pinf.r1(:,:,5);
        else
            flux(:,:,1) = B_1cpM.*pinf.r1(:,:,1) + pinf.B_2csM.*(pinf.csP.*(2*pinf.muM).*pinf.nx.*(sn1-ndsn_T_nx) + pinf.muP.*(2*pinf.muM).*pinf.nx.*(dvx-ndv_T_nx));
            flux(:,:,2) = B_1cpM.*pinf.r1(:,:,2) + pinf.B_2csM.*(pinf.csP.*(2*pinf.muM).*pinf.ny.*(sn2-ndsn_T_ny) + pinf.muP.*(2*pinf.muM).*pinf.ny.*(dvy-ndv_T_ny));
            flux(:,:,3) = B_1cpM.*pinf.r1(:,:,3) + pinf.B_2csM.*(pinf.csP.*(2*pinf.muM).*(pinf.nx/2.*(sn2-ndsn_T_ny) + pinf.ny/2.*(sn1-ndsn_T_nx)) + pinf.muP.*(2*pinf.muM).*(pinf.nx/2.*(dvy-ndv_T_ny) + pinf.ny/2.*(dvx-ndv_T_nx)));
            flux(:,:,4) = B_1cpM.*pinf.r1(:,:,4) + pinf.B_2csM.*(csP_T_csM.*(sn1-ndsn_T_nx) + muP_T_csM.*(dvx-ndv_T_nx));
            flux(:,:,5) = B_1cpM.*pinf.r1(:,:,5) + pinf.B_2csM.*(csP_T_csM.*(sn2-ndsn_T_ny) + muP_T_csM.*(dvy-ndv_T_ny));
        end

        %% %%%%%%%%%%%%%%%%%%%%%
        %%% Derivative Terms %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate local derivatives of fields
        ds11dr = pinf.Dr*s11N; ds11ds = pinf.Ds*s11N;
        ds12dr = pinf.Dr*s12N; ds12ds = pinf.Ds*s12N;
        ds22dr = pinf.Dr*s22N; ds22ds = pinf.Ds*s22N;
        dvxdr   = pinf.Dr*vxN;   dvxds   = pinf.Ds*vxN;
        dvydr   = pinf.Dr*vyN;   dvyds   = pinf.Ds*vyN;
        
        % Compute physical derivatives of fields
        ds11dy = pinf.ry.*ds11dr + pinf.sy.*ds11ds;
        ds11dx = pinf.rx.*ds11dr + pinf.sx.*ds11ds;
        ds12dy = pinf.ry.*ds12dr + pinf.sy.*ds12ds;
        ds12dx = pinf.rx.*ds12dr + pinf.sx.*ds12ds;
        ds22dy = pinf.ry.*ds22dr + pinf.sy.*ds22ds;
        ds22dx = pinf.rx.*ds22dr + pinf.sx.*ds22ds;
        dvxdy   = pinf.ry.*dvxdr   + pinf.sy.*dvxds;
        dvxdx   = pinf.rx.*dvxdr   + pinf.sx.*dvxds;
        dvydy   = pinf.ry.*dvydr   + pinf.sy.*dvyds;
        dvydx   = pinf.rx.*dvydr   + pinf.sx.*dvyds;
      
        %% %%%%%%%%%%%%%%%%%%
        %%% Compute qvelo %%%
        %%%%%%%%%%%%%%%%%%%%%
        %Defining matrix A and B entries
        A14_D_rho = pinf.A14_D_rho;
        A24_D_rho = pinf.A24_D_rho;
        B15_D_rho = A24_D_rho;
        B25_D_rho = A14_D_rho;
        A35_D_rho = pinf.mu_D_rho;
        B34_D_rho = pinf.mu_D_rho;
        
        rhs_s11(ids) = A14_D_rho.*dvxdx + B15_D_rho.*dvydy - pinf.LIFT*(pinf.Fscale.*flux(:,:,1));
        rhs_s22(ids) = A24_D_rho.*dvxdx + B25_D_rho.*dvydy - pinf.LIFT*(pinf.Fscale.*flux(:,:,2));
        rhs_s12(ids) = A35_D_rho.*dvxdy + B34_D_rho.*dvydx - pinf.LIFT*(pinf.Fscale.*flux(:,:,3));
        rhs_vx(ids)  = ds11dx + ds12dy - pinf.LIFT*(pinf.Fscale.*flux(:,:,4));
        rhs_vy(ids)  = ds12dx + ds22dy - pinf.LIFT*(pinf.Fscale.*flux(:,:,5));
    end
end             