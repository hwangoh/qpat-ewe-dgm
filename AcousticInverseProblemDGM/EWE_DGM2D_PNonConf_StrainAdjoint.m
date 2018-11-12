function [rhs_e11Adj,rhs_e22Adj,rhs_e12Adj,rhs_vxAdj,rhs_vyAdj]=EWE_DGM2D_PNonConf_StrainAdjoint(RunOptions,pinfoI,e11Adj,e22Adj,e12Adj,vxAdj,vyAdj,e11DerivObjFunct,e22DerivObjFunct,e12DerivObjFunct,vxDerivObjFunct,vyDerivObjFunct,NpI,KI)

% EWE_DGM2D_PNonConf_StrainAdjoint.m is the main script file for the discontinuous Galerkin method
% implementation of the adjoint of the acoustic forward problem using the elastic wave
% equation in strain-velocity form. This code is enabled for p-nonconforming meshes and is written in a
% style more similar to Timo Lahivaara's code and also
% utilizes codes from the textbook "Nodal Discontinous Galerkin Methods, Algorithms, Analysis and
% Applications" by Jan Hesthaven and Tim Warburton, 2007
%
% Inputs: 
%   RunOptions:
%          BC - Set boundary conditions
%   pinfoI: structural data with respect to the differing element orders
%   e11Adj: Adjoint variable of strain Tensor entry
%   e12Adj: Adjoint variable of strain Tensor entry
%   e22Adj: Adjoint variable of strain Tensor entry
%   vxAdj: Adjoint variable of x component of the velocity of the displacement
%   vyAdj: Adjoint variable of y component of the velocity of the displacement
%   e11ObjFunct - e11 component of derivative of objective functional
%   e22ObjFunct - e22 component of derivative of objective functional
%   e12ObjFunct - e12 component of derivative of objective functional
%   vxObjFunct -  vx component of derivative of objective functional
%   vyObjFunct -  vy component of derivative of objective functional
%   NpI - Number of nodes per element in inverse DGM grid
%   KI - Number of elements in inverse DGM grid
%          
%  Outputs:
%   rhs_e11Adj: Velocity of strain tensor component e11
%   rhs_e22Adj: Velocity of strain tensor component e22
%   rhs_e12ADj: Velocity of strain tensor component e12
%   rhs_vxAdj: Velocity of the x component of the velocity of the displacement
%   rhs_vyAdj: Velocity of the y component of the velocity of the displacement
%
% Hwan Goh, 12/09/2017, University of Auckland, New Zealand
% Last Edited: 16/11/2017 - removed reliance on Globals2D and save everything into structures
% Note: NEED TO FIX pinfo FOR INVERSE MESH!

Nfaces = 3;

%% %%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Structures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== rhs Storage ===% %Needs to be revised! K defined before K=pinf.K is only fine for a uniform mesh
rhs_e11Adj = zeros(NpI*KI,1);
rhs_e22Adj = zeros(NpI*KI,1);
rhs_e12Adj = zeros(NpI*KI,1);
rhs_vxAdj   = zeros(NpI*KI,1);
rhs_vyAdj  = zeros(NpI*KI,1);

for N1=1:length(pinfoI)
    pinf = pinfoI(N1);
    NpI = pinf.Np;
    KI = pinf.K;
    Nfp = pinf.Nfp;
    mapO = pinf.mapO;
    mapD = pinf.mapD;
    mapN = pinf.mapN;
    V = pinf.V;
    if (KI>0)
        ids = pinf.ids;
        Fmask = pinf.Fmask;
        e11AdjN = e11Adj(ids);
        e22AdjN = e22Adj(ids);
        e12AdjN = e12Adj(ids);
        vxAdjN = vxAdj(ids);
        vyAdjN = vyAdj(ids);
        %=== Interior and Exterior Terms ===%
        e11AdjM = e11AdjN(Fmask(:),:); 
        e22AdjM = e22AdjN(Fmask(:),:);
        e12AdjM = e12AdjN(Fmask(:),:);
        vxAdjM   = vxAdjN(Fmask(:),:);
        vyAdjM   = vyAdjN(Fmask(:),:);
        e11AdjP = zeros(Nfaces*Nfp,KI);
        e22AdjP = zeros(Nfaces*Nfp,KI);
        e12AdjP = zeros(Nfaces*Nfp,KI);
        vxAdjP = zeros(Nfaces*Nfp,KI);
        vyAdjP = zeros(Nfaces*Nfp,KI);
        for N2 = 1:length(pinfoI)
            if ~isempty(pinf.fmapM{N2})
                interpMAT = pinf.interpP{N2};
                fmapM = pinf.fmapM{N2};
                vmapP = pinf.vmapP{N2};
                e11AdjP(fmapM) = interpMAT*e11Adj(vmapP);
                e12AdjP(fmapM) = interpMAT*e12Adj(vmapP);
                e22AdjP(fmapM) = interpMAT*e22Adj(vmapP);
                vxAdjP(fmapM)   = interpMAT*vxAdj(vmapP);
                vyAdjP(fmapM)   = interpMAT*vyAdj(vmapP);
            end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Set Boundary Conditions %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bcvece11Adj = ones(Nfaces*Nfp,KI);
        bcvece22Adj = ones(Nfaces*Nfp,KI);
        bcvece12Adj = ones(Nfaces*Nfp,KI);
        bcvecvxAdj = ones(Nfaces*Nfp,KI);
        bcvecvyAdj = ones(Nfaces*Nfp,KI);
        
        %Out/Absorbing
        bcvece11Adj(mapO) = 0;
        bcvece22Adj(mapO) = 0;
        bcvece12Adj(mapO) = 0;
        bcvecvxAdj(mapO) = 0;
        bcvecvyAdj(mapO) = 0;
        %Dirichlet
        bcvece11Adj(mapD) = 1;
        bcvece22Adj(mapD) = 1;
        bcvece12Adj(mapD) = 1;
        bcvecvxAdj(mapD) = -1;
        bcvecvyAdj(mapD) = -1;
        %Neumann
        bcvece11Adj(mapN) = -1;
        bcvece22Adj(mapN) = -1;
        bcvece12Adj(mapN) = 1;
        bcvecvxAdj(mapN) = 1;
        bcvecvyAdj(mapN) = 1;
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Some Precomputed Values %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % _T_ means .* and _P_ means + and _D_ means ./
        %=== Sums of Strain Values ===%
        e11M_P_e22M = e11AdjM+e22AdjM;
        bcvece11_T_e11P_P_bcvece22_T_e22P = bcvece11Adj.*e11AdjP+bcvece22Adj.*e22AdjP;
        
        %=== Sums of Wave Speeds and Medium Densities ===%
        rhoP_T_csM = pinf.rhoP.*pinf.csM;
        %=== [[v_e]] ===%
        ndv  = pinf.nx.*vxAdjM+pinf.ny.*vyAdjM-(pinf.nx.*bcvecvxAdj.*vxAdjP+pinf.ny.*bcvecvyAdj.*vyAdjP);
        %=== [v_e] ===%
        dvx = vxAdjM-bcvecvxAdj.*vxAdjP;
        dvy = vyAdjM-bcvecvyAdj.*vyAdjP;
        %=== [[S]] ===%
        sn1 = 2*pinf.muM.*(e11AdjM.*pinf.nx+pinf.ny.*e12AdjM)+pinf.lambdaM.*e11M_P_e22M.*pinf.nx-...
             (2*pinf.muP.*(bcvece11Adj.*e11AdjP.*pinf.nx+pinf.ny.*bcvece12Adj.*e12AdjP)+pinf.lambdaP.*bcvece11_T_e11P_P_bcvece22_T_e22P.*pinf.nx); %First row of [[S]]
        sn2 = 2*pinf.muM.*(e12AdjM.*pinf.nx+pinf.ny.*e22AdjM)+ pinf.lambdaM.*e11M_P_e22M.*pinf.ny-...
             (2*pinf.muP.*(bcvece12Adj.*e12AdjP.*pinf.nx+pinf.ny.*bcvece22Adj.*e22AdjP)+pinf.lambdaP.*bcvece11_T_e11P_P_bcvece22_T_e22P.*pinf.ny); %Second row of [[S]]
        %=== n'[[S]] ===%
        ndsn = 2*pinf.muM.*(pinf.nxnx.*e11AdjM + 2*pinf.nxny.*e12AdjM + pinf.nyny.*e22AdjM)+pinf.lambdaM.*e11M_P_e22M-...
              (2*pinf.muP.*(pinf.nxnx.*bcvece11Adj.*e11AdjP + 2*pinf.nxny.*bcvece12Adj.*e12AdjP + pinf.nyny.*bcvece22Adj.*e22AdjP) + pinf.lambdaP.*bcvece11_T_e11P_P_bcvece22_T_e22P);
        %=== Some other values for the flux ===%
        ndsn_T_nx = ndsn.*pinf.nx;
        ndv_T_nx = ndv.*pinf.nx;
        ndsn_T_ny = ndsn.*pinf.ny;
        ndv_T_ny = ndv.*pinf.ny;
        %=== B1cpM ===%
        B_1cpM = pinf.cpM.*(pinf.d11.*ndsn - pinf.d13.*ndv);
        
        %% %%%%%%%%%
        %%% Flux %%%
        %%%%%%%%%%%%
        flux = zeros(Nfaces*Nfp, KI, 5);        
        flux(:,:,1) = -B_1cpM.*pinf.r1(:,:,1) + pinf.B_2csM.*(-pinf.nx.*(sn1-ndsn_T_nx)                                         + rhoP_T_csM.*pinf.nx.*(dvx-ndv_T_nx));
        flux(:,:,2) = -B_1cpM.*pinf.r1(:,:,2) + pinf.B_2csM.*(-pinf.ny.*(sn2-ndsn_T_ny)                                         + rhoP_T_csM.*pinf.ny.*(dvy-ndv_T_ny));
        flux(:,:,3) = -B_1cpM.*pinf.r1(:,:,3) + pinf.B_2csM.*(-(pinf.nx/2.*(sn2-ndsn_T_ny) + pinf.ny/2.*(sn1-ndsn_T_nx))        + rhoP_T_csM.*(pinf.nx/2.*(dvy-ndv_T_ny) + pinf.ny/2.*(dvx-ndv_T_nx)));
        flux(:,:,4) =  B_1cpM.*pinf.r1(:,:,4) + pinf.B_2csM.*(pinf.csM.*(sn1-ndsn_T_nx)                                         - rhoP_T_csM.*pinf.csM.*(dvx-ndv_T_nx));
        flux(:,:,5) =  B_1cpM.*pinf.r1(:,:,5) + pinf.B_2csM.*(pinf.csM.*(sn2-ndsn_T_ny)                                         - rhoP_T_csM.*pinf.csM.*(dvy-ndv_T_ny));
        
        %% %%%%%%%%%%%%%%%%%%%%%
        %%% Derivative Terms %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate local derivatives of fields
        de11dr = pinf.Dr*e11AdjN; de11ds = pinf.Ds*e11AdjN;
        de12dr = pinf.Dr*e12AdjN; de12ds = pinf.Ds*e12AdjN;
        de22dr = pinf.Dr*e22AdjN; de22ds = pinf.Ds*e22AdjN;
        dvxdr  = pinf.Dr*vxAdjN; dvxds  = pinf.Ds*vxAdjN;
        dvydr  = pinf.Dr*vyAdjN; dvyds  = pinf.Ds*vyAdjN;
        
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
        
        %Inverse of Mass Matrix
        Minv = (V*V');
        
        %Computing RHS,(V*V') is the inverse of the mass matrix
        rhs_e11Adj(ids) = dvxdx - pinf.LIFT*(pinf.Fscale.*flux(:,:,1)) - Minv*e11DerivObjFunct;
        rhs_e22Adj(ids) = dvydy - pinf.LIFT*(pinf.Fscale.*flux(:,:,2)) - Minv*e22DerivObjFunct;
        rhs_e12Adj(ids) = (1/2)*dvxdy + (1/2)*dvydx - pinf.LIFT*(pinf.Fscale.*flux(:,:,3)) - Minv*e12DerivObjFunct;
        rhs_vxAdj(ids)  = A41_D_rho.*de11dx + A42_D_rho.*de22dx + B43_D_rho.*de12dy - pinf.LIFT*(pinf.Fscale.*flux(:,:,4)) - Minv*vxDerivObjFunct;
        rhs_vyAdj(ids)  = B51_D_rho.*de11dy + B52_D_rho.*de22dy + A53_D_rho.*de12dx - pinf.LIFT*(pinf.Fscale.*flux(:,:,5)) - Minv*vyDerivObjFunct;
    end
end