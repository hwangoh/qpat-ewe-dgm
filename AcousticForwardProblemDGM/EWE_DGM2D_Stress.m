function [rhs_s11,rhs_s22,rhs_s12,rhs_vx,rhs_vy]=EWE_DGM2D_Stress(RunOptions,DataVrbls,s11,s22,s12,vx,vy)

% EWE_DGM2D_Stress.m is the main script file for the discontinuous Galerkin method
% implementation of the acoustic forward problem using the elastic wave equation in stress-velocity form. 
% Utilizes codes from the textbook "Nodal Discontinous Galerkin Methods, Algorithms, Analysis and
% Applications" by Jan Hesthaven and Tim Warburton, 2007
%
% Inputs: 
%   RunOptions:
%          BC - Set boundary conditions
%   DataVrbls:
%          rho_e - Mass density per unit volume
%          mu_e - First Lame parameter, also known as the shear modulus. Setting the shear modulus to zero yields the acoustic fluid wave equation
%          lambda_e = Second Lame parameter
%   s_11: Stress Tensor entry
%   s_12: Stress Tensor entry
%   s_22: Stress Tensor entry
%   vx: x component of the velocity of the displacement
%   vy: y component of the velocity of the displacement
%          
%
%  Outputs:
%   rhs_s11: Velocity of stress tensor component s_11
%   rhs_s22: Velocity of stress tensor component s_22
%   rhs_s12: Velocity of stress tensor component s_12
%   rhs_vx: Velocity of the x component of the velocity of the displacement
%   rhs_vy: Velocity of the y component of the velocity of the displacement
%
% Notes: For the case of polynomial order 2, the 1st and 3rd rows of my vmapP 
%        are switched when I do not use Timo's version of BuildPNonCon2D. 
%        Hence why I have included lines 77 to 79 which switches the rows. 
%        This allows the use of interpMAT in lines 94 to 100 which originates from
%        BuildPNonCon2D. Alternatively, lines 87 and 89 can be uncommented
%        with lines 77 to 79 and lines 94 to 100 commented to directly use Timo's
%        values for the purposes of debugging. Also note that the textbook's
%        version of BuildPNonCon2D creates indexing which, for only some elements,
%        differs from Timo's version of BuildPNonCon2D which I have called
%        BuildPNonCon2DTimo.
%
% Hwan Goh, 27/07/2017 University of Auckland, New Zealand

Globals2D;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Structures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Parameters ===%
rho = DataVrbls.rho_e;
mu = DataVrbls.mu_e;
if RunOptions.UseAcousticViaElasticDGMFormulation == 1;
    mu = 0*mu;
end
lambda = DataVrbls.lambda_e;
c_p = sqrt((lambda + 2*mu)./rho);
c_s = sqrt(mu./rho);

%=== Reshaping Inputs ===%
s11 = reshape(s11,Np,K);
s22 = reshape(s22,Np,K);
s12 = reshape(s12,Np,K);
vx = reshape(vx,Np,K);
vy = reshape(vy,Np,K);

%=== rhs Storage ===% %Needs to be revised! K defined before K=pinf.K is only fine for a uniform mesh
rhs_s11 = zeros(Np*K,1);
rhs_s22 = zeros(Np*K,1);
rhs_s12 = zeros(Np*K,1);
rhs_vx   = zeros(Np*K,1);
rhs_vy  = zeros(Np*K,1);

%=== Interior and Exterior Terms ===%
muM = mu(vmapM);
lambdaM = lambda(vmapM);
s11M = s11(vmapM);
s22M = s22(vmapM);
s12M = s12(vmapM);
vxM = vx(vmapM);
vyM = vy(vmapM);
c_pM = c_p(vmapM);
c_sM = c_s(vmapM);

muP = mu(vmapP);
lambdaP = lambda(vmapP);
s11P = s11(vmapP);
s22P = s22(vmapP);
s12P = s12(vmapP);
vxP = vx(vmapP);
vyP = vy(vmapP);
c_pP = c_p(vmapP);
c_sP = c_s(vmapP);

s11P = s11P(:);
s22P = s22P(:);
s12P = s12P(:);
vxP = vxP(:);
vyP = vyP(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Boundary Conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bcvecs11 = ones(3*Nfp*K,1);
bcvecs22 = ones(3*Nfp*K,1);
bcvecs12 = ones(3*Nfp*K,1);
bcvecvx = ones(3*Nfp*K,1);
bcvecvy = ones(3*Nfp*K,1);

%Out/Absorbing
    bcvece11(mapO) = 0;
    bcvece22(mapO) = 0;
    bcvece12(mapO) = 0;
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

%%%%%%%%%%%%%%%%%%%%%%
%%% Numerical Flux %%%
%%%%%%%%%%%%%%%%%%%%%%
%=== [[S]] ===%
sn1 = s11M.*nx(:) + s12M.*ny(:) + bcvecs11.*s11P.*-nx(:) + bcvecs12.*s12P.*-ny(:);
sn2 = s12M.*nx(:) + s22M.*ny(:) + bcvecs12.*s12P.*-nx(:) + bcvecs22.*s22P.*-ny(:);

%=== [[v_e]] ===%
ndv = nx(:).*vxM + ny(:).*vyM + -nx(:).*bcvecvx.*vxP + -ny(:).*bcvecvy.*vyP;
%=== [v_e] ===%
dvx = vxM - bcvecvx.*vxP;
dvy = vyM - bcvecvy.*vyP;

%=== n'[[S]] ===% 
ndsn = nx(:).*sn1 + ny(:).*sn2;
%=== B_1c_pMr_1M ===%
B_1c_pMr_1M = sparse(5*3*Nfp,K);
B_1c_pM = (c_pM.*c_pP.*ndsn + c_pM.*(lambdaP+2*muP).*ndv)./...
                 (c_pP.*(lambdaM+2*muM) + c_pM.*(lambdaP+2*muP));
B_1c_pMr_1M(1:3*Nfp*K) = B_1c_pM(:).*(lambdaM + 2*muM.*nx(:).*nx(:));
B_1c_pMr_1M(1*3*Nfp*K+1:2*3*Nfp*K) = B_1c_pM(:).*(lambdaM + 2*muM.*ny(:).*ny(:));  
B_1c_pMr_1M(2*3*Nfp*K+1:3*3*Nfp*K) = B_1c_pM(:).*(2*muM).*nx(:).*ny(:);
B_1c_pMr_1M(3*3*Nfp*K+1:4*3*Nfp*K) = B_1c_pM(:).*c_pM.*nx(:);
B_1c_pMr_1M(4*3*Nfp*K+1:5*3*Nfp*K) = B_1c_pM(:).*c_pM.*ny(:);

%=== B_2c_sMr_2M ===%
B_2c_sMr_2M = sparse(5*3*Nfp,K);
if RunOptions.UseAcousticViaElasticDGMFormulation ~= 1;
ndsnn_x = ndsn.*nx(:);
ndsnn_y = ndsn.*ny(:);
ndvn_x = ndv.*nx(:);
ndvn_y = ndv.*ny(:);
B_2c_sM = 1./(c_sP.*muM+c_sM.*muP);
B_2c_sMr_2M(1:3*Nfp*K) =  B_2c_sM(:).*(c_sM.*c_sP).*(2*muM.*nx(:).*(sn1 - ndsnn_x)) +...
                          B_2c_sM(:).*(c_sM.*muP).*(2*muM.*nx(:).*(dvx - ndvn_x));
B_2c_sMr_2M(1*3*Nfp*K+1:2*3*Nfp*K) = B_2c_sM(:).*((c_sM.*c_sP).*(2*muM.*ny(:).*(sn2 - ndsnn_y)) +...
                                                  (c_sM.*muP).*(2*muM.*ny(:).*(dvy - ndvn_y)));
B_2c_sMr_2M(2*3*Nfp*K+1:3*3*Nfp*K) = B_2c_sM(:).*((c_sM.*c_sP).*(2*muM).*(nx(:)/2.*(sn2 - ndsnn_y) + ny(:)/2.*(sn1 - ndsnn_x))+...
                                                  (c_sM.*muP).*(2*muM).*(nx(:)/2.*(dvy - ndvn_y) + ny(:)/2.*(dvx -  ndvn_x)));
B_2c_sMr_2M(3*3*Nfp*K+1:4*3*Nfp*K) = B_2c_sM(:).*((c_sM.*c_sP).*(c_sM.*(sn1 - ndsnn_x) +...
                                                  (c_sM.*muP).*(c_sM.*(dvx - ndvn_x))));
B_2c_sMr_2M(4*3*Nfp*K+1:5*3*Nfp*K) = B_2c_sM(:).*((c_sM.*c_sP).*(c_sM.*(sn2 - ndsnn_y) +...
                                                  (c_sM.*muP).*(c_sM.*(dvy - ndvn_y))));
end


%=== Flux ===%
flux = (B_1c_pMr_1M(:) + B_2c_sMr_2M(:));
flux = [reshape(flux(1:3*Nfp*K,1),3*Nfp,K);reshape(flux(1*3*Nfp*K+1:2*3*Nfp*K,1),3*Nfp,K);reshape(flux(2*3*Nfp*K+1:3*3*Nfp*K,1),3*Nfp,K);reshape(flux(3*3*Nfp*K+1:4*3*Nfp*K,1),3*Nfp,K);reshape(flux(4*3*Nfp*K+1:5*3*Nfp*K,1),3*Nfp,K)];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating q_velo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs_s11 = (lambda + 2*mu).*(rx.*(Dr*vx)+sx.*(Ds*vx)) + lambda.*(ry.*(Dr*vy)+sy.*(Ds*vy)) - LIFT*(Fscale.*flux(1:3*Nfp,1:K));
rhs_s22 = lambda.*(rx.*(Dr*vx)+sx.*(Ds*vx)) + (lambda + 2*mu).*(ry.*(Dr*vy)+sy.*(Ds*vy)) - LIFT*(Fscale.*flux(1*3*Nfp+1:2*3*Nfp,1:K));
rhs_s12 = mu.*(rx.*(Dr*vy)+sx.*(Ds*vy)) + mu.*(ry.*(Dr*vx)+sy.*(Ds*vx)) - LIFT*(Fscale.*flux(2*3*Nfp+1:3*3*Nfp,1:K));
rhs_vx = (1./rho).*((rx.*(Dr*s11)+sx.*(Ds*s11)) + (ry.*(Dr*s12)+sy.*(Ds*s12))) - LIFT*(Fscale.*flux(3*3*Nfp+1:4*3*Nfp,1:K));
rhs_vy = (1./rho).*((rx.*(Dr*s12)+sx.*(Ds*s12)) + (ry.*(Dr*s22)+sy.*(Ds*s22))) - LIFT*(Fscale.*flux(4*3*Nfp+1:5*3*Nfp,1:K));

rhs_s11 = rhs_s11(:);
rhs_s22 = rhs_s22(:);
rhs_s12 = rhs_s12(:);
rhs_vx = rhs_vx(:);
rhs_vy = rhs_vy(:);
return;
