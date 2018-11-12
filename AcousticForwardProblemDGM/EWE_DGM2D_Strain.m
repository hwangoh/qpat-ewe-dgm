function [rhs_e11,rhs_e22,rhs_e12,rhs_vx,rhs_vy]=EWE_DGM2D_Strain(RunOptions,DataVrbls,e11,e22,e12,vx,vy)

% EWE_DGM2D_Strain.m is the main script file for the discontinuous Galerkin method
% implementation of the acoustic forward problem using the elastic wave equation in the strain-velocity form. 
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
%   e_11: Strain Tensor entry
%   e_12: Strain Tensor entry
%   e_22: Strain Tensor entry
%   vx: x component of the velocity of the displacement
%   vy: y component of the velocity of the displacement
%          
%  Outputs:
%   rhs_e11: Velocity of strain tensor component e_11
%   rhs_e22: Velocity of strain tensor component e_22
%   rhs_e12: Velocity of strain tensor component e_12
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
% Hwan Goh, 15/06/2016 (sunny day in Kuopio!), University of Auckland, New Zealand

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
e11 = reshape(e11,Np,K);
e22 = reshape(e22,Np,K);
e12 = reshape(e12,Np,K);
vx = reshape(vx,Np,K);
vy = reshape(vy,Np,K);

%=== rhs Storage ===% %Needs to be revised! K defined before K=pinf.K is only fine for a uniform mesh
rhs_e11 = zeros(Np*K,1);
rhs_e22 = zeros(Np*K,1);
rhs_e12 = zeros(Np*K,1);
rhs_vx   = zeros(Np*K,1);
rhs_vy  = zeros(Np*K,1);

%=== Interior and Exterior Terms ===%
muM = mu(vmapM);
lambdaM = lambda(vmapM);
e11M = e11(vmapM);
e22M = e22(vmapM);
e12M = e12(vmapM);
vxM = vx(vmapM);
vyM = vy(vmapM);
c_pM = c_p(vmapM);
c_sM = c_s(vmapM);

muP = mu(vmapP);
lambdaP = lambda(vmapP);
e11P = e11(vmapP);
e22P = e22(vmapP);
e12P = e12(vmapP);
vxP = vx(vmapP);
vyP = vy(vmapP);
c_pP = c_p(vmapP);
c_sP = c_s(vmapP);

e11P = e11P(:);
e22P = e22P(:);
e12P = e12P(:);
vxP = vxP(:);
vyP = vyP(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Boundary Conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bcvece11 = ones(3*Nfp*K,1);
bcvece22 = ones(3*Nfp*K,1);
bcvece12 = ones(3*Nfp*K,1);
bcvecvx = ones(3*Nfp*K,1);
bcvecvy = ones(3*Nfp*K,1);

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

%%%%%%%%%%%%%%%%%%%%%%
%%% Numerical Flux %%%
%%%%%%%%%%%%%%%%%%%%%%
%=== [[S]] ===%
sn1 = 2*muM.*(e11M.*nx(:) + e12M.*ny(:)) + lambdaM.*(e11M+e22M).*nx(:) +...
         (2*muP.*(bcvece11.*e11P.*-nx(:)+bcvece12.*e12P.*-ny(:)) + lambdaP.*(bcvece11.*e11P+bcvece22.*e22P).*-nx(:)); %First row of [[S]]
sn2 = 2*muM.*(e12M.*nx(:) + e22M.*ny(:)) + lambdaM.*(e11M+e22M).*ny(:) +...
         (2*muP.*(bcvece12.*e12P.*-nx(:)+bcvece22.*e22P.*-ny(:)) + lambdaP.*(bcvece11.*e11P+bcvece22.*e22P).*-ny(:)); %Second row of [[S]]

%=== [[v_e]] ===%
ndv = nx(:).*vxM + ny(:).*vyM + -nx(:).*bcvecvx.*vxP + -ny(:).*bcvecvy.*vyP;
%=== [v_e] ===%
dvx = vxM - bcvecvx.*vxP;
dvy = vyM - bcvecvy.*vyP;

%=== n'[[S]] ===% 
ndsn = 2*muM.*(nx(:).*nx(:).*e11M + 2*nx(:).*ny(:).*e12M + ny(:).*ny(:).*e22M) + lambdaM.*(e11M+e22M) +...
         (2*muP.*(-nx(:).*nx(:).*bcvece11.*e11P + -2*nx(:).*ny(:).*bcvece12.*e12P + -ny(:).*ny(:).*bcvece22.*e22P) + -lambdaP.*(bcvece11.*e11P+bcvece22.*e22P));

%=== B_1c_pMr_1M ===%
B_1c_pMr_1M = sparse(5*3*Nfp,K);
B_1c_pM = (c_pM.*c_pP.*ndsn + c_pM.*(lambdaP+2*muP).*ndv)./...
                 (c_pP.*(lambdaM+2*muM) + c_pM.*(lambdaP+2*muP));
B_1c_pMr_1M(1:3*Nfp*K) = B_1c_pM(:).*nx(:).*nx(:);
B_1c_pMr_1M(1*3*Nfp*K+1:2*3*Nfp*K) = B_1c_pM(:).*ny(:).*ny(:);  
B_1c_pMr_1M(2*3*Nfp*K+1:3*3*Nfp*K) = B_1c_pM(:).*nx(:).*ny(:);
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
B_2c_sMr_2M(1:3*Nfp*K) =  B_2c_sM(:).*(c_sM.*c_sP).*(nx(:).*(sn1 - ndsnn_x)) +...
                          B_2c_sM(:).*(c_sM.*muP).*(nx(:).*(dvx - ndvn_x));
B_2c_sMr_2M(1*3*Nfp*K+1:2*3*Nfp*K) = B_2c_sM(:).*((c_sM.*c_sP).*(ny(:).*(sn2 - ndsnn_y)) +...
                                                  (c_sM.*muP).*(ny(:).*(dvy - ndvn_y)));
B_2c_sMr_2M(2*3*Nfp*K+1:3*3*Nfp*K) = B_2c_sM(:).*((c_sM.*c_sP).*(nx(:)/2.*(sn2 - ndsnn_y) + ny(:)/2.*(sn1 - ndsnn_x))+...
                                                  (c_sM.*muP).*(nx(:)/2.*(dvy - ndvn_y) + ny(:)/2.*(dvx - ndvn_x)));
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
rhs_e11 = (rx.*(Dr*vx)+sx.*(Ds*vx)) - LIFT*(Fscale.*flux(1:3*Nfp,1:K));
rhs_e22 = (ry.*(Dr*vy)+sy.*(Ds*vy)) - LIFT*(Fscale.*flux(1*3*Nfp+1:2*3*Nfp,1:K));
rhs_e12 = (1/2)*(rx.*(Dr*vy)+sx.*(Ds*vy)) + (1/2)*(ry.*(Dr*vx)+sy.*(Ds*vx)) - LIFT*(Fscale.*flux(2*3*Nfp+1:3*3*Nfp,1:K));
rhs_vx = (1./rho).*((2*mu+lambda).*(rx.*(Dr*e11)+sx.*(Ds*e11)) + lambda.*(rx.*(Dr*e22)+sx.*(Ds*e22)) + 2*mu.*(ry.*(Dr*e12)+sy.*(Ds*e12))) - LIFT*(Fscale.*flux(3*3*Nfp+1:4*3*Nfp,1:K));
rhs_vy = (1./rho).*(2*mu.*(rx.*(Dr*e12)+sx.*(Ds*e12)) + lambda.*(ry.*(Dr*e11) + sy.*(Ds*e11)) + (2*mu+lambda).*(ry.*(Dr*e22)+sy.*(Ds*e22))) - LIFT*(Fscale.*flux(4*3*Nfp+1:5*3*Nfp,1:K));

rhs_e11 = rhs_e11(:);
rhs_e22 = rhs_e22(:);
rhs_e12 = rhs_e12(:);
rhs_vx = rhs_vx(:);
rhs_vy = rhs_vy(:);

return;
