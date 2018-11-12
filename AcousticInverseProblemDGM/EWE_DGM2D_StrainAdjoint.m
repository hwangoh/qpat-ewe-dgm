function [rhs_e11Adj,rhs_e22Adj,rhs_e12Adj,rhs_vxAdj,rhs_vyAdj]=EWE_DGM2D_StrainAdjoint(RunOptions,NpI,KI,rho,mu,lambda,e11Adj,e22Adj,e12Adj,vxAdj,vyAdj,e11Data,e22Data,e12Data,vxData,vyData,e11Recon,e22Recon,e12Recon,vxRecon,vyRecon,c)

% EWE_DGM2D_StrainAdjoint.m is the main script file for the discontinuous Galerkin method
% implementation of the adjoint of the acoustic forward problem using the elastic wave equation in the strain-velocity form. 
% Utilizes codes from the textbook "Nodal Discontinous Galerkin Methods, Algorithms, Analysis and
% Applications" by Jan Hesthaven and Tim Warburton, 2007
%
% Inputs: 
%   RunOptions:
%          BC - Set boundary conditions
%   NpI - Number of nodes per element in inverse DGM grid
%   KI - Number of elements in inverse DGM grid
%   rho_e - Mass density per unit volume
%   mu_e - First Lame parameter, also known as the shear modulus. Setting the shear modulus to zero yields the acoustic fluid wave equation
%   lambda_e = Second Lame parameter
%   pinfo: structural data with respect to the differing element orders
%   e11Adj: Adjoint variable of strain Tensor entry
%   e12Adj: Adjoint variable of strain Tensor entry
%   e22Adj: Adjoint variable of strain Tensor entry
%   vxAdj: Adjoint variable of x component of the velocity of the displacement
%   vyAdj: Adjoint variable of y component of the velocity of the displacement
%   e11Data - Data for e11 at current time step
%   e22Data - Data for e22 at current time step
%   e12Data - Data for e12 at current time step
%   vxData -  Data for vx at current time step
%   vyData -  Data for vy at current time step
%   e11Recon: Current reconstruction of strain Tensor entry at current time step
%   e12Recon: Current reconstruction of strain Tensor entry at current time step
%   e22Recon: Current reconstruction of strain Tensor entry at current time step
%   vxRecon: Current reconstruction of x component of the velocity of the displacement at current time step
%   vyRecon: Current reconstruction of y component of the velocity of the displacement at current time step
%   c: error model constant
%          
%  Outputs:
%   rhs_e11Adj: Velocity of strain tensor component e11
%   rhs_e22Adj: Velocity of strain tensor component e22
%   rhs_e12ADj: Velocity of strain tensor component e12
%   rhs_vxAdj: Velocity of the x component of the velocity of the displacement
%   rhs_vyAdj: Velocity of the y component of the velocity of the displacement
%
% Hwan Goh, 12/09/2017, University of Auckland, New Zealand

Globals2D;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Structures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Parameters ===%
if RunOptions.UseAcousticViaElasticDGMFormulation == 1;
    mu = 0*mu;
end
c_p = sqrt((lambda + 2*mu)./rho);
c_s = sqrt(mu./rho);

%=== Reshaping Inputs ===%
e11Adj = reshape(e11Adj,NpI,KI);
e22Adj = reshape(e22Adj,NpI,KI);
e12Adj = reshape(e12Adj,NpI,KI);
vxAdj = reshape(vxAdj,NpI,KI);
vyAdj = reshape(vyAdj,NpI,KI);

%=== rhs Storage ===% %Needs to be revised! K defined before K=pinf.K is only fine for a uniform mesh
rhs_e11Adj = zeros(NpI*KI,1);
rhs_e22Adj = zeros(NpI*KI,1);
rhs_e12Adj = zeros(NpI*KI,1);
rhs_vxAdj  = zeros(NpI*KI,1);
rhs_vyAdj  = zeros(NpI*KI,1);

%=== Interior and Exterior Terms ===%
muM = mu(vmapM);
lambdaM = lambda(vmapM);
e11AdjM = e11Adj(vmapM);
e22AdjM = e22Adj(vmapM);
e12AdjM = e12Adj(vmapM);
vxAdjM = vxAdj(vmapM);
vyAdjM = vyAdj(vmapM);
c_pM = c_p(vmapM);
c_sM = c_s(vmapM);

muP = mu(vmapP);
lambdaP = lambda(vmapP);
e11AdjP = e11Adj(vmapP);
e22AdjP = e22Adj(vmapP);
e12AdjP = e12Adj(vmapP);
vxAdjP = vxAdj(vmapP);
vyAdjP = vyAdj(vmapP);
c_pP = c_p(vmapP);
c_sP = c_s(vmapP);

e11AdjP = e11AdjP(:);
e22AdjP = e22AdjP(:);
e12AdjP = e12AdjP(:);
vxAdjP = vxAdjP(:);
vyAdjP = vyAdjP(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Boundary Conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bcvece11Adj = ones(3*Nfp*KI,1);
bcvece22Adj = ones(3*Nfp*KI,1);
bcvece12Adj = ones(3*Nfp*KI,1);
bcvecvxAdj = ones(3*Nfp*KI,1);
bcvecvyAdj = ones(3*Nfp*KI,1);

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

%%%%%%%%%%%%%%%%%%%%%%
%%% Numerical Flux %%%
%%%%%%%%%%%%%%%%%%%%%%
%=== [[S]] ===%
sn1 = 2*muM.*(e11AdjM.*nx(:) + e12AdjM.*ny(:)) + lambdaM.*(e11AdjM+e22AdjM).*nx(:) +...
         (2*muP.*(bcvece11Adj.*e11AdjP.*-nx(:)+bcvece12Adj.*e12AdjP.*-ny(:)) + lambdaP.*(bcvece11Adj.*e11AdjP+bcvece22Adj.*e22AdjP).*-nx(:)); %First row of [[S]]
sn2 = 2*muM.*(e12AdjM.*nx(:) + e22AdjM.*ny(:)) + lambdaM.*(e11AdjM+e22AdjM).*ny(:) +...
         (2*muP.*(bcvece12Adj.*e12AdjP.*-nx(:)+bcvece22Adj.*e22AdjP.*-ny(:)) + lambdaP.*(bcvece11Adj.*e11AdjP+bcvece22Adj.*e22AdjP).*-ny(:)); %Second row of [[S]]

%=== [[v_e]] ===%
ndv = nx(:).*vxAdjM + ny(:).*vyAdjM + -nx(:).*bcvecvxAdj.*vxAdjP + -ny(:).*bcvecvyAdj.*vyAdjP;
%=== [v_e] ===%
dvx = vxAdjM - bcvecvxAdj.*vxAdjP;
dvy = vyAdjM - bcvecvyAdj.*vyAdjP;

%=== n'[[S]] ===% 
ndsn = 2*muM.*(nx(:).*nx(:).*e11AdjM + 2*nx(:).*ny(:).*e12AdjM + ny(:).*ny(:).*e22AdjM) + lambdaM.*(e11AdjM+e22AdjM) +...
         (2*muP.*(-nx(:).*nx(:).*bcvece11Adj.*e11AdjP + -2*nx(:).*ny(:).*bcvece12Adj.*e12AdjP + -ny(:).*ny(:).*bcvece22Adj.*e22AdjP) + -lambdaP.*(bcvece11Adj.*e11AdjP+bcvece22Adj.*e22AdjP));

%=== B_1c_pMr_1M ===%
B_1c_pMr_1M = sparse(5*3*Nfp,KI);
B_1c_pM = (c_pM.*c_pP.*ndsn - c_pM.*(lambdaP+2*muP).*ndv)./...
                 (c_pP.*(lambdaM+2*muM) + c_pM.*(lambdaP+2*muP));
B_1c_pMr_1M(1:3*Nfp*KI) = -B_1c_pM(:).*nx(:).*nx(:);
B_1c_pMr_1M(1*3*Nfp*KI+1:2*3*Nfp*KI) = -B_1c_pM(:).*ny(:).*ny(:);  
B_1c_pMr_1M(2*3*Nfp*KI+1:3*3*Nfp*KI) = -B_1c_pM(:).*nx(:).*ny(:);
B_1c_pMr_1M(3*3*Nfp*KI+1:4*3*Nfp*KI) = B_1c_pM(:).*c_pM.*nx(:);
B_1c_pMr_1M(4*3*Nfp*KI+1:5*3*Nfp*KI) = B_1c_pM(:).*c_pM.*ny(:);

%=== B_2c_sMr_2M ===%
B_2c_sMr_2M = sparse(5*3*Nfp,KI);
if RunOptions.UseAcousticViaElasticDGMFormulation ~= 1;
ndsnn_x = ndsn.*nx(:);
ndsnn_y = ndsn.*ny(:);
ndvn_x = ndv.*nx(:);
ndvn_y = ndv.*ny(:);
B_2c_sM = 1./(c_sP.*muM+c_sM.*muP);
B_2c_sMr_2M(1:3*Nfp*KI) =  B_2c_sM(:).*(-(c_sM.*c_sP).*(nx(:).*(sn1 - ndsnn_x)) +...
                                       (c_sM.*muP).*(nx(:).*(dvx - ndvn_x)));
B_2c_sMr_2M(1*3*Nfp*KI+1:2*3*Nfp*KI) = B_2c_sM(:).*(-(c_sM.*c_sP).*(ny(:).*(sn2 - ndsnn_y)) +...
                                                  (c_sM.*muP).*(ny(:).*(dvy - ndvn_y)));
B_2c_sMr_2M(2*3*Nfp*KI+1:3*3*Nfp*KI) = B_2c_sM(:).*(-(c_sM.*c_sP).*(nx(:)/2.*(sn2 - ndsnn_y) + ny(:)/2.*(sn1 - ndsnn_x))+...
                                                  (c_sM.*muP).*(nx(:)/2.*(dvy - ndvn_y) + ny(:)/2.*(dvx - ndvn_x)));
B_2c_sMr_2M(3*3*Nfp*KI+1:4*3*Nfp*KI) = B_2c_sM(:).*((c_sM.*c_sP).*(-c_sM.*(sn1 - ndsnn_x) +...
                                                  -(c_sM.*muP).*(c_sM.*(dvx - ndvn_x))));
B_2c_sMr_2M(4*3*Nfp*KI+1:5*3*Nfp*KI) = B_2c_sM(:).*((c_sM.*c_sP).*(-c_sM.*-(sn2 - ndsnn_y) +...
                                                  -(c_sM.*muP).*(c_sM.*(dvy - ndvn_y))));
end


%=== Flux ===%
flux = (B_1c_pMr_1M(:) + B_2c_sMr_2M(:));
flux = [reshape(flux(1:3*Nfp*KI,1),3*Nfp,KI);reshape(flux(1*3*Nfp*KI+1:2*3*Nfp*KI,1),3*Nfp,KI);reshape(flux(2*3*Nfp*KI+1:3*3*Nfp*KI,1),3*Nfp,KI);reshape(flux(3*3*Nfp*KI+1:4*3*Nfp*KI,1),3*Nfp,KI);reshape(flux(4*3*Nfp*KI+1:5*3*Nfp*KI,1),3*Nfp,KI)];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating q_velo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Error Model Constant
c1 = c(1);
c2 = c(NpI*KI+1);
c3 = c(2*NpI*KI+1);
c4 = c(3*NpI*KI+1);
c5 = c(4*NpI*KI+1);

if RunOptions.UseSensorData == 1
    rhs_e11Adj = (rx.*(Dr*vxAdj)+sx.*(Ds*vxAdj)) - LIFT*(Fscale.*flux(1:3*Nfp,1:KI)) - c1*(V*V')*(e11Data - e11Recon);
    rhs_e22Adj = (ry.*(Dr*vyAdj)+sy.*(Ds*vyAdj)) - LIFT*(Fscale.*flux(1*3*Nfp+1:2*3*Nfp,1:KI)) - c2*(V*V')*(e22Data - e22Recon);
    rhs_e12Adj = (1/2)*(rx.*(Dr*vyAdj)+sx.*(Ds*vyAdj)) + (1/2)*(ry.*(Dr*vxAdj)+sy.*(Ds*vxAdj)) - LIFT*(Fscale.*flux(2*3*Nfp+1:3*3*Nfp,1:KI)) - c3*(V*V')*(e12Data - e12Recon);
    rhs_vxAdj = (1./rho).*((2*mu+lambda).*(rx.*(Dr*e11Adj)+sx.*(Ds*e11Adj)) + lambda.*(rx.*(Dr*e22Adj)+sx.*(Ds*e22Adj)) + 2*mu.*(ry.*(Dr*e12Adj)+sy.*(Ds*e12Adj))) - LIFT*(Fscale.*flux(3*3*Nfp+1:4*3*Nfp,1:KI)) - c4*(V*V')*(vxData - vxRecon);
    rhs_vyAdj = (1./rho).*(2*mu.*(rx.*(Dr*e12Adj)+sx.*(Ds*e12Adj)) + lambda.*(ry.*(Dr*e11Adj) + sy.*(Ds*e11Adj)) + (2*mu+lambda).*(ry.*(Dr*e22Adj)+sy.*(Ds*e22Adj))) - LIFT*(Fscale.*flux(4*3*Nfp+1:5*3*Nfp,1:KI)) - c5*(V*V')*(vyData - vyRecon);
else
    rhs_e11Adj = (rx.*(Dr*vxAdj)+sx.*(Ds*vxAdj)) - LIFT*(Fscale.*flux(1:3*Nfp,1:KI)) - c1*(e11Data - e11Recon);
    rhs_e22Adj = (ry.*(Dr*vyAdj)+sy.*(Ds*vyAdj)) - LIFT*(Fscale.*flux(1*3*Nfp+1:2*3*Nfp,1:KI)) - c2*(e22Data - e22Recon);
    rhs_e12Adj = (1/2)*(rx.*(Dr*vyAdj)+sx.*(Ds*vyAdj)) + (1/2)*(ry.*(Dr*vxAdj)+sy.*(Ds*vxAdj)) - LIFT*(Fscale.*flux(2*3*Nfp+1:3*3*Nfp,1:KI)) - c3*(e12Data - e12Recon);
    rhs_vxAdj = (1./rho).*((2*mu+lambda).*(rx.*(Dr*e11Adj)+sx.*(Ds*e11Adj)) + lambda.*(rx.*(Dr*e22Adj)+sx.*(Ds*e22Adj)) + 2*mu.*(ry.*(Dr*e12Adj)+sy.*(Ds*e12Adj))) - LIFT*(Fscale.*flux(3*3*Nfp+1:4*3*Nfp,1:KI)) - c4*(vxData - vxRecon);
    rhs_vyAdj = (1./rho).*(2*mu.*(rx.*(Dr*e12Adj)+sx.*(Ds*e12Adj)) + lambda.*(ry.*(Dr*e11Adj) + sy.*(Ds*e11Adj)) + (2*mu+lambda).*(ry.*(Dr*e22Adj)+sy.*(Ds*e22Adj))) - LIFT*(Fscale.*flux(4*3*Nfp+1:5*3*Nfp,1:KI)) - c5*(vyData - vyRecon);
end

rhs_e11Adj = rhs_e11Adj(:);
rhs_e22Adj = rhs_e22Adj(:);
rhs_e12Adj = rhs_e12Adj(:);
rhs_vxAdj = rhs_vxAdj(:);
rhs_vyAdj = rhs_vyAdj(:);

return;
