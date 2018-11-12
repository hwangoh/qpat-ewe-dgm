function int = pwlIntPosT2ndDerivMass(p0,VerticesCoords,MassCoeffs,detJacT,PrecompCoeffMatricesInv)

% pwlIntPosT2ndDerivMass computes the components of the second derivative
% of the positivity transformed mass matrix elementwise for the piecewise constant case
%
% Inputs:
%    p0 - Input vector of size N_Nodes with respect to FEM Mesh
%    TiiAdjTime0 - Adjoint variables, to be used for computing Newton direction
%    VerticesCoords - a 3 by 2 matrix where each row contains coordinates of the
%                     vertices for the element
%    mu - a 3 by 1 vector where the entries represents the parameter
%         values for that element.
%    detJacT - Determinant of the transpose of the Jacobian
%    PrecompCoeffMatricesInv - N_Elm by 1 cell where each cell contains
%                             the required matrix for computing the coefficients for that element's
%                             interpolating function of two variables. These matrices will be required for
%                             computing the Jacobian of this interpolation function
%
% Outputs:
%    int - the 3 by 3 matrix which represents the entries of the mass matrix  
%          for the element
%
% Hwan Goh 06/07/2018, University of Auckland, New Zealand

A = zeros(3,3);

%=== Quadrature Objects ===%
Weights = [-27/6,25/96,25/96,25/96]; %Gaussian Quadrature Weights
QuadratureCoords = [1/3,1/3;
                    1/5,1/5;
                    3/5,1/5;
                    1/5,3/5];
MappedQuadratureCoords = zeros(4,2);
for ii=1:4
    MappedQuadratureCoords(ii,:) = VerticesCoords(1,:)*(1-QuadratureCoords(ii,1)-QuadratureCoords(ii,2)) + VerticesCoords(2,:)*QuadratureCoords(ii,1) + VerticesCoords(3,:)*QuadratureCoords(ii,2);
end

%=== Function of two variables coefficients ===%
p0InterpCoeffs = PrecompCoeffMatricesInv*p0;
MassCoeffsInterpCoeffs = PrecompCoeffMatricesInv*MassCoeffs;
phiInterpCoeffs = zeros(3,3); %Each column represents coefficients for phi_h
for h=1:3;
    if h==1; %phi_1
        phiInterpCoeffs(:,h)= PrecompCoeffMatricesInv*[1;0;0];
    end
    if h==2; %phi_2
        phiInterpCoeffs(:,h)= PrecompCoeffMatricesInv*[0;1;0];
    end
    if h==3; %phi_3
        phiInterpCoeffs(:,h)= PrecompCoeffMatricesInv*[0;0;1];
    end
end

%=== Computations ===%
for ii=1:3
    for jj=1:3
        for ll=1:4
            A(ii,jj) = A(ii,jj) + Weights(ll)*LogExpT2ndDerivp0phi_jphi_i(p0InterpCoeffs,MassCoeffsInterpCoeffs,MappedQuadratureCoords(ll,:),phiInterpCoeffs,ii,jj);
        end
    end
end
int = detJacT*A;

%% =======================================================================%
%                     Updating Using Current p0Recon
%=========================================================================%
function value = LogExpT2ndDerivp0phi_jphi_i(p0InterpCoeffs,MassCoeffsInterpCoeffs,MappedQuadratureCoords,phiInterpCoeffs,ii,jj)

%=== Interpolated Values ===%
p0InterpValue = p0InterpCoeffs(1) + p0InterpCoeffs(2)*MappedQuadratureCoords(1) + p0InterpCoeffs(3)*MappedQuadratureCoords(2);
MassCoeffsInterpValue = MassCoeffsInterpCoeffs(1) + MassCoeffsInterpCoeffs(2)*MappedQuadratureCoords(1) + MassCoeffsInterpCoeffs(3)*MappedQuadratureCoords(2);
phi_jInterpValue = phiInterpCoeffs(1,jj) + phiInterpCoeffs(2,jj)*MappedQuadratureCoords(1) + phiInterpCoeffs(3,jj)*MappedQuadratureCoords(2);
phi_iInterpValue = phiInterpCoeffs(1,ii) + phiInterpCoeffs(2,ii)*MappedQuadratureCoords(1) + phiInterpCoeffs(3,ii)*MappedQuadratureCoords(2);

%=== Quadrature Value ===%
value = p0InterpValue*MassCoeffsInterpValue*phi_jInterpValue*phi_iInterpValue;
