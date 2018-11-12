function int = pwlIntStiffness(VerticesCoords,kappa,L,JacT,detJacT)

% pwlIntStiffness computes the components of the stiffness matrix
% elementwise for the piecewise linear case
%
% Inputs:
%    VerticesCoords - a 3 by 2 matrix where each row contains coordinates of the
%        vertices for the element
%    kappa - a 3 by 1 vector where the entries represents the parameter
%            values for that element.
%    L - 2 by 3 matrix with the gradients of the nodal functions for the
%        reference triangle as its columns.
%    JacT - The transpose of the Jacobian
%    detJacT - Determinant of the transpose of the Jacobian
%    
% Outputs:
%    int - the 3 by 3 matrix which represents the entries of the stiffness matrix 
%          for the element
%
% Hwan Goh 05/05/2013, University of Auckland, New Zealand
% Adapted from P.J. Hadwin 07/07/2010, University of Auckland, New Zealand

%Construction of S
invJacTL = JacT\L; % The inverse of the Jacobian multiplied by L, the "\" uses
              % Gaussian elimination to avoid the need of explicitly
              % constructing the inverse of Jac_T.

%Crudely constructing the matrix S to avoid multiplication and inversion
%Note: I have tried this method and it seems to be slower than the one above
%Jac_T_inv = zeros(2,2);
%Jac_T_inv(1,1) = E(3,2) - E(1,2);
%Jac_T_inv(1,2) = -(E(2,2) - E(1,2));
%Jac_T_inv(2,1) = -(E(3,1) - E(1,1));
%Jac_T_inv(2,2) = E(2,1) - E(1,1);
%S = Jac_T_inv*L;

int = 1/6*invJacTL'*invJacTL*detJacT*(kappa(1)+kappa(2)+kappa(3));

