function partialEMass = pwlpartialEMassmu_a(E,pos)

% partialEMassmu_a computes the partial derivative of the mass matrix 
% with respect to mu_(a_i) for the nodes corresponding to the jth element,
% assuming a piecewise linear approximation of the parameter
%
% Inputs:
%    E - a 3 by 2 matrix where each row contains coordinates of the
%        vertices for the jth element
%    pos - finds whether i is the first, second or third listed node of a
%          row of the Elements matrix
%    
% Outputs:
%    partialEMass - the 3 by 3 matrix which represents the entries of the
%                 matrix delMass/delmu_(a_i)
%
% Hwan Goh 29/08/2013, University of Auckland, New Zealand

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Jacobian of Global to Reference element mapping %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = [-1,1,0;-1,0,1];
Jac_T = L*E; % The transpose of the Jacobian
detJac = abs(det(Jac_T)); % Determinant of the Jacobian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Jacobian of Global to Reference element mapping %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(3,3);

if pos == 1
A(1,1) = 1/20;
A(1,2) = 1/60;
A(1,3) = 1/60;

A(2,1) = 1/60;
A(2,2) = 1/60;
A(2,3) = 1/120;

A(3,1) = 1/60;
A(3,2) = 1/120;
A(3,3) = 1/60;
end

if pos == 2
A(1,1) = 1/60;
A(1,2) = 1/60;
A(1,3) = 1/120;

A(2,1) = 1/60;
A(2,2) = 1/20;
A(2,3) = 1/60;

A(3,1) = 1/120;
A(3,2) = 1/60;
A(3,3) = 1/60;
end

if pos == 3
A(1,1) = 1/60;
A(1,2) = 1/120;
A(1,3) = 1/60;

A(2,1) = 1/120;
A(2,2) = 1/60;
A(2,3) = 1/60;

A(3,1) = 1/60;
A(3,2) = 1/60;
A(3,3) = 1/20;
end

partialEMass = detJac*A;
    