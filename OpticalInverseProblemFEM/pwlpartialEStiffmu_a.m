function partialEStiff = pwlpartialEStiffmu_a(Coords,kappa)

% pwlpartialEStiffmu_a computes the partial derivative of the stiffness matrix 
% with respect to mu_(a_i) for the nodes corresponding to the jth element,
% assuming a piecewise linear approximation of the parameter
%
% Inputs:
%    Coords - a 3 by 2 matrix where each row contains coordinates of the
%        vertices for the jth element
%    kappa - value of kappa at the ith node
%    
% Outputs:
%    partialEStiff - the 3 by 3 matrix which represents the entries of the
%                 matrix delStiff/delmu_(a_i)
%
% Hwan Goh 29/08/2013, University of Auckland, New Zealand

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Jacobian of Global to Reference element mapping %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = [-1,1,0;-1,0,1];
Jac_T = L*Coords; % The transpose of the Jacobian
S = Jac_T\L; % The Inverse of the Jacobian multiplied by L
detJac = abs(det(Jac_T)); % Determinant of the Jacobian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Partial Derivative of Stiffness Matrix %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partialEStiff = 1/6*S'*S*detJac*(-2*kappa^2);



