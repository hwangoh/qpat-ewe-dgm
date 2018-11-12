function int = pwlIntMass(VerticesCoords,mu,detJacT)

% pwcIntMass computes the components of the mass matrix
% elementwise for the piecewise constant case
%
% Inputs:
%    VerticesCoords - a 3 by 2 matrix where each row contains coordinates of the
%        vertices for the element
%    mu - a 3 by 1 vector where the entries represents the parameter
%         values for that element.
%    detJacT - Determinant of the transpose of the Jacobian
%
% Outputs:
%    int - the 3 by 3 matrix which represents the entries of the mass matrix  
%          for the element
%
% Hwan Goh 05/05/2013, University of Auckland, New Zealand
% Adapted from P.J. Hadwin 07/07/2010, University of Auckland, New Zealand

A = zeros(3,3);
A(1,1) = mu(1)/20 + mu(2)/60 + mu(3)/60;
A(1,2) = mu(1)/60 + mu(2)/60 + mu(3)/120;
A(1,3) = mu(1)/60 + mu(2)/120 + mu(3)/60;

A(2,1) = mu(1)/60 + mu(2)/60 + mu(3)/120;
A(2,2) = mu(1)/60 + mu(2)/20 + mu(3)/60;
A(2,3) = mu(1)/120 + mu(2)/60 + mu(3)/60;

A(3,1) = mu(1)/60 + mu(2)/120 + mu(3)/60;
A(3,2) = mu(1)/120 + mu(2)/60 + mu(3)/60;
A(3,3) = mu(1)/60 + mu(2)/60 + mu(3)/20;

int = detJacT*A;

                