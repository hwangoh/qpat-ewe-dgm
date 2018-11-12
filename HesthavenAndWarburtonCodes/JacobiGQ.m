function [x,w] = JacobiGQ(GNalpha,beta,N);

% function [x,w] = JacobiGQ(GNalpha,beta,N)
% Purpose: Compute the N'th order Gauss quadrature points, x, 
%          and weights, w, associated with the Jacobi 
%          polynomial, of type (GNalpha,beta) > -1 ( <> -0.5).

if (N==0) x(1)= -(GNalpha-beta)/(GNalpha+beta+2); w(1) = 2; return; end;

% Form symmetric matrix from recurrence.
J = zeros(N+1);
h1 = 2*(0:N)+GNalpha+beta;
J = diag(-1/2*(GNalpha^2-beta^2)./(h1+2)./h1) + ...
    diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+GNalpha+beta).*...
    ((1:N)+GNalpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
if (GNalpha+beta<10*eps) J(1,1)=0.0;end;
J = J + J';

%Compute quadrature by eigenvalue solve
[V,D] = eig(J); x = diag(D);
w = (V(1,:)').^2*2^(GNalpha+beta+1)/(GNalpha+beta+1)*gamma(GNalpha+1)*...
    gamma(beta+1)/gamma(GNalpha+beta+1);
return;
