function [x] = JacobiGL(GNalpha,beta,N);

% function [x] = JacobiGL(GNalpha,beta,N)
% Purpose: Compute the N'th order Gauss Lobatto quadrature 
%          points, x, associated with the Jacobi polynomial,
%          of type (GNalpha,beta) > -1 ( <> -0.5). 

x = zeros(N+1,1);
if (N==1) x(1)=-1.0; x(2)=1.0; return; end;

[xint,w] = JacobiGQ(GNalpha+1,beta+1,N-2);
x = [-1, xint', 1]';
return;
