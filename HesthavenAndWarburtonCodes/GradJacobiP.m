function [dP] = GradJacobiP(r, GNalpha, beta, N);

% function [dP] = GradJacobiP(r, GNalpha, beta, N);
% Purpose: Evaluate the derivative of the Jacobi polynomial of type (GNalpha,beta)>-1,
%	       at points r for order N and returns dP[1:length(r))]        

dP = zeros(length(r), 1);
if(N == 0)
  dP(:,:) = 0.0; 
else
  dP = sqrt(N*(N+GNalpha+beta+1))*JacobiP(r(:),GNalpha+1,beta+1, N-1);
end;
return
