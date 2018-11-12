function Jacmu_aphi = ConstructJacmuaphiAbs(Mesh,Prmtrs,phi,A)

% ConstructJacmuaphiAbs constructs the Jacobian of mu_aphi for Gauss-Newton
% iterations when the scattering coefficient is known
%
% Inputs:
%   Mesh:
%      Mesh.N_Nodes - Number of nodes
%      Mesh.Nodes - Coordinates of Nodes
%      Mesh.N_Elm - Number of elements
%      Mesh.Elements - Matrix where row number corresponds to the element
%      number and the entries of the row are the vertices of the element
%   Prmtrs:
%      Prmtrs.kappa - N_Nodes by 1 vector for the diffusion coefficient
%      Prmtrs.mu_a - N_Nodes by 1 vector for the absorption coefficient
%   phi - ODTForward solution using the current value of mu_a
%   A - S+M+R where S,M and R are the system matrices from ODTForward
%
% Outputs:
%   Jacmu_aphi - Jacobian of mu_aphi; an N by N matrix 
%
% Hwan Goh 01/9/2013, University of Auckland, New Zealand

%Shortening labels
N_Nodes = Mesh.N_Nodes;
Nodes = Mesh.Nodes;
Elements = Mesh.Elements;
kappa = Prmtrs.kappa;
mu_a = Prmtrs.mu_a;

%%%%%%%%%%%%%%%%%%%%%
%%% Gradient of A %%%
%%%%%%%%%%%%%%%%%%%%%
gradA = sparse(N_Nodes^2,N_Nodes);
for i=1:N_Nodes;
    partialA = sparse(N_Nodes,N_Nodes);
    [F,~]=find(i==Elements);
    for j=1:size(F,1); 
        Ver=Elements(F(j),:); 
        Coords=Nodes(Ver,:); 
        %Stiffness Matrix
        partialA(Ver,Ver) = partialA(Ver,Ver)+ pwlpartialEStiffmu_a(Coords,kappa(i));
        %Mass Matrix
        pos = find(Ver==i); %finds position of i in a row of the 
                            %Elements matrix
        partialA(Ver,Ver) = partialA(Ver,Ver) + pwlpartialEMassmu_a(Coords,pos);
    end
    gradA(:,i) = partialA(:);
end

%%%%%%%%%%%%%%%%%%%
%%% (gradA)*phi %%%
%%%%%%%%%%%%%%%%%%%
%The creation of gradA*phi ensures that we perform matrix inversion only 
%once, instead of an N_Node number of times
gradAphi = gradA*phi;
gradAphi = reshape(gradAphi,N_Nodes,N_Nodes);
% L = chol(A);
Ainv_gradAphi = A\gradAphi; 

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gradient of mu_a %%%
%%%%%%%%%%%%%%%%%%%%%%%%
gradmu_a = eye(N_Nodes,N_Nodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Jacobian of mu_aphi %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jacmu_aphi = zeros(N_Nodes,N_Nodes);
for i=1:N_Nodes;
    Jacmu_aphi(:,i) = -mu_a.*Ainv_gradAphi(:,i) + gradmu_a(:,i).*phi;
end