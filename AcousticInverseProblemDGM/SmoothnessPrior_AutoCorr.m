function [invL,traceCov,Cov,X_s] = SmoothnessPrior_AutoCorr(Nodes,Exp,Var,Corr,N_Samples,PLOT)

% SmoothnessPrior_AutoCorr constructs the blocks needed to form the 
% matrix square root of the covariance matrix.
%
% Inputs:
%   Nodes - Coordinates of Nodes
%   Elements - Matrix where row number corresponds to the element number and the entries of the row are the vertices of the element
%   Exp - Expected Value
%   Var - Variance
%   Corr - Correlation Length
%
% Outputs:
%   L - Matrix square root of the covariance matrix
%   traceCov - Trace of covariance matrix
%   X_s - N_Nodes by N_Samples matrix where each column is a drawn sample
%
% Hwan Goh 01/03/2018, University of Auckland, New Zealand

N_Nodes = length(Nodes); %No. of nodes
Cov = zeros(N_Nodes,N_Nodes);

%=== Covariance Matrix ===%
for ii=1:N_Nodes
    for jj=1:N_Nodes
        Cov(ii,jj) = Var*exp(-(norm(Nodes(ii,:) - Nodes(jj,:),2))^2/(2*Corr^2));
        Cov(jj,ii) = Cov(ii,jj);
    end
end

%=== Normalizing Covariance ===%
% if Normalize == 1
%     Cov_Diag_Inv = sparse(diag(diag(1./Cov).^(1/2)));
%     Cov = Cov_Diag_Inv*Cov*Cov_Diag_Inv;
%     Std_Diag = sqrt(Var)*speye(N_Nodes);
%     Cov = Std_Diag*Cov*Std_Diag;
% end

%=== Cholesky of Covariance ===%
traceCov = trace(Cov);

% L = chol(Cov); %L'*L = Cov
% invCov = inv(Cov);
% invL = chol(invCov);

L = chol(Cov+10^10*eps*eye(size(Cov)));
invL = inv(L);
invL = invL';

% L = chol(Cov+10^10*eps*eye(size(Cov)));
% invCov = inv(Cov+1e-2*eye(size(Cov)));
% invL = chol(invCov);

%=== Plot Samples ===%
if N_Samples==0,
    X_s=[];
else
    X_s=bsxfun(@plus,Exp,(L.')*randn(N_Nodes,N_Samples));
    X_s(X_s<0) = 0; %Positivity constraint
    PLOT.TRI=delaunay(Nodes(:,1),Nodes(:,2));
    for ii=1:min(N_Samples,10)
        if PLOT.PriorSamples==1;
            figure
            trisurf(PLOT.TRI,Nodes(:,1),Nodes(:,2),X_s(:,ii));
            shading interp
            view(2)
            colorbar
            colormap(jet(256))
            caxis([0 500])
            zlim([0 500])
%             caxis([0 0.3])
%             zlim([0,0.25])
        end
    end
%     k = 0.02;
%     X_sT=(1/k)*log(exp(k*X_s)+1);
%     PLOT.TRI=delaunay(Nodes(:,1),Nodes(:,2));
%     for ii=1:min(N_Samples,10)
%         if PLOT.PriorSamples==1;
%             figure
%             trisurf(PLOT.TRI,Nodes(:,1),Nodes(:,2),X_sT(:,ii));
%             shading interp
%             view(2)
%             colorbar
%             colormap(jet(256))
%             caxis([0 500])
% %             caxis([0 0.25])
%         end
%     end
end














