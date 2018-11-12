function [mu,rho] = SmoothnessPrior_AutoCorr_EL(RunOptions,Elements,Nodes,DomainIndices,Exp_rho,Var_rho,Corr_rho,ExpShiftExp_rho,ExpShiftVar_rho,Exp_cs,Var_cs,Corr_cs,ExpShiftExp_cs,ExpShiftVar_cs,N_Samples,PLOT)

% SmoothnessPrior_AutoCorr_EL constructs the blocks needed to form the 
% matrix square root of the covariance matrix for the elastic layer
%
% Inputs:
%   Nodes - Coordinates of Nodes
%   DomainIndices - Indicating whether a node is in the elastic layer or not
%   Elements - Matrix where row number corresponds to the element number and the entries of the row are the vertices of the element
%   Exp - Expected Value
%   Var - Variance
%   Corr - Correlation Length
%   ExpShiftExp - Expected value of expected value shift
%   ExpShiftVar - Variance of expected value shift
%   N_Samples - Number of samples to be drawn
%   PLOT - To plot or not to plot, that is the question
%
% Outputs:
%   mu - N_Nodes by N_Samples matrix where each column is a drawn sample of mu
%
% Hwan Goh 30/07/2018, University of Auckland, New Zealand

N_Nodes = length(Nodes); %Number of nodes
N_Elements = size(Elements,1); %Number of elements
Cov_rho = zeros(N_Nodes,N_Nodes);
Cov_cs = zeros(N_Nodes,N_Nodes);

%=== Expected Value and Variance of rho and cs ===%
Exp_rho_Vec = Exp_rho*zeros(N_Nodes,1);
Var_rho_Vec = Var_rho*zeros(N_Nodes,1);
Exp_cs_Vec = Exp_cs*zeros(N_Nodes,1);
Var_cs_Vec = Var_cs*zeros(N_Nodes,1);

for k=1:N_Elements;
    if DomainIndices(k) == 2 %If not a node in the elastic layer
        for jj=1:3
            Exp_rho_Vec(Elements(k,jj)) = Exp_rho;
            Var_rho_Vec(Elements(k,jj)) = Var_rho;
            Exp_cs_Vec(Elements(k,jj)) = Exp_cs;
            Var_cs_Vec(Elements(k,jj)) = Var_cs;
        end
    end
end

%=== Expected Value and Variance of Mean Shifts ===%
ExpShiftExp_rho = ExpShiftExp_rho*ones(N_Nodes,1);
Cov_ExpShift_rho = ExpShiftVar_rho^2*ones(N_Nodes,1)*ones(N_Nodes,1)';
ExpShiftExp_cs = ExpShiftExp_cs*ones(N_Nodes,1);
Cov_ExpShift_cs = ExpShiftVar_rho^2*ones(N_Nodes,1)*ones(N_Nodes,1)';

%% =======================================================================%
%                 Covariance and Samples of rho and cs
%=========================================================================%
%=== Covariance Matrix ===%
for ii=1:N_Nodes
    for jj=1:N_Nodes
        Cov_rho(ii,jj) = Var_rho_Vec(jj)*exp(-(norm(Nodes(ii,:) - Nodes(jj,:),2))^2/(2*Corr_rho^2));
        Cov_rho(jj,ii) = Cov_rho(ii,jj);
        Cov_cs(ii,jj) = Var_cs_Vec(jj)*exp(-(norm(Nodes(ii,:) - Nodes(jj,:),2))^2/(2*Corr_cs^2));
        Cov_cs(jj,ii) = Cov_cs(ii,jj);
    end
end

%=== Cholesky of Covariance ===%
L_rho = chol(Cov_rho + Cov_ExpShift_rho + 10^12*eps*eye(size(Cov_rho)));
L_cs = chol(Cov_cs + Cov_ExpShift_cs + 10^12*eps*eye(size(Cov_cs)));

%=== Plot Samples ===%
if N_Samples==0,
    rho=[];
else
    rho=bsxfun(@plus,Exp_rho_Vec,(L_rho.')*randn(N_Nodes,N_Samples));
    PLOT.TRI=delaunay(Nodes(:,1),Nodes(:,2));
    for ii=1:min(N_Samples,10)
        if PLOT.PriorSamples==1;
            figure
            trisurf(PLOT.TRI,Nodes(:,1),Nodes(:,2),rho(:,ii));
            shading interp
            view(2)
            colorbar
            colormap(jet(256))
            caxis([0 2000])
        end
    end
    cs=bsxfun(@plus,Exp_cs_Vec,(L_cs.')*randn(N_Nodes,N_Samples));
    for ii=1:min(N_Samples,10)
        if PLOT.PriorSamples==1;
            figure
            trisurf(PLOT.TRI,Nodes(:,1),Nodes(:,2),cs(:,ii));
            shading interp
            view(2)
            colorbar
            colormap(jet(256))
            caxis([0 2000])
        end
    end
end

%% =======================================================================%
%                             Samples of mu
%=========================================================================%
%=== Computing Mu ===%
mu = (cs.^2).*rho;

%=== Plot Samples ===%
if N_Samples==0,
    mu=[];
else
    for ii=1:min(N_Samples,10)
        if PLOT.PriorSamples==1;
            figure
            trisurf(PLOT.TRI,Nodes(:,1),Nodes(:,2),mu(:,ii));
            shading interp
            view(2)
            colorbar
            colormap(jet(256))
            caxis([0 5e9])
        end
    end
end

%=== Comparing Samples With True mu ===%
True_mu = zeros(N_Nodes,1);

for k=1:N_Elements;
    if DomainIndices(k) == 2 %If not a node in the elastic layer
        for jj=1:3
            True_mu(Elements(k,jj)) = RunOptions.ElasticLamemu;
        end
    end
end

Diff = zeros(N_Nodes,N_Samples);

for ii=1:N_Samples;
    Diff(:,ii) = abs(mu(:,ii) - True_mu);
    figure
    trisurf(PLOT.TRI,Nodes(:,1),Nodes(:,2),Diff(:,ii));
    shading interp
    view(2)
    colorbar
    colormap(jet(256))
%     caxis([0 1e9])
end














