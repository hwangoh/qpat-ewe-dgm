function [Exp, Cov, invCov, traceCov, X_s]=SmoothnessPrior_Informative(Nodes,Elements,MeshDimensions,MinMax,Correlation,Exp,Normalize,N_Samples,PLOT)

% SmoothnessPrior_Informative Implements a proper smoothness prior for Bayesian inverse 
% problems over a rectangle. It is constructed to the conditional
% marginalization of an improper smoothing preprior.
%
%Inputs:
%   Nodes - Co-ordinate matrix
%   Elements - Topology matrix
%   MeshDimensions - Dimensions of Rectangle [length, width]
%   MinMax - 2*1 vector = [sigma_min, sigma_max]
%   Correlation - 2*1 vector = [correlation_x, correlation_y]. The larger the number,
%          the more points near to a marginalisation point are correlated
%          in the x and y direction.
%   Exp -  Expectation of the parameter
%   N_Samples - Number of samples to be drawn and saved in X_s
%
%Outputs:
%   Exp - Expected value of the prior model
%   Cov - Covariance matrix of the parameter 
%   invCov - Inverse of the covariance matrix of the parameter (precision matrix)
%   traceCov - Trace of covariance matrix
%   X_s - N_Nodes by N_Samples matrix where each column is a drawn sample
%
% P. J. Hadwin, University of Auckland, New Zealand
%    26/08/2013 - Adapted from SmoothnessPrior for Hwan
%    24/10/2013 - Modified circular code to obtain this rectangular code
%
% Thanks Paul! 30/01/2018- Edited line 336 to keep variances contained

%Finding marginalisation points for Jari's PrePrior
MargIndices=FindMargPoints(Nodes,Correlation,MeshDimensions);

%Forming the Reg matrix
R=RegularizationMatrix(Nodes,Elements,Correlation);

%Setting parameters for smoothness prior
StdDeviation=[MinMax(2)-Exp,Exp-MinMax(1)];
Variance=[StdDeviation(1)/2.5 StdDeviation(1)/2.5 StdDeviation(2)/5].^2; %A somewhat ad hoc means of defining the variance. By default, just want there to be 2.5 standard deviations from the expected value to the max/min. We also want, by default, 5 standard deviations from the min for the expected value for the background variance.
invCovJoint = R'*R; %Inverse of the joint Covariance Matrix
CovMarg=Variance(1)*eye(length(MargIndices))+Variance(3); %Note: the background variance Vars(3) is only included in here for 'computational reasons' related to sparsity of a matrix. We still wish for the marginalisation points to be independent and identically distributed
N_Nodes=length(Nodes);
Exp=mean(Exp(:))*ones(N_Nodes,1); %

%Smoothness Prior
[Exp, Cov, invCov, traceCov]=preprior(N_Nodes,MargIndices,Exp,invCovJoint,CovMarg,Variance);

disp(['  Variance of Parameters = ' num2str(Variance(1))])
disp(['  Variance of Background = ' num2str(Variance(3))])
disp(['  Trace of Covariance = ' num2str(traceCov)])

%Normalizing Covariance
if Normalize == 1
    Cov_Diag_Inv = sparse(diag(diag(1./Cov).^(1/2)));
    Cov = Cov_Diag_Inv*Cov*Cov_Diag_Inv;
    Std_Diag = sqrt(Variance(1))*speye(N_Nodes);
    Cov = Std_Diag*Cov*Std_Diag;
    traceCov = trace(Cov);
    invCov = inv(Cov);
end

%Plot Marginalisation Points
if PLOT.PriorMargPoints==1;
    figure(PLOT.Figure_PriorMargPoints)
    FemPlot2D(Nodes,Elements),hold on
    plot(Nodes(MargIndices,1),Nodes(MargIndices,2),'ro')
    drawnow
%     title(PLOT.Figure_Prior_Title,'FontWeight','bold') 
    pause(0.1)
end

%Plot Samples
if N_Samples==0,
    X_s=[];
else
    X_s=bsxfun(@plus,Exp,(chol(Cov).')*randn(length(Nodes),N_Samples));
    X_s(X_s<0) = 0; %Positivity constraint
    if PLOT.PriorSamples==1;
        PLOT.TRI=delaunay(Nodes(:,1),Nodes(:,2));
        for ii=1:N_Samples
            figure
            trisurf(PLOT.TRI,Nodes(:,1),Nodes(:,2),X_s(:,ii));
            shading interp
            view(2)
            colorbar
            colormap(jet(256))
            caxis([0 600])
            zlim([0 600])
        end
    end
%     if PLOT.PriorSamples==1;
%         kappa = 0.08;
%         LogExpTX_s = (1/kappa)*log(exp(kappa*X_s)+1);
%         for ii=1:N_Samples
%             figure
%             trisurf(PLOT.TRI,Nodes(:,1),Nodes(:,2),LogExpTX_s(:,ii));
%             shading interp
%             view(2)
%             colorbar
%             colormap(jet(256))
%             caxis([0 600])
%             zlim([0 600])
%         end
%     end
%    if PLOT.PriorSamples==1;
%         figure
%         alpha = (220 + 7)/2;
%         beta = (220 - 7)/2;
%         kappa = 1/beta;
%         SigTX_s = alpha + beta*tanh(kappa*(X_s - alpha));
%         FemPlot2D(Nodes,Elements,SigTX_s(:,1:min(N_Samples,10))) %Only plots up to 10
%     end
end

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%

function MargIndices=FindMargPoints(Nodes,Correlation,MeshDimensions)

%Defines the marginalisation points.
%Inputs:
%   Nodes - Co-ordinate matrix
%   Correlation - 2*1 vector = [correlation_x, correlation_y]. The larger the number,
%          the more points near to a marginalisation point are correlated
%          in the x and y direction. For example, if corr(1) is large, we 
%          will obtain less marginalisation points in the x direction since
%          more points between 2 marginalisation points will be correlated
%          and we wish for marginalisation points to be independent of each
%          other.
%   MeshDimensions - Dimensions of Rectangle [length, width]
%
%Outputs:
%   MargIndices - index of nodes that act as the marginalisation points
%
% P. J. Hadwin, this preamble was written by Hwan Goh who does not know the
% exact date P wrote this code. University of Auckland

corr_x=Correlation(1);
corr_y=Correlation(2);
xm=min(Nodes(:,1)); xM=max(Nodes(:,1)); indx=xm:corr_x:xM;
ym=min(Nodes(:,2)); yM=max(Nodes(:,2)); indy=ym:corr_y:yM;

if ~(indx(end)==xM) %if the last index is not equal to the x-edge the domain
    indx=indx+0.5*(xM-indx(end)); %shifts indx so that |xm - min(Nodes(:,1))| = |xM - max(Nodes(:,2))|
end 
if ~(indy(end)==yM) %if the last index is not equal to the y-edge the domain
    indy=indy+0.5*(yM-indy(end)); %shifts indy so that |ym - min(Nodes(:,1))| = |yM - max(Nodes(:,2))|
end 

[X,Y]=meshgrid(indx,indy); %Creates mesh with quadrilateral elements. The grid points have coordinates formed by combinations of indx and indy
x=X(:); 
y=Y(:); 

x_ok_ind = find(abs(x)<(0.96*MeshDimensions(1))); %Find grid points that are within 0.96 times the boundaries of the mesh
y_ok_ind = find(abs(y)<(0.96*MeshDimensions(2))); %Find grid points that are within 0.96 times the boundaries of the mesh
x = x(x_ok_ind); y = y(y_ok_ind);

MargIndices = NodeSearch(Nodes,[x y]);

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
function MargIndices=NodeSearch(Nodes,X)

% NodeSearch associates the nodes on the computational mesh with the
% selected marginalization points on the quadrilateral grid. Specifically,
% it finds the index of the nodes best suited to serve as a marginalisation
% point.
% NodeSearch does the same thing as dsearch however
% dsearch is being removed and the replacement you need a
% Delaunay object. Whereas NodeSearch only needs the co-ordinates. 
%
%Inputs:
%   Nodes - Coordinates of Nodes where each row defines a new node
%   X - Coordinates of Points of interest where each row is a new point
%
%Outputs:
%   MargIndices - indices of the nodes closest to the points in X
%
% P. J. Hadwin, University of Auckland, New Zealand
%    14/02/2012 - Original

if ~(size(Nodes,2)==size(X,2))
    error('NodeSearch:dims1',['Dimension missmatch with Nodes and X. Make sure ' ...
                            'they have the same number of columns.'])
end
if size(Nodes,1)<size(X,1)
	printf('\nNodeSearch Note:')
    printf(['You are trying to find more Nodes that are avalible. ' ...
        '\nCheck the inputs are in the correct order.'])
end
MargIndices=zeros(size(X,1),1);
for ii=1:size(X,1)
    dist=sqrt(sum(abs(bsxfun(@plus, -X(ii,:),Nodes)).^2,2)); %Distances of each node on the mesh with each marginalization points: sqrt((n_x - X_ii)^2 + (n_y - X_ii)^2)
    [~, MargIndices(ii)]=min(dist);
end

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%

function R=RegularizationMatrix(Nodes,Elements,Correlation)

% Forms a 2*N_Elm by N_Nodes regularization matrix R. Later, we use this to
% form the N_Nodes by N_Nodes inverse of the joint covariance matrix
% invCovJoint = R'*R.
%
% Rows of the matrix R correspond to components of the gradient
% of a linear function in all the elements of the grid.
%
% The gradients are weighted with the area of the corresponding element
% so that the product R*x gives the integrals of the components of the
% gradient in each element
%
% Hwan's note - "smaller triangles are given larger weight. Nodes that are
%                closer together are more correlated"
%
% Aku Seppänen 5.10.2001
% Aku Seppänen 6.12.2001 % weighting with the areas of the triagles

Lambda = diag([Correlation(1),Correlation(2)]);

R = zeros(2*length(Elements),length(Nodes));
L = [-1 0 1;-1 1 0];
for ii = 1:length(Elements)
    x1 = Nodes(Elements(ii,1),1); y1 = Nodes(Elements(ii,1),2);
    x2 = Nodes(Elements(ii,2),1); y2 = Nodes(Elements(ii,2),2);
    x3 = Nodes(Elements(ii,3),1); y3 = Nodes(Elements(ii,3),2);
    Jac_T = [x3-x1 y3-y1; x2-x1 y2-y1];
    Q = Jac_T\L; %This is J_{F_K}^-T \nabla N_j
    Q = sqrt(abs(.5*det([x1 y1 1; x2 y2 1; x3 y3 1])))*Q;
    R(2*ii-1:2*ii,Elements(ii,:)) = Lambda*Q; %Fills in 2x3 entry of R. It is 2*ii-1:2*ii instead of ii:ii+1 because we skip rows in pairs
end

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
function [Exp, Cov, invCov, traceCov]=preprior(N_Nodes, MargIndices, Exp, invCovJoint, CovMarg, Variance)

%preprior Implements a proper smoothness prior for Bayesian inverse problems. 
%           Uses Variance equalisation, background variation and mean GNalpha
%           criteria. If Vars(2) = 0 then Variance equalisation not used.
%           If Vars(3) = 0 then background variation not used. 
%        
%Inputs:
%   N_Node - Number of Nodes in the Mesh
%   MargIndices - Indicies of the marginalisation points
%   Exp - Prior guess for the estimate - \mu_X
%   invCovJoint - Inverse of the joint Covariance Matrix - D'D
%   CovMarg - Covariance of the marginisation points - \Gamma_Y
%   Variance - vector of different variances [Margin Points, Other Points, Background]
%
%Outputs:
%   Exp - Updated Prior guess of the estimate - \mu_{X_{pre}}
%   Cov - PrePrior Covariance matrix \Gamma_{X_{pre}}. 
%   invCov - inverse of the PrePrior Cov Matrix
%
% P. J. Hadwin, University of Auckland, New Zealand
%    23/08/2011 - Original. Method in Jari's Book and Arridge2006
%    31/10/2011 - added variance equalisation, background variation and
%                  mean alpha criteria

notMargIndices=setdiff((1:N_Nodes)',MargIndices);

%blocks of squared reg matrix L'L
B12=invCovJoint(notMargIndices,MargIndices);
B11=invCovJoint(notMargIndices,notMargIndices);
invB11B12 = (B11\B12);

%mean GNalpha criterion
ss=Variance(2); 
if ss==0 
    ss=Variance(1); 
end
m1 = trace(inv(B11)); 
m2 = trace(invB11B12*CovMarg*invB11B12'); 
nd = length(notMargIndices);
alpha=min([abs(m1/(nd*ss - m2)),20000]); %Changed this from 20000, not sure what it does but it made my variances go into the thousands
disp(['  alpha = ' num2str(alpha)])

%inverse cov matrix for preprior
invCov=sparse(N_Nodes,N_Nodes);
invCov(notMargIndices,notMargIndices)=alpha*B11;
invCov(notMargIndices,MargIndices)=alpha*B12;
invCov(MargIndices,notMargIndices)=alpha*(B12.');
invCov(MargIndices,MargIndices)=alpha*(B12.')*invB11B12+inv(CovMarg);

%expectation
Exp(notMargIndices)=-invB11B12*Exp(MargIndices);

%Outputs
Cov=inv(invCov);
traceCov = trace(Cov);

