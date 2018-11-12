% Jari's cute little trick for letting the means vary. Suppose we have a
% random variable z with mean 0 and covariance Gamma_Z. Now suppose we want
% to take draws of z from z but with varying means. We then instead
% consideer a random variable x = z + q where q is cleverly designed to
% allow this behaviour

close all
clear all
clc

%=== Discretizing Domain ===%
N = 100;
t = [0:N-1]/(N-1);
N_Draws = 5;

%=== Defining and Drawing From Z ===%
lambda = 3;
gamma = exp(-lambda*t); %Ornstein Uhlenbeck process. Correlation decays the further away the points are
Gamma_Z = toeplitz(gamma);
L_Z = chol(Gamma_Z);
Z = L_Z'*randn(N,N_Draws);
figure
plot(t,Z);
title('Draws of Z')

%=== Defining and Drawing From Q ===%
mu_Q = 4*ones(N,1);
sigma_Q = 2;
Gamma_Q = sigma_Q^2*ones(N,1)*ones(N,1)';

[V,D] = eig(Gamma_Q); %Since Gamma_Q rank 1 and not positive definite, we can't chol(Gamma_Q) to get draws. Instead, we simply get another matrix square root via eigendecomposition
L_Q = sqrt(D)*V';
Q = bsxfun(@plus,mu_Q,L_Q'*randn(N,N_Draws));
figure
plot(t,Q);
set(gca,'ylim',[-3, 10]);
title('Draws of Q')

%=== Defining and Drawing From Z ===%
Gamma_X = Gamma_Z + Gamma_Q;
L_X = chol(Gamma_X);
X = bsxfun(@plus,mu_Q,L_X'*randn(N,N_Draws));
figure
plot(t,X);
set(gca,'ylim',[-3, 10]);
title('Draws of X')
