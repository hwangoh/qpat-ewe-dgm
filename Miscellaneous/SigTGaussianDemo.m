close all
clear all
clc

%=== Gaussian Parameters ===%
% for i = 0:10;
mu = 0.15;
sigma = 0.15; 

%=== Sigmoid Transform Parameters ===%
Prior.SigTalpha = (0+3*0.1)/2; %Parameter of the sigmoid transform
Prior.SigTbeta = (0+3*0.1)/2; %Parameter of the sigmoid transform
Prior.SigTkappa = 1/Prior.SigTbeta; %Parameter of the sigmoid transform
Prior.SigTdelta = (0+3*0.1)/2; %Parameter of the sigmoid transform

%% =======================================================================%
%                         Gaussian Distribution
%=========================================================================
%=== Probability Density Function ===%
x=mu-3*sigma:0.01:mu+3*sigma;
PDF_N=1/(sigma*sqrt(2*pi))*exp(-1/2*((x-mu)/sigma).^2);

%=== Drawing Samples and Plotting ===%
n_Samples=2^20;
Gaussian_Samples=mu+sigma*randn(n_Samples,1);
[Nbins,Nc]=hist(Gaussian_Samples,25);
Ndx=Nc(2)-Nc(1);
normalised_Nbins=Nbins/(n_Samples*Ndx);
figure(1)
bar(Nc, normalised_Nbins)
hold on
plot(x,PDF_N,'r','Linewidth',2)
drawnow
hold off
title(n_Samples)
xlabel('x')
grid on
movegui(figure(1),'west');

%% =======================================================================%
%               Sigmoidal Transformed Gaussian Distribution
%=========================================================================%
%=== Plotting Sigmoidal Transform ===%
y = Prior.SigTalpha + Prior.SigTbeta*tanh(Prior.SigTkappa*(x - Prior.SigTdelta));

figure(2)
plot(x,y,'r')
grid on
hold on
plot(x,x)
hold off
ylim([mu-3*sigma,mu+3*sigma])
title('Sigmoid Transform')

%=== Plotting Sigmoidal Transformed Gaussian Distribution ===%
SigT_Gaussian_Samples = Prior.SigTalpha + Prior.SigTbeta*tanh(Prior.SigTkappa*(Gaussian_Samples - Prior.SigTalpha));
[SigTNbins,SigTNc]=hist(SigT_Gaussian_Samples,25);
SigTNdx=SigTNc(2)-SigTNc(1);
normalised_SigTNbins=SigTNbins/(n_Samples*Ndx);
figure(3)
bar(SigTNc, normalised_SigTNbins)
xlabel('x')
grid on
movegui(figure(3),'east');
title(n_Samples)
hold on
plot(x,PDF_N,'r','Linewidth',2)
hold off

keyboard
% end