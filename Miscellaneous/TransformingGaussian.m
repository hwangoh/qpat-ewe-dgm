close all
clear all
clc

%=== Transformation ===%
LogExpTransform = 1;
SigmoidalTransform = 0;

%=== Initial Guess of p0 ===%
initialguess = 10;

%=== Mean and Standard Deviation ===%
Prior.Exp_p0 = 10;
Prior.InformSmooth_Bounds_p0 = [0,450];
Prior.InformSmooth_STD_p0 = (Prior.InformSmooth_Bounds_p0(2) - Prior.Exp_p0)/2.5; 
% Prior.InformSmooth_STD_p0 = 150;

%=== Log-Exp Transform Parameters ===%
Prior.LogExpTkappa = 0.02;
k = Prior.LogExpTkappa;

%=== Sigmoid Transform Parameters ===%
Prior.SigTalpha = (Prior.Exp_p0+3*Prior.InformSmooth_STD_p0)/2; %Parameter of the sigmoid transform
Prior.SigTbeta = (Prior.Exp_p0+3*Prior.InformSmooth_STD_p0)/2; %Parameter of the sigmoid transform
Prior.SigTkappa = 1/Prior.SigTbeta; %Parameter of the sigmoid transform
Prior.SigTdelta = (Prior.Exp_p0+3*Prior.InformSmooth_STD_p0)/2; %Parameter of the sigmoid transform

%% =======================================================================%
%                         Gaussian Distribution
%=========================================================================%
%=== Probability Density Function ===%
x=Prior.Exp_p0-3*Prior.InformSmooth_STD_p0:0.01:Prior.Exp_p0+3*Prior.InformSmooth_STD_p0;
PDF_N=1/(Prior.InformSmooth_STD_p0*sqrt(2*pi))*exp(-1/2*((x-Prior.Exp_p0)/Prior.InformSmooth_STD_p0).^2);

%=== Drawing Samples and Plotting ===%
n_Samples=2^20;
Gaussian_Samples=Prior.Exp_p0+Prior.InformSmooth_STD_p0*randn(n_Samples,1);
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
grid on
movegui(figure(1),'west');

%% =======================================================================%
%                    Log-Exp Gaussian Distribution
%=========================================================================%
if LogExpTransform == 1
%=== Plotting Sigmoidal Transform ===%
y = (1/k)*log(exp(k*x)+1);
InvLogExpTinitialguess = (1/k)*log(exp(k*initialguess)-1)
Diff = Prior.Exp_p0 - InvLogExpTinitialguess
Slope = exp(k*InvLogExpTinitialguess)./(exp(k*InvLogExpTinitialguess)+1)

figure(2)
plot(x,y,'r')
grid on
hold on
plot(x,x)
hold off
title('Log-Exp Transform')

%=== Plotting Log-Exp Transformed Gaussian Distribution ===%
LogExpT_Gaussian_Samples = (1/k)*log(exp(k*Gaussian_Samples)+1);
[LogExpTNbins,LogExpTNc]=hist(LogExpT_Gaussian_Samples,25);
LogExpTNdx=LogExpTNc(2)-LogExpTNc(1);
normalised_LogExpTNbins=LogExpTNbins/(n_Samples*Ndx);
figure(3)
bar(LogExpTNc, normalised_LogExpTNbins)
grid on
movegui(figure(3),'east');
title(n_Samples)
hold on
plot(x,PDF_N,'r','Linewidth',2)
hold off
end

%% =======================================================================%
%               Sigmoidal Transformed Gaussian Distribution
%=========================================================================%
if SigmoidalTransform == 1
%=== Plotting Sigmoidal Transform ===%
y = Prior.SigTalpha + Prior.SigTbeta*tanh(Prior.SigTkappa*(x - Prior.SigTdelta));
InvSigTinitialguess = (atanh((initialguess - Prior.SigTalpha)/Prior.SigTbeta))/Prior.SigTkappa + Prior.SigTdelta
Diff = Prior.Exp_p0 - InvSigTinitialguess
Slope = Prior.SigTbeta*((sech(Prior.SigTkappa*(InvSigTinitialguess - Prior.SigTdelta))).^2)*Prior.SigTkappa

figure(2)
plot(x,y,'r')
grid on
hold on
plot(x,x)
hold off
title('Sigmoid Transform')

%=== Plotting Sigmoidal Transformed Gaussian Distribution ===%
SigT_Gaussian_Samples = Prior.SigTalpha + Prior.SigTbeta*tanh(Prior.SigTkappa*(Gaussian_Samples - Prior.SigTalpha));
[SigTNbins,SigTNc]=hist(SigT_Gaussian_Samples,25);
SigTNdx=SigTNc(2)-SigTNc(1);
normalised_SigTNbins=SigTNbins/(n_Samples*Ndx);
figure(3)
bar(SigTNc, normalised_SigTNbins)
grid on
movegui(figure(3),'east');
title(n_Samples)
hold on
plot(x,PDF_N,'r','Linewidth',2)
hold off
end