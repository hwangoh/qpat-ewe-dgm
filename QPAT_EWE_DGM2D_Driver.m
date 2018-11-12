% QPAT_EWEDGM2D_Driver.m is the main script file for the quantitative photoacoustic
% tomography inverse problem using the elastic wave equation and discontinuous Galerkin method.
%
% % Hwan Goh, 17/09/2017, University of Auckland, New Zealand

%=== Some useful command lines for running this code over the server ===%
% aitken.math.auckland.ac.nz
% ssh maclaurin.math
% addpath('QPAT_EWE_DGM2D_Codes_Hwan')
% run QPAT_EWEDGM2D_Driver
% cd /gpfs1m/projects/nesi00452/QPAT_EWE_DGM2D_Codes_Hwan

close all
clear all
clc
restoredefaultpath
format long
warning('off','all')
addpath(genpath('../QPAT_EWE_DGM2D_Codes_Hwan'))

%adding paths one by one because the servers don't like the genpath command 
addpath('QPAT_EWE_DGM2D_Codes_Hwan')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/AcousticForwardProblemDGM')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/AcousticInverseProblemDGM')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/HesthavenAndWarburtonCodes')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/MATFiles')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/MeshGenerationAndInterpolation')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/Miscellaneous')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/OpticalForwardProblemFEM')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/OpticalInverseProblemFEM')
addpath('QPAT_EWE_DGM2D_Codes_Hwan/ParameterGeneration')

%% =======================================================================%
%                        Demonstration Properties
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Parameters and Mesh or Load Existing? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Parameters ===%
RunOptions.GeneratePrmtrs = 0; %Whether or not to generate mesh
RunOptions.GenerateCentralGaussian = 0; %Generate centrally located blob or random blobs
RunOptions.GenerateAwayFromBoundary = 0; %Generate one chromophore away from boundary
RunOptions.CentralGaussianSize = 4e-6; %Diameter of the central blob, see GenPrmtrsI line 54. Default is 4e-6
RunOptions.ScalingOptical = 1; %10^3; %Scale that the optical mesh is in, default is [m] so e.g. 10^3 is [mm]. Problem with marginalisation points of prior when optical mesh is in [m] (ScalingOptical = 1)
%=== Mesh ===%
RunOptions.GenerateMesh = 0; %Whether or not to generate mesh

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Obtained Reconstruction %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunOptions.SaveReconstruction = 0;

%% %%%%%%%%%%%%%%%%%%%%
%%% Mesh Properties %%%
%%%%%%%%%%%%%%%%%%%%%%%
%=== Which Mesh ===%
RunOptions.UseTrelisMesh = 1; %Use Trelis generated mesh, set these properties in Trelis_Skull2D_JournalGenerator.m
RunOptions.UseMyRectangularMesh = 0; %Use my mesh generation code. If using pre-generated parameters of my mesh, make sure this is 1!

%=== Sensor and Lights Source Properties ===%
RunOptions.NumberofSensorsOnOneBoundaryEdge = 10; %Number of sensors on one boundary edge of the domain
RunOptions.UseFullDomainData = 0; %Output of forward problem sensory data consisting of the full domain; still set RunOptions.FullqVectorData = 1
RunOptions.FullqVectorData = 0; %Output of forward problem sensory data of full q vector
RunOptions.VelocitiesData = 1; %Output of forward problem sensory data of velocities, only works when sensors are placed on the boundary
RunOptions.OptI_s=1; %Diffuse boundary condition - "Input strength of 1"- 'The finite element method for Propogation...,' M.Schweiger

%=== Trelis Mesh Properties ===%
RunOptions.InverseCrime = 0; %For using the Trelis mesh, whether we are commiting an inverse crime must be manually specified since I have not found a clear metric by which this can be determined (Not number of elements or nodes).
RunOptions.TrelisMeshDElementSize = '0003'; %Entry for generating Trelis data mesh, main purpose is for the file name when saving
RunOptions.TrelisMeshIElementSize = '0007'; %Entry for generating Trelis inverse mesh, main purpose is for the file name when saving
RunOptions.TrelisMeshAElementSize = '00032'; %Entry for generating Trelis accurate mesh, main purpose is for the file name when saving
RunOptions.BoundaryCondition = 'Neumann';

%=== Rectangular Mesh Generation Properties ===%
if RunOptions.UseMyRectangularMesh == 1
%Defining Nodes
RunOptions.DefineNodesVertically = 1; %Nodes are indexed vertically upwards from bottom left corner (vertical or horizontal doesn't seem to make a difference)
RunOptions.DefineNodesHorizontally = 0; %Nodes are indexed horizontally to the right from bottom left corner (vertical or horizontal doesn't seem to make a difference)
RunOptions.RotateCornerElementsToOneElement = 0; %Rotate the corner elements so that every corner contains only one element
RunOptions.RotateCornerElementsToTwoElements = 0; %Rotate the corner elements so that every corner is shared by two elements
RunOptions.BndExtndx = 0; %total extension of boundary in +/-x direction [m] default=0.015
RunOptions.BndExtndy = 0; %total extension of boundary in +/-y direction [m] default=0.015
%Boundary Condition
RunOptions.TopBC = 'Neumann'; %Homogeneous Boundary conditions for DGM, top edge
RunOptions.LeftBC = 'Neumann'; %Homogeneous Boundary conditions for DGM, left edge
RunOptions.RightBC = 'Neumann'; %Homogeneous Boundary conditions for DGM, right edge
RunOptions.BottomBC = 'Neumann'; %Homogeneous Boundary conditions for DGM, bottom edge
if RunOptions.GenerateMesh == 1
%Data Mesh Properties
MeshD.Dimensns = (1/2)*[0.02 0.02]; %Width and height of the rectangle, default is in metres
MeshD.N_BndElmx = 5; %Number of boundary elements in 0.005m for top and bottom edges, default = 5
MeshD.N_BndElmy = 5; %Number of boundary elements in 0.005m for left and right edges, default = 5
%Inversion Mesh Properties
MeshI.Dimensns=[MeshD.Dimensns(1) MeshD.Dimensns(2)]; %Width and height of the rectangle from center  
MeshI.N_BndElmx = 3; %Number of boundary elements in 0.005m for top and bottom edges, default = 3
MeshI.N_BndElmy = 3; %Number of boundary elements in 0.005m for left and right edges, default = 3
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Optical Forward Problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%=== Acoustic Parameters ===% %More optical parameters in QPATPrmtrs2D, but these do not need to be changed as often
PrmtrsI.mu_sEst = 2; %Constant Value for the Scattering Coefficient, [mm^-1]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acoustic Forward Problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Domain Type ===%
RunOptions.FluidDomainWithSolidLayerMeshD = 1; %Use Fluid domain with solid layer representing the skull for data mesh
RunOptions.FluidDomainWithSolidLayerMeshI = 0; %Use Fluid domain with solid layer representing the skull for inverse mesh
RunOptions.FluidMeshD = 0; %Use purely fluid domain for data mesh
RunOptions.FluidMeshI = 1; %Use purely fluid domain for inverse mesh
RunOptions.SolidMeshD = 0; %Use purely solid domain for data mesh
RunOptions.SolidMeshI = 0; %Use purely solid domain for inverse mesh

%=== Acoustic Parameters ===% (Mitsuhasi: In fluid: cp = 1500 m/s, cs = 0. In skull, cp = 3000 m/s and cs = 1500 m/s.)
RunOptions.TestAcousticParameters = 0; %Use test acoustic parameters
RunOptions.RealAcousticParameters = 1; %Use realistic acoustic parameters
if RunOptions.TestAcousticParameters == 1;
RunOptions.AcousticMediumWaveSpeed_cp = sqrt(2);
RunOptions.ElasticMediumWaveSpeed_cp = sqrt(2);
RunOptions.ElasticMediumWaveSpeed_cs = sqrt(1/2);
RunOptions.AcousticMediumDensity = 1;
RunOptions.ElasticMediumDensity = 1;
RunOptions.ElasticLamemu = (RunOptions.ElasticMediumWaveSpeed_cs^2)*RunOptions.ElasticMediumDensity; %also known as the shear modulus, cs = sqrt(mu/rho) => mu = cs^2*rho
RunOptions.ElasticLamelambda = RunOptions.ElasticMediumWaveSpeed_cp^2*RunOptions.ElasticMediumDensity - 2*RunOptions.ElasticLamemu; %cp = sqrt{(lambda + 2mu)/rho)} => lambda = cp^2*rho - 2*mu
end
if RunOptions.RealAcousticParameters == 1;
RunOptions.AcousticMediumWaveSpeed_cp = 1500; %default is 1500 m/s: "Mitsuhashi: A forward-adjoint pair based on the elastic wave equation for use in transcranial photoacoustic computed tomography"
RunOptions.ElasticMediumWaveSpeed_cp = 3000; %default is 3000 m/s: "Mitsuhashi: A forward-adjoint pair based on the elastic wave equation for use in transcranial photoacoustic computed tomography"
RunOptions.ElasticMediumWaveSpeed_cs = 1500; %default is 1500 m/s: "Mitsuhashi: A forward-adjoint pair based on the elastic wave equation for use in transcranial photoacoustic computed tomography"
RunOptions.AcousticMediumDensity = 1000; %default is 1000 kg/m^3: "Mitsuhashi: A forward-adjoint pair based on the elastic wave equation for use in transcranial photoacoustic computed tomography"
RunOptions.ElasticMediumDensity = 1850; %default is 1850 kg/m^3: "Mitsuhashi: A forward-adjoint pair based on the elastic wave equation for use in transcranial photoacoustic computed tomography"
RunOptions.ElasticLamemu = (RunOptions.ElasticMediumWaveSpeed_cs^2)*RunOptions.ElasticMediumDensity; %also known as the shear modulus, cs = sqrt(mu/rho) => mu = cs^2*rho
RunOptions.ElasticLamelambda = RunOptions.ElasticMediumWaveSpeed_cp^2*RunOptions.ElasticMediumDensity - 2*RunOptions.ElasticLamemu; %cp = sqrt{(lambda + 2mu)/rho)} => lambda = cp^2*rho - 2*mu
end
RunOptions.LinearHeatExpansionCoeff = 1;
RunOptions.SpecificHeatCoeff = 1;

%=== DGM Properties ===%
%=== Initial Conditions ===%
RunOptions.StrainVelocityForm = 1; %Strain-Velocity form of the elastic wave equation
RunOptions.StressVelocityForm = 0; %Stress-Velocity form of the elastic wave equation
RunOptions.InterplteGaussianp0 = 0; %Test by Gaussian initial condition interpolated on DGM Mesh
RunOptions.p0Test1 = 0; %Test case 1
RunOptions.p0Test2 = 0; %Test case 2
RunOptions.GaussianWidth = 1e5;%(1e11); %Increase to decrease width of the Gaussian
%=== RHS Computation ===%
RunOptions.DGMPolyOrder = 2; %Polynomial order used for approximation

%=== Time-Stepping Properties ===%
RunOptions.TimeLSERK4 = 1; %Use low-storage five-stage explicit fourth order Runge-Kutta method
RunOptions.TimeStepUseCFLTimo = 1/2; %Set CFL condition for Timo's time step size, set to 0 if you don't want to use
RunOptions.TimeStepSizeLSERK4 = 1e-6; %Set own value for time step size for LSERK4 time stepping, set to 0 if you don't want to use
if RunOptions.TestAcousticParameters == 1;
    RunOptions.FinalTime = 0.008; %Final Time
end
if RunOptions.RealAcousticParameters == 1;
    RunOptions.FinalTime = 0.000004; %Final Time
end

%% %%%%%%%%%%%%%%%%%
%%% Adding noise %%%
%%%%%%%%%%%%%%%%%%%%
RunOptions.AddNoise = 1; %Add noise?
RunOptions.NoiseLevel = 0.001; %Scalar multiplier of noise draws
RunOptions.NoiseMinMax = 0; %Use max minus min of data
RunOptions.NoiseMinMaxS = 0; %Use max minus min of data at each sensor
RunOptions.NoiseMax = 1; %Use max of all data
RunOptions.NoiseMaxS = 0; %Use max of data at each sensor
RunOptions.LoadNoisyData = 1; %Load noisy data, otherwise generates and saves

%% %%%%%%%%%%%%%%%%%%%%%%
%%% Inversion Method  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%=== Uncertainty Estimates ===%
RunOptions.ApproximatePosteriorCovariance = 0; %Form approximation of the posterior covariance
%=== Acoustic Parameters ===%
RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensity = 1; %Reconstruct absorbed energy density
RunOptions.EWE_DGM2D_ReconstructAbsorbedEnergyDensityHyperPrior = 0; %Reconstruct absorbed energy density using hyperprior
%=== Noise Model ===%
RunOptions.Cov_ENoiseLevel = RunOptions.NoiseLevel; %Scalar multiplier of diagonal for conventional noise model
RunOptions.AE_EWE = 1; %Noise model accounts for approximation error due to discretization
RunOptions.AE_DA = 0; %Noise model accounts for approximation error due to discretization
RunOptions.AEN_Samples = 1000; %Number of samples for approximation error for acoustic inverse problem. Default is 2500 from "Approximation errors and model reduction with an application in optical diffusion tomography - Arridge et al 2006"
RunOptions.AE_DA_N_Samples = 1000; %Number of samples for approximation error for optical inverse problem.
RunOptions.GenerateAndSave_AESamples = 0; %Generate and save samples
RunOptions.AE_VarySolidLayer = 0; %Draw elastic layer from prior models
RunOptions.AE_ConstructStats_SVD = 1; %Generate and save approximation error statistics using SVD
RunOptions.AE_ConstructStats_QR = 0; %Generate and save rank reduced approximation error statistics using QR
RunOptions.Save_AEStats = 1; %Generate and save approximation error statistics
RunOptions.Load_AESampleMeans = 0; %Load approximation error sample means
RunOptions.Load_AESampleCovs = 0; %Load approximation error sample covariances
RunOptions.Save_AEErrorModel = 1; %Generate and save AE enchanced error model
RunOptions.Load_AEErrorModel = 0; %Load AE enchanced error model
%=== Acoustic Inverse Line Search Options ===%
RunOptions.EWE_DGM2D_UseSteepestDescent = 0; %Use steepest descent direction
RunOptions.EWE_DGM2D_UseNewton = 0; %Use Newton direction
RunOptions.EWE_DGM2D_UseGaussNewton = 1; %Use Gauss-Newton direction
RunOptions.EWE_LS_hInitialGuess = 10; %Initial guess constant
RunOptions.EWE_LS_MaxItrtn = 10; %Maximum number of LS iterations
RunOptions.EWE_LS_RemoveTimesSteps = 0; %Number of the first times steps to neglect
RunOptions.EWE_LS_ErrTol = 1e-10; %Default 1e-10
RunOptions.EWE_LS_LogExpTransform = 1; %Use Log-Exp transform positivity constraint
RunOptions.EWE_LS_SigmoidTransform = 0; %Use sigmoid transform positivity constraint
RunOptions.EWE_LS_ExponentialTransform = 0; %Use exponential transform positivity constraint
RunOptions.EWE_LS_Analytical = 0; %Photoacoustic wave propagation is linear and so we can analytically compute the stepsize
RunOptions.EWE_LS_fminsearch = 1; %Use fminsearch to compute stepsize
%=== Optical Inverse Gauss-Newton Options ===%
RunOptions.DA_GN_CoupledErrorModel = 1; %Use Coupled Error Model
RunOptions.DA_GN_MaxItrtn = 100; %Maximum number of GN iterations
RunOptions.DA_GN_ErrTol = 1e-10; %Default 1e-10
RunOptions.DA_GN_LS_fminsearch = 1; %Use fminsearch to compute stepsize
RunOptions.DA_GN_LS_PaulsFunction = 0; %Use Paul's LS to compute stepsize for optical inverse problem
RunOptions.DA_GN_LS_StepSize = 0.5; %Initial stepsize for Paul's line search function
RunOptions.DA_GN_LS_MaxStepSize = 2; %Max step size permitted for Paul's line search function
RunOptions.DA_GN_LS_MaxItrtn = 8; %Max number of iterations for Paul's line search function
RunOptions.DA_GN_LS_PostivityConstraint = 1; %Whether positivity constraints are implemented for Paul's line search function
RunOptions.DA_GN_LS_mu_aMin = 1e-8;  %Minimum mu_a permitted for Paul's line search function

%% %%%%%%%%%%%%%%%%%%%%%
%%% Prior Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%
Prior.N_Samples = 5;
%=== Absorbed Energy Density ===%
RunOptions.UseInformSmoothPrior = 1; %Use informative smoothness prior
RunOptions.UseACPrior = 0; %Use autocorrelation prior
Prior.Hyperprior_alpha = 1e-3; %Initial guess of alpha
Prior.Hyperprior_Exp_alpha = 1e-4; %expected value of alpha
Prior.Hyperprior_Std_alpha = sqrt(-(2*(Prior.Hyperprior_Exp_alpha-Prior.Hyperprior_alpha))/(125/Prior.Hyperprior_alpha)); %standard deviation of alpha
Prior.Exp_h = 70; %Expected Value of h, Default = 70
Prior.InformSmooth_Normalize = 0; %Normalize the covariance
Prior.InformSmooth_Bounds_h = [0,450]; % approximately [h_min, h_max], Default = [0, 450]
Prior.InformSmooth_Corr_h = [0.008 0.008]*RunOptions.ScalingOptical; %[correlation_x, correlation_y]
Prior.InformSmooth_STD_h = (Prior.InformSmooth_Bounds_h(2) - Prior.Exp_h)/2.5; 
Prior.AC_Var_h = 150^2; %Variance of p_0
Prior.AC_Corr_h = 0.001; %Correlation Length of h, Default: 0.0001;
%=== Prior for Elastic Layer ===%
Prior.Exp_rho_e = 1850; %Expected Value of rho_e
Prior.AC_Var_rho_e = 0; %Variance of rho_e
Prior.AC_Corr_rho_e = 0.00001; %Correlation Length of rho_e
Prior.AC_ExpShiftExp_rho_e = 10; %Mean of mean shift of rho_e
Prior.AC_ExpShiftVar_rho_e = 10^2; %Variance of mean shift of rho_e
Prior.Exp_cs = 1500; %Expected Value of cs
Prior.AC_Var_cs = 0; %Variance of cs
Prior.AC_Corr_cs = 0.00001; %Correlation Length of cs
Prior.AC_ExpShiftExp_cs = 10; %Mean of mean shift of rho_e
Prior.AC_ExpShiftVar_cs = 10^2; %Variance of mean shift of rho_e
%=== Positivity Transform ===%
Prior.LogExpTkappa = 0.08; %Parameter of the log-exp transform
Prior.SigTalpha = (Prior.Exp_h+3*Prior.InformSmooth_STD_h)/2; %Parameter of the sigmoid transform
Prior.SigTbeta = (Prior.Exp_h+3*Prior.InformSmooth_STD_h)/2; %Parameter of the sigmoid transform
Prior.SigTkappa = 1/Prior.SigTbeta; %Parameter of the sigmoid transform
Prior.SigTdelta = (Prior.Exp_h+3*Prior.InformSmooth_STD_h)/2; %Parameter of the sigmoid transform
%=== Absorption Coefficient ===%
Prior.Exp_mu_a = 0.05; %Expected value of mu_a
Prior.InformSmooth_Bounds_mu_a = [0.01,0.35]; % approximately [mu_a_min, mu_a_max] 
Prior.InformSmooth_Corr_mu_a = [0.008 0.008]*RunOptions.ScalingOptical; %[correlation_x, correlation_y], Default = [8 8]
if RunOptions.ScalingOptical == 1 %When the optical forward mesh is in metres (RunOptions.ScalingOptical = 1), we convert mm^-1 to m^-1. Note that m^-1 = (10^-3mm)^-1 = 10^3mm^-1
    Prior.InformSmooth_Bounds_mu_a = Prior.InformSmooth_Bounds_mu_a*10^3;
    Prior.Exp_mu_a = Prior.Exp_mu_a*10^3;
end
Prior.AC_Var_mu_a = 0.96*10^2; %Variance of p_0
Prior.AC_Corr_mu_a = 0.002; %Correlation Length of mu_a
Prior.Exp_mu_s = 2; %Initial guess of mu_s
Prior.InformSmooth_Bounds_mu_s = [2,2.3]; % approximately [mu_s_min, mu_s_max] 
Prior.InformSmooth_Corr_mu_s = [0.008 0.008]*RunOptions.ScalingOptical; %[correlation_x, correlation_y], Default = [8 8]
if RunOptions.ScalingOptical == 1 %When the optical forward mesh is in metres (RunOptions.ScalingOptical = 1), we convert mm^-1 to m^-1. Note that m^-1 = (10^-3mm)^-1 = 10^3mm^-1
    Prior.InformSmooth_Bounds_mu_s = Prior.InformSmooth_Bounds_mu_s*10^3;
    Prior.Exp_mu_s = Prior.Exp_mu_s*10^3;
end

%% %%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%
PLOT.MeshD=1; %Data Mesh
PLOT.MeshI=1; %Inverse Mesh
PLOT.mu_aD=1; %Absorption Coefficient
PLOT.mu_sD=0; %Scattering Coefficient on Data Mesh
PLOT.mu_sI=0; %Scattering Coefficient on Inverse Mesh
PLOT.phi=1; %Light Fluence, Forward Problem Solution
PLOT.AbsorbedEnergyDensity=0; %Plot absorbed energy density
PLOT.WaveDGMForwardMesh=0; %Plot Mesh for Acoustic Forward Problem
PLOT.WaveDGMInverseMesh=0; %Plot Mesh for Acoustic Forward Problem
PLOT.WaveDGMSensorNormals=0; %Plot normal vectors on mesh for inverse inverse problem
PLOT.IntrplteAbsorbedEnergyDensity=1; %Absorbed Energy Density interpolated
PLOT.DGMForward=0; %Plot DGM Generated Forward Acoustic Data
PLOT.DGMForwardQuiver=0; %Plot DGM Generated Forward Acoustic Data using Quiver
PLOT.DGMForwardSensorData=0; %Plot DGM Generated Forward Acoustic Data
PLOT.Noise=0; %Plot noisy data
PLOT.PriorMargPoints=1; %Plot marginalisation points for prior
PLOT.PriorSamples=1; %Plot draws from prior
PLOT.PriorHist=0; %Plot prior histogram
PLOT.hRecon=1; %Absorbed Energy Density Reconstruction
PLOT.cRecon=1; %Wave Speed Reconstruction
PLOT.AdjWaveField = 0; %Plot adjoint wave field propagation
PLOT.AdjStSearchDirection = 1; %Plot Search Direction obtained using the adjoint-state method
PLOT.hGN=1; %Plot estimated absorbed energy density at each Gauss-Newton Iteration
PLOT.mu_aGN=1; %Plot estimated absorption coefficient at each Gauss-Newton Iteration
PLOT.mu_aRecon=1; %Absorption Coefficient Reconstruction
PLOT.ColourBar='EastOutside'; %For FEMPLOT of parameters
PLOT.ColourScale=[0.01, 0.03]; %For FEMPLOT of parameters
PLOT.DGMPlotBirdsEyeView = 1; %Use birds eye view for DGM plots
PLOT.DGMPlotUsezLim = 0; %Use z axis limits for DGM plots
PLOT.AbsorbedEnergyDensityzAxis = [0 300]; %z axis limits for plotting QPAT elastic wave propogation
PLOT.DGMPlotzAxis = [0 1100]; %z axis limits for plotting QPAT elastic wave propogation
PLOT.DGMAdjointPlotzAxis = [0 20]; %z axis limits for plotting QPAT elastic wave propogation
PLOT.HoldColourAxis = 1; %Hold colour axis
PLOT.ColourAxis = [0 600]; %Colour axis

%=== Pre-define Figures ===%
DefineFigures

%% =======================================================================%
%                 Generate and Plot Mesh and Parameters
%=========================================================================%
if RunOptions.GenerateMesh == 1
    %=== Using Trelis Generated Mesh ===%
    if RunOptions.UseTrelisMesh == 1
        MeshD = QPAT_TrelisMesh2DGenerator(RunOptions,'Trelis_QPATMesh_Data.inp',3,1);
        MeshI = QPAT_TrelisMesh2DGenerator(RunOptions,'Trelis_QPATMesh_Inverse.inp',3,1);
    end
    %=== Generating the Rectangular Mesh Using My Code ===%
    if RunOptions.UseMyRectangularMesh == 1
        QPAT_RectangularMesh2DGenerator
    end
else
    %=== Using Trelis Generated Mesh ===%
    RunOptions.LoadFileNameMeshD = sprintf('Mesh-%s',RunOptions.TrelisMeshDElementSize);
    load(RunOptions.LoadFileNameMeshD);
    MeshD = Mesh; 
    clear Mesh
    RunOptions.LoadFileNameMeshI = sprintf('Mesh-%s',RunOptions.TrelisMeshIElementSize);
    load(RunOptions.LoadFileNameMeshI);
    MeshI = Mesh; 
    clear Mesh
end
%keyboard %keyboard here to save newly generated mesh

%=== Plotting FEM Mesh ===%
PlotFEMMesh

%=== Generate the Parameters ===%
QPAT_EWE_DGM2D_Prmtrs 

%=== Check if Committing Inverse Crime ===%
if RunOptions.UseMyRectangularMesh == 1 && (MeshD.N_BndElmx == MeshI.N_BndElmx); 
    RunOptions.InverseCrime = 1;
else
    RunOptions.InverseCrime = 0;
end

%% =======================================================================%
%                         Display Selected Options
%=========================================================================%
if RunOptions.UseTrelisMesh == 1
    FilenamesofRunOptions  
    %=== Displayed Text ===%
    printf(['MeshD Number of Elements: ' num2str(MeshD.N_Elm)]);
    printf(['MeshD Number of Nodes: ' num2str(MeshD.N_Nodes)]);
    printf(['MeshI Number of Elements: ' num2str(MeshI.N_Elm)]);
    printf(['MeshI Number of Nodes: ' num2str(MeshI.N_Nodes)]);
    printf(['Domain Types: ' RunOptions.SaveFileNameDomain]);
    printf(['Noise Level: 0.' num2str(RunOptions.SaveFileNameNoiseLevel) '%']);
    printf(['Number of Sensors On One Boundary Edge: ' num2str(RunOptions.NumberofSensorsOnOneBoundaryEdge)]);
    printf(['Final Time: ' num2str(RunOptions.FinalTime)]);
    printf(['Number of Approximation Error Samples: ' num2str(RunOptions.AEN_Samples)]);
    printf(['\nSave File Name: ' RunOptions.SaveFileName]);
end
%keyboard %keyboard here before Diagnosing BAE
%% =======================================================================%
%                          Testing Prior Samples
%=========================================================================%
% close all
% if PLOT.PriorMargPoints == 1;
%     PLOT.Figure_PriorMargPoints = figure;
%     PLOT.Figure_Prior_Title = 'Marginalisation Points';
%     drawnow
% end
% if RunOptions.EWE_LS_ExponentialTransform == 1;
%     Prior.Exp_h = log(Prior.Exp_h);
%     Prior.InformSmooth_Bounds_h(2) = log(220);
%     Prior.InformSmooth_Bounds_h(1) = log(7);
% end
% [Prior.Exp_h,Cov,invCov,Prior.traceCov_h,~] = SmoothnessPrior_Informative(MeshI.Nodes,MeshI.Elements,2*MeshI.Dimensns,Prior.InformSmooth_Bounds_h,Prior.InformSmooth_Corr_h,Prior.Exp_h,Prior.InformSmooth_Normalize,2,PLOT);
% % [Prior.L_pr,Prior.traceCov_h,Cov_pr,~] = SmoothnessPrior_AutoCorr(MeshI.Nodes,Prior.Exp_h,Prior.AC_Var_h,Prior.AC_Corr_h,2,PLOT);
% % [~,~] = SmoothnessPrior_AutoCorr_EL(RunOptions,MeshI.Elements,MeshI.Nodes,MeshI.DomainIndices,Prior.Exp_rho_e,Prior.AC_Var_rho_e,Prior.AC_Corr_rho_e,Prior.AC_ExpShiftExp_rho_e,Prior.AC_ExpShiftVar_rho_e,Prior.Exp_cs,Prior.AC_Var_cs,Prior.AC_Corr_cs,Prior.AC_ExpShiftExp_cs,Prior.AC_ExpShiftVar_cs,1,PLOT);
% % [~,~,~,~] = SmoothnessPrior_AutoCorr(MeshI.Nodes,Prior.Exp_mu_a,Prior.AC_Var_mu_a,Prior.AC_Corr_mu_a,2,PLOT);
% % [Prior.Exp_mu_a,Cov,invCov,Prior.traceCov_mu_a,~] = SmoothnessPrior_Informative(MeshI.Nodes,MeshI.Elements,2*MeshI.Dimensns,Prior.InformSmooth_Bounds_mu_a,Prior.InformSmooth_Corr_mu_a,Prior.Exp_mu_a,Prior.InformSmooth_Normalize,2,PLOT);
% keyboard

%% =======================================================================%
%                           Forward Problem
%=========================================================================%
                          %==================%
                          %    Setting Up    %
                          %==================%                            
%%%%%%%%%%%%%%%%%%
%%% Set Up DGM %%%
%%%%%%%%%%%%%%%%%%
[DGMMeshD,PrecomputedIntrplteObjectsD] = EWE_DGM2D_Setup(RunOptions,MeshD,RunOptions.FluidMeshD,RunOptions.SolidMeshD,RunOptions.FluidDomainWithSolidLayerMeshD);
[DGMMeshI,PrecomputedIntrplteObjectsI] = EWE_DGM2D_Setup(RunOptions,MeshI,RunOptions.FluidMeshI,RunOptions.SolidMeshI,RunOptions.FluidDomainWithSolidLayerMeshI);
ConstructDGMSensors
PlotDGMMesh
AcousticForwardTimeStepsize

%=== For Diagnosing BAE ===%
% RunOptions.SaveFileNameMeshAndParameters = sprintf('MeshAndParameters-%sD-%sI-%dSensors-%sFinalTime',RunOptions.TrelisMeshDElementSize,RunOptions.TrelisMeshIElementSize,RunOptions.NumberofSensorsOnOneBoundaryEdge,RunOptions.SaveFileNameFinalTime);
% save(RunOptions.SaveFileNameMeshAndParameters);
% keyboard

                         %====================%
                         %    Computations    %
                         %====================%    
%===============================%
%    Optical Forward Problem    %
%===============================%
QPAT_DA_FEM2D_OpticalForward

%================================%
%    Acoustic Forward Problem    %
%================================%                 
QPAT_EWE_DGM2D_AcousticForward

%% =======================================================================%
%                            Inverse Problem
%=========================================================================%
                   %================================%
                   %    Acoustic Inverse Problem    %
                   %================================%
[DataVrblsOptical.hRecon,AcousticInverseItrtnInfo,Trmntn,PosteriorCovariance] = QPAT_EWE_DGM2D_AcousticInverse(RunOptions,DataVrblsWave,MeshD,MeshI,Prior,DGMMeshI,dt,PrecomputedIntrplteObjectsI,PLOT);

                %===============================================%
                %    Loading Acoustic Inverse Reconstruction    %
                %===============================================%
SaveFileNameReconstructions = sprintf('Reconstructions-%s',RunOptions.SaveFileName);
load(SaveFileNameReconstructions,'AcousticInverseItrtnInfo','Trmntn')
hRecon = AcousticInverseItrtnInfo.hRecon;
DataVrblsOptical.hRecon = hRecon(:,Trmntn);

%=== Plotting Reconstruction ===%
figure(PLOT.Figure_CurrenthRecon)
PLOT.TRI_MeshI=delaunay(MeshI.Nodes(:,1),MeshI.Nodes(:,2));
trisurf(PLOT.TRI_MeshI,MeshI.Nodes(:,1),MeshI.Nodes(:,2),full(DataVrblsOptical.hRecon));
if PLOT.DGMPlotBirdsEyeView == 1;
    view(2)
end
if PLOT.DGMPlotUsezLim == 1;
    zlim(PLOT.AbsorbedEnergyDensityzAxis);
end
shading interp %thanks Ru!
colorbar
caxis([0 200])
colormap(jet(256))
title(PLOT.Figure_CurrenthRecon_Title,'FontWeight','bold')

                    %===============================%
                    %    Optical Inverse Problem    %
                    %===============================%
if RunOptions.AE_DA == 0;
    DGMMeshI = 'Not Required';
    DataVrblsWave.SensorsI = 'Not Required';
    dt = 'Not Required';
end
[mu_aRecon,OpticalInverseItrtnInfo,Trmntn] = QPAT_DA_FEM2D_OpticalInverse(RunOptions,DataVrblsOptical,MeshD,MeshI,PrmtrsD,PrmtrsI,PrmtrsPrp,Prior,DGMMeshI,PrecomputedIntrplteObjectsI,DataVrblsWave.SensorsI,dt,PLOT);    
