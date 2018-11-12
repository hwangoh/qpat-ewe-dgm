% QPAT_EWE_DGM2D_Prmtrs defines the target parameters
%
% Hwan Goh, University of Auckland, New Zealand - 6/7/2015

%=========================================================================%
%                           Target Generation
%=========================================================================%
% The following two statements motivate this code:
%       "...most common is to use the diffustion and absorption coefficeints as image parameters", 
%       "The coefficients (kappa, mu_a) are in most cases assumed to be mutally independent"
%            - 'Novel approaches to image reconstruction in diffusion tomography', Ville Kolehmainen, Ph.D. thesis

PrmtrsPrp.g = 0.8; %mean of the cosine of the scattering angle - 'Reconstructing absorption and scattering...,' T.Tarvainen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Background Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PrmtrsPrp.bckgrnd_mu_a = 0.01; %ambient absorption coefficient [mm^-1] - 'Reconstructing absorption and scattering...,' T.Tarvainen
PrmtrsPrp.bckgrnd_mu_s = 2; %ambient scattering coefficient [mm^-1] - 'Reconstructing absorption and scattering...,' T.Tarvainen
PrmtrsPrp.bckgrnd_mu_rs = (1-PrmtrsPrp.g)*PrmtrsPrp.bckgrnd_mu_s; %ambient reduced scattering coefficient [mm^-1]
PrmtrsPrp.bckgrnd_kappa = 1/(2*(PrmtrsPrp.bckgrnd_mu_a + PrmtrsPrp.bckgrnd_mu_rs)); %ambient diffusion coefficient [mm^2/s]

%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter Range %%%
%%%%%%%%%%%%%%%%%%%%%%%
PrmtrsPrp.scalar_mu_a = 0.5; %control the range of the absorption coefficient
PrmtrsPrp.scalar_mu_s = 0.5; %control the range of the scattering coefficient
PrmtrsPrp.scalar_kappa = 0.5; %control the range of the diffusion coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters Maximums %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
PrmtrsPrp.max_mu_a = 0.1; %maximum absorption coefficient [mm^-1] - 'Reconstructing absorption and scattering...,' T.Tarvainen
PrmtrsPrp.max_mu_s = 10; %maximum scattering coefficient [mm^-1] - 'Reconstructing absorption and scattering...,' T.Tarvainen
PrmtrsPrp.max_mu_rs = (1-PrmtrsPrp.g)*PrmtrsPrp.max_mu_s; %maximum reduced scattering coefficient [mm^-1]
PrmtrsPrp.max_kappa = 1/(2*(PrmtrsPrp.max_mu_a + PrmtrsPrp.max_mu_rs)); %maximum diffusion coefficient [mm^2/s]

                      %=============================%
                      %    Generating Parameters    %
                      %=============================%

%%%%%%%%%%%%%%%%%%%%%%%
%%% Generating mu_a %%%
%%%%%%%%%%%%%%%%%%%%%%% 
[PrmtrsD.mu_a,~,PrmtrsPrp.R_mu_a] = GenPrmtrsI(RunOptions,MeshD,MeshI,PrmtrsPrp.bckgrnd_mu_a,PrmtrsPrp.scalar_mu_a,'absorption'); %The absorption coefficient
%Plotting mu_a
PrmtrsD.mu_a_elmts = FEM_Construct3ByN_ElmArray(MeshD.Nodes,MeshD.Elements,PrmtrsD.mu_a); % Creates a 3 by N_Elm matrix representing the parameter values on each element. 
if PLOT.mu_aD == 1;
    figure(PLOT.Figure_mu_aD)
    PLOT.TRI_MeshD=delaunay(MeshD.Nodes(:,1),MeshD.Nodes(:,2));
    trisurf(PLOT.TRI_MeshD,MeshD.Nodes(:,1),MeshD.Nodes(:,2),PrmtrsD.mu_a);
    view(2)
    shading interp %thanks Ru!
    colorbar
    caxis([0,0.25])
    colormap(jet(256))
    axis 'image'
    title(PLOT.Figure_mu_aD_Title,'FontWeight','bold')
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Generating mu_s %%%
%%%%%%%%%%%%%%%%%%%%%%%
[PrmtrsD.mu_s,PrmtrsI.mu_s,PrmtrsPrp.R_kappa] = GenPrmtrsI(RunOptions,MeshD,MeshI,PrmtrsPrp.bckgrnd_mu_s,PrmtrsPrp.scalar_mu_s,'scattering'); %The diffusion coefficient
%plotting mu_s
PrmtrsD.mu_s_elmts = FEM_Construct3ByN_ElmArray(MeshD.Nodes,MeshD.Elements,PrmtrsD.mu_s); % Creates a 3 by N_Elm matrix representing the parameter values on each element. 
if PLOT.mu_sD == 1;
    figure(PLOT.Figure_mu_sD)
    trisurf(PLOT.TRI_MeshD,MeshD.Nodes(:,1),MeshD.Nodes(:,2),PrmtrsD.mu_s);
    view(2)
    shading interp %thanks Ru!
    colorbar
    colormap(jet(256))
    axis 'image'
    title(PLOT.Figure_mu_sD_Title,'FontWeight','bold')
end
PrmtrsI.mu_s_elmts = FEM_Construct3ByN_ElmArray(MeshI.Nodes,MeshI.Elements,PrmtrsI.mu_s); % Creates a 3 by N_Elm matrix representing the parameter values on each element. 
if PLOT.mu_sI == 1;     
    figure(PLOT.Figure_mu_sI)
    PLOT.TRI_MeshI=delaunay(MeshI.Nodes(:,1),MeshI.Nodes(:,2));
    trisurf(PLOT.TRI_MeshI,MeshI.Nodes(:,1),MeshI.Nodes(:,2),PrmtrsI.mu_s);
    view(2)
    shading interp %thanks Ru!
    colorbar
    colormap(jet(256))
    axis 'image'
    title(PLOT.Figure_mu_sI_Title,'FontWeight','bold')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scaling Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%When the optical forward mesh is in metres (RunOptions.ScalingOptical = 1), we convert mm^-1 to m^-1. Note that m^-1 = (10^-3mm)^-1 = 10^3mm^-1
if RunOptions.ScalingOptical == 1
    PrmtrsD.mu_a = PrmtrsD.mu_a*10^3;
    PrmtrsD.mu_a_elmts = PrmtrsD.mu_a_elmts*10^3;
    PrmtrsD.mu_s = PrmtrsD.mu_s*10^3;
    PrmtrsD.mu_s_elmts = PrmtrsD.mu_s_elmts*10^3;%     PrmtrsI.mu_s = PrmtrsI.mu_s*10^3; %for testing with true mu_s
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generating kappa %%%
%%%%%%%%%%%%%%%%%%%%%%%%                    
PrmtrsD.kappa = zeros(MeshD.N_Nodes,1);
for i=1:MeshD.N_Nodes;
    PrmtrsD.kappa(i) = 1/(2*(PrmtrsD.mu_a(i) + ((1-PrmtrsPrp.g)*PrmtrsD.mu_s(i))));
end
PrmtrsD.kappa_elmts = zeros(3,MeshD.N_Elm);
PrmtrsD.kappa_elmts = reshapeNodes(MeshD,PrmtrsD.kappa);


%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set "known" mu_s %%%
%%%%%%%%%%%%%%%%%%%%%%%%
PrmtrsI.mu_s = PrmtrsI.mu_sEst*ones(MeshI.N_Nodes,1);
PrmtrsI.mu_s_elmts = PrmtrsI.mu_sEst*ones(3,MeshI.N_Nodes);
if RunOptions.ScalingOptical == 1
    PrmtrsI.mu_s = PrmtrsI.mu_s*10^3;
    PrmtrsI.mu_s_elmts = PrmtrsI.mu_s_elmts*10^3;
end
