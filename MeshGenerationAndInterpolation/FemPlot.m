function C = FemPlot(Mesh,Prmtr,PLOT)

% FemPlot Plots 2D functions with a piecewise linear basis when 
% given a FEM mesh. It also reorganizes the Prmtr vector into an elementwise
% representation of its values.
%
% Inputs:
%   Mesh - FEM Mesh Structure
%        Nodes - Coordinates of Nodes [x y z] needs to be at least 2-D
%        Elements - Nodes of Elements, assumed to be linear Basis
%   Prmtr - 1 by N_Nodes vector where each entry reperesents the
%       coefficient's value at that node
%   PLOT - Option Structure (optional)
%        colorbar - string defining location of the Colorbar, supported 
%                   locations can be found in the help of colorbar. Must be 
%                   same location for all plots.
%        Range - N by 2 matrix of Colorbar limits. Can be just 2 by 1, i.e.
%                same for plots
%   DEFAULTS: 'EastOutside', [min(F)' max(F)']
%
% Outputs:
%   C = 3 by N_elm matrix where each column represents an element and each
%       entry of the column represents the value of the parameter at the
%       node
%
% Hwan Goh, University of Auckland, New Zealand 15/06/2013
% Adapted from P. J. Hadwin, University of Auckland, New Zealand
%    24/11/11 - Original

[X Y C cbar N Range]=CheckInput(Mesh,Prmtr,PLOT);

%Plots
hold on
for ii=1:N
    view(2),
    fp=patch(X,Y,C(:,:,ii)); set(fp,'edgecolor','none');
    axis('image'),caxis(Range(ii,:)),colorbar(cbar),colormap(jet(256)),drawnow
end


%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%

function [X Y C cbar N Range]=CheckInput(Mesh,Prmtr,PLOT)

%CheckInput Checks and manipulates the various inputs for FemPlot
%
%Inputs: 
%   Same as those above
%
%Outputs:
%   X,Y - x and y co-ordinates for the scattering coefficient prmtr
%   C - [size(X) N] matrix of colour levels for plots
%   N - number of functions to be plotted
%   cbar - colorbar location
%   R - N by 2 matrix of ranges for the colorbar
%
%Hwan Goh, University of Auckland 15/06/2013 
%Adapted from P. J. Hadwin, University of Auckland, New Zealand, 24/11/2011

%Colourbar
cbar = PLOT.ColourBar;
if ~isfield(PLOT,'Range');
    Range=[0.98*min(Prmtr)' 1.02*max(Prmtr)'];
else
    Range = PLOT.Range;
end

%Input Manipulation 
N=size(Prmtr,2);
g=Mesh.Nodes; 
H=Mesh.Elements;
[I1 I2]=size(H);

X=reshape(g(H(:,1:I2),1),I1,I2)';
Y=reshape(g(H(:,1:I2),2),I1,I2)';
C=zeros([size(X) N]);

for ii=1:N
    C(:,:,ii)=reshape(Prmtr(H(:,1:I2),ii),I1,I2)';
end