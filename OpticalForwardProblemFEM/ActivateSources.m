function [ActiveSources_Nodes] = ActivateSources(RunOptions,Mesh)

% ActivateSources controls which light sources are involved in illuminating
% for a given illumination pattern.
%
% Inputs:
%    Mesh:
%       Mesh.Elements - Matrix where row number corresponds to the element
%                        number and the entries of the row are the vertices of the element
%       Mesh.N_Nodes - Number of nodes
%       Mesh.N_Lght - Number of light sources. For rectangular mesh this is
%                     Mesh.N_Lght \in {1,2,3,4} that represents how many domain are light sources. 
%                     This is no longer relevant for my PhD codes since I removed the option of different illumination
%                     patterns. Doesn't take much to re-implement this though.
%
% Outputs:
%    ActiveSources_Nodes - 1 by N_Nodes vector containing a 1 when a node
%                          corresponds with a node of an element under an
%                          active light source.
%
% Hwan Goh, 18/08/2013, University of Auckland, New Zealand
% Last Edited 25/01/2019: Now supports circular meshes 

%% Shortening Labels
Elements = Mesh.Elements;
N_Nodes = Mesh.N_Nodes;
Lght_Nelm = Mesh.Lght_Nelm;
Lght_ElmtInd = Mesh.Lght_ElmtInd;

ActiveSources_Nodes = zeros(1,N_Nodes);

%% Activating Sources
SourceNodes = zeros(Lght_Nelm,3); %Lght_Nelm by 3 matrix where each row contains the nodes of the elements under the active sources; filled in by the for loop below.
for k=1:Lght_Nelm;
    SourceNodes(k,:) = Elements(Lght_ElmtInd(k),:);
end
nonzeroNodes = find(SourceNodes);
for l=1:size(nonzeroNodes,1);
    ActiveSources_Nodes(SourceNodes(nonzeroNodes(l))) = 1;
end
