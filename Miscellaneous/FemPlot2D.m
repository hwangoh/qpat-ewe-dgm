function FemPlot2D(g,H,sig,range,labs,cond)

%FemPlot2d Plots a function over a 2D piecewise linear basis FEM mesh.
%Originally part of SmoothnessPrior.m. I think Paul wrote this so it does
%not require a Delaunay triangulation object
%
%Inputs:
%   g - Nodal coordinates
%   H - Topology matrix
%   sig - pwl image
%   range - color range
%   labs - to put labels or not boolean
%   cond - condition for plotting
%
%Outputs:
%   Image of sig
%
% P. J. Hadwin, University of Auckland, New Zealand
%    24/07/2013 - Adapted from FemPlot for simpler meshes

%Processing Inputs
switch nargin
    case 1
        N_n=size(g,1); sig=zeros(N_n,1); range=[]; cond=''; labs=1;
        H=unique(sort(delaunay(g(:,1),g(:,2))')','rows');
    case 2
        N_n=size(g,1); sig=zeros(N_n,1); range=[]; cond=''; labs=1;
    case 3
        N_n=size(g,1); range=[]; cond=''; labs=1;
    case 4
        cond=''; labs=1;
    case 5
        cond='';
    case 6
    otherwise
        error('FemPlot2d:nargin2','Not correct number of Inputs')
end
if isempty(sig), F=zeros(N_n,1); end
if isempty(H), H=unique(sort(delaunay(g(:,1),g(:,2))')','rows'); end
labs=any(labs);

% Range
N=size(sig,2);
if isempty(range), range=[min(sig(:))' max(sig(:))']; end
if (sum(size(range))==3), range=[range(1)*ones(N,1) range(2)*ones(N,1)]; end

%Input Manipulation
if ~isempty(cond),
    ind=find(eval(cond));
    H=H(sum(ismember(H,ind),2)==3,:);
end
I1=size(H,1);
X=reshape(g(H,1),I1,3)';
Y=reshape(g(H,2),I1,3)';
C=zeros([size(X) N]);
for ii=1:N
    C(:,:,ii)=reshape(sig(H,ii),I1,3)';
end

%Plots
if N>1
    for ii=1:N
        figure
        fp=patch(X,Y,C(:,:,ii)); set(fp,'edgecolor','none');
        axis('image'),caxis(range(1,:));
%         if labs,
%             xlabel('X'),ylabel('Y'), colorbar('EastOutside')
%         else
%             axis off, box on, colorbar('off')
%         end
        drawnow
        colorbar
        caxis([0 600])
        zlim([0 600])
        colormap(jet(256))
        pause(0.03)
    end
else
    fp=patch(X,Y,C(:,:,1)); set(fp,'edgecolor','none');
    axis('image'),caxis(range(1,:));
    if labs,
        xlabel('X'),ylabel('Y'), colorbar('EastOutside')
    else
        axis off, box on, colorbar('off')
    end
    colormap(jet(256))
end
if nnz(C)==0
    set(fp,'facecolor', [1 1 1],'edgecolor','g')
    colorbar('off')
end
drawnow
hold off, grid on,