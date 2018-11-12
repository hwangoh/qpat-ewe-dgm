% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; NODETOL = 1e-12;
% Compute nodal set
[x,y] = Nodes2D(N); [r,s] = xytors(x,y);

% Build reference element matrices
V = Vandermonde2D(N,r,s); invV = inv(V);
MassMatrix = invV'*invV;
[Dr,Ds] = Dmatrices2D(N, r, s, V);

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

% find all the nodes that lie on each edge. Note: the extra two bits of
% code on each line was used by Timo in "BuildPNonCon2D" to reorganize the
% indexing of the element boundary nodes
fmask1   = find( abs(s+1) < NODETOL)'; [~,fmask] = sort(r(fmask1), 'ascend');  fmask1 = fmask1(fmask);
fmask2   = find( abs(r+s) < NODETOL)'; [~,fmask] = sort(r(fmask2), 'descend');  fmask2 = fmask2(fmask);
fmask3   = find( abs(r+1) < NODETOL)'; [~,fmask] = sort(s(fmask3), 'descend');  fmask3 = fmask3(fmask);
Fmask  = [fmask1;fmask2;fmask3]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);

% Create surface integral terms
LIFT = Lift2D();

% calculate geometric factors
[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);

% calculate geometric factors
[nx, ny, sJ] = Normals2D();
Fscale = sJ./(J(Fmask,:));

% Build connectivity matrix
%[EToE, EToF] = tiConnect2D(EToV);
[EToE, EToF] = Connect2D(EToV);

% Build connectivity maps
BuildMaps2D();

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs] = GradVandermonde2D(N, r, s);
Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');

% Return nodes outside of mesh to within the mesh
if exist('kgrid','var') == 1 || exist('kgridD','var') == 1 || exist('kgridI','var') == 1;
    [x,y] = ReturnNodesWithinMesh(x,y,Mesh);
end