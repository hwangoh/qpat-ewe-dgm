
function [EToV, VX, VY, Nel, Nv, BCType, Bnd_NodeInd, doms] = ...
    ReadGrid_trelis(filename, Ndoms, Nbc, visualization)

% function reads the computational grid exported from Trelis software
% Inputs: filename (meshfile), Ndoms (number of domains), Nbc (number of
% boundary conditions)
% modified: Timo L. Jan. 28, 2016
% Edited: Hwan Goh, Jan 30 2018, allow output of list of boundary nodes (lines 47 to 56)

% opens the meshfile
fid = fopen(filename,'r');

%%% read header %%%
for ii = 1:9
    str = fgets(fid);
end

% read vertices
disp('  reading vertex data')
[rawdata,cnt] = fscanf(fid,'%d, %f, %f', inf);
rawdata = reshape(rawdata, 3, cnt/3);
VX = rawdata(2, :);
VY = rawdata(3, :);
clear rawdata;
Nv = length(VX); % number of vertices

% skip info lines
for ii = 1:3
    str = fgets(fid);
end

% read element data (topology)
EToV = []; doms = [];
for ii = 1:Ndoms
    disp(['  reading element data ' num2str(ii)])
    [rawdata,cnt] = fscanf(fid,'%d, %d, %d, %d', inf);
    Htmp = reshape(rawdata, 4, cnt/4)';
    EToV = [EToV; Htmp(:,2:end)];
    doms = [doms; zeros(length(Htmp(:,1)),1)+ii];
    str = fgets(fid);
    disp(['  Nelements: ' num2str(length(Htmp(:,1)))])    
end
clear rawdata;
Nel = length(EToV(:,1)); % number of elements

% skip info lines
for ii = 1:2
    str = fgets(fid);
end

% read node data - I added this
disp(['  reading nodeset data ' num2str(ii)])
ingroup= sscanf(str, '*NSET, NSET=NS%d');
[Bnd_NodeInd,~] = fscanf(fid,'%d,', inf);
str = fgets(fid); %Need this to move on to next section? I saw it in line 41

% skip info lines
for ii = 1:2
    str = fgets(fid);
end

% read boundary condition info
bc = [];
for ii = 1:Nbc
    ts = 0;
    disp(['  reading bc data ' num2str(ii)])
    ingroup= sscanf(str, '*ELSET, ELSET=SS%d_E%d');
    [A,cnt] = fscanf(fid,'%d,', inf);
    bc = [bc; [A(:) zeros(cnt,1)+ingroup(2) zeros(cnt,1)+ii]];
    ts = ts+1;
    str = fgets(fid);
    str = fgets(fid);
    while strcmp(str(2), 'S') ~= 1
        ingroup= sscanf(str, '*ELSET, ELSET=SS%d_E%d');
        [A,cnt] = fscanf(fid,'%d,', inf);
        bc = [bc; [A(:) zeros(cnt,1)+ingroup(2) zeros(cnt,1)+ii]];
        ts = ts+1;
        str = fgets(fid);
        str = fgets(fid);
    end
    for ii = 1:ts+1
        str = fgets(fid);
    end    
end
% close meshfile
fclose(fid);

% at point point, we know elements (+ faces) which are on the external
% boundary  
% now we can compute the corresponding vertices to these boundary faces
or = [1 2;2 3;1 3]; % this are node points that correspond to face index 
                    %--> face index 1 is defined with vertices with indices 1 2, etc
tmp = zeros(length(bc(:,1)), 4);
% tmp: [elementindex, vertex1, vertex2, boundarycondition]
for ii = 1:length(bc(:,1))
    elm = bc(ii,1);
    tmp(ii,:) = [elm EToV(elm, or(bc(ii,2), :)) bc(ii,3)];
end

% visualize boundary conditions
if visualization.initialization
    markertypes = {'*','o','x','d'};
    figure    
    set(gca, 'fontsize', 15)
    [~, loopids] = unique(bc(:,3));
    for ii = 1:length(loopids)%1:length(bc(:,1))
        randmarker = randi(4,1,1);
       myMap = rand(1, 3); % randomize color
       rows = tmp(:,end) == bc(loopids(ii), 3); % collect rows that contain corresponding bc vertices
       xbc = VX(tmp(rows,[2 3])); % x-coordinates
       ybc = VY(tmp(rows,[2 3])); % y-coordinates
       plot(xbc(:), ybc(:), char(markertypes(randmarker)), 'Color', myMap)
       hold on
    end
    axis equal tight 
    tmpval = axis;
    scax = 0.02*tmpval(2);
    scay = 0.02*tmpval(4);
    axis([tmpval(1)-scax tmpval(2)+scax tmpval(3)-scay tmpval(4)+scay])
    grid on
    drawnow
    print -dpng -r300 boundaryconditions.png
    close all
end

% next we check that grid has faces oriented anticlockwise
disp(' check H')
for ii = 1:Nel
    iv = EToV(ii,:); % node indices for ii'th element
    x1 = VX(iv); % x-vertices
    y1 = VY(iv); % y-vertices
    rpx = mean(x1); % x-coordinate (center point)
    rpy = mean(y1); % y-coordinate (center point)
    angles = zeros(3,1); % angle between two vectors 
    for jj = 1:3 % loop over vertices
        angles(jj) = atan2((y1(jj)-rpy),(x1(jj)-rpx));%/pi*180;
    end
    [~, id] = sort(angles);
    EToV(ii,:) = EToV(ii,id); % sort 
end

% BCType matrix defines the boundary condition for each face of each
% element -> if value is 0, face is on the internal boundary and ~= 0 on
% the external boundary
BCType = 0*EToV;
for ii = 1:length(bc(:,1))
    elm = tmp(ii,1); %element
    vec = ismember(EToV(elm,:), tmp(ii,2:3));%vertex id
    if sum(vec(1:2)) == 2%find correct face
        val = 1;
    elseif sum(vec([2 3])) == 2
        val = 2;
    else
        val = 3;
    end
    BCType(tmp(ii,1), val) = tmp(ii,4);%set bc
end
return
