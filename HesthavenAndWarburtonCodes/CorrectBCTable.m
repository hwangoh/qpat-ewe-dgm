function [BCType,TopBndInfo,LeftBndInfo,RightBndInfo,BttmBndInfo] = CorrectBCTable(Nfaces,K,VX,VY,EToV,BCType,topBndVerts,leftBndVerts,rightBndVerts,bttmBndVerts,TopBC,LeftBC,RightBC,BttmBC)

% Purpose: Setup BCType for boundary conditions in 2D
% By Allan P. Engsig-Karup
%
% Edited to also find top, left, right, bottom boundary indices

VNUM = [1 2;2 3;3 1]; % face orientations

TopBndInfo = zeros(K,3);
LeftBndInfo = zeros(K,3);
RightBndInfo = zeros(K,3);
BttmBndInfo = zeros(K,3);

for k = 1:K    
    % Test for each edge
    for l = 1:Nfaces 
        m = EToV(k,VNUM(l,1)); n = EToV(k,VNUM(l,2));    
        % if both points are on the boundary then it is a boundary point!
        %top edge
        ok=sum(ismember([m n],topBndVerts));
        if ok==2
            if (VX(m) == VX(n)) || (VY(m) == VY(n)) %This is to avoid the possibility that the face has one node on a different boundary to its other node (i.e. diagonal face from left domain boundary to bottom domain boundary)
            BCType(k,l)=TopBC;
            TopBndInfo(k,l) = 1;
            end
        end
        %left edge
        ok=sum(ismember([m n],leftBndVerts));
        if ok==2
            if (VX(m) == VX(n)) || (VY(m) == VY(n)) %This is to avoid the possibility that the face has one node on a different boundary to its other node (i.e. diagonal face from left domain boundary to bottom domain boundary)
            BCType(k,l)=LeftBC;
            LeftBndInfo(k,l) = 2;
            end
        end
        %right edge
        ok=sum(ismember([m n],rightBndVerts));
        if ok==2
            if (VX(m) == VX(n)) || (VY(m) == VY(n)) %This is to avoid the possibility that the face has one node on a different boundary to its other node (i.e. diagonal face from left domain boundary to bottom domain boundary)
                BCType(k,l)=RightBC;
                RightBndInfo(k,l) = 3;
            end
        end
        %bottom edge
        ok=sum(ismember([m n],bttmBndVerts));
        if ok==2
            if (VX(m) == VX(n)) || (VY(m) == VY(n)) %This is to avoid the possibility that the face has one node on a different boundary to its other node (i.e. diagonal face from left domain boundary to bottom domain boundary)
                BCType(k,l)=BttmBC;
                BttmBndInfo(k,l) = 4;
            end
        end
    end
end

return

%==============================Debugging==================================%
%=== Plotting ===% %Copy to line 11
% figure
% triplot(EToV,VX,VY);
% hold on
% for i=1:size(mapnodes)
%     mapnodes(i)
%     [VX(mapnodes(i)),VY(mapnodes(i))]
%     plot(VX(mapnodes(i)),VY(mapnodes(i)),'ro');
%     hold on
%     pause(0.01)  
% end

%=== Check ===% %Copy to line 24 (after Plotting has also been copied)
% Check = zeros(3*K,10);
% h = 1;
% BndFace = 0;

% Check(h,1) = k; %Copy to line 32 (after Plotting has also been copied)
% Check(h,2) = m;
% Check(h,3) = n;
% Check(h,4:5) = ismember([m n],mapnodes);
% Check(h,6) = sum(ismember([m n],mapnodes));
% Check(h,7:8) = [VX(m),VY(m)];
% Check(h,9:10) = [VX(n),VY(n)];
% if Check(h,6) == 2
%     BndFace = BndFace+1;
% end
% h = h+1;
% figure
% triplot(EToV,VX,VY);
% hold on
% plot(VX(m),VY(m),'ro');
% plot(VX(n),VY(n),'ro');
% keyboard
% close all

% Diff = BCType - BCType1; %Copy to line 56
% [r,~] = find(Diff~=0);
%=========================================================================%