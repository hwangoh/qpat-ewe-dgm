function [r,s] = EWE_DGM2D_Getrs(nodes,xs,ys)
  
% EWE_DGM2D_ConstructSensorArray forms the sensor array by finding elements that 
% contains and receivers.
%
% Inputs: 
%    Nodes: 1 by 6 array containing the three x coordinates and three y coordinates of the vertices of the element containing the sensor
%    xs: x-coordinate of the sensor
%    ys: y-coordinate of the sensor
%
% Output: 
%    r: r-coordinate of the reference element such that x^k(r,s) = [xs,ys]
%    s: s-coordinate of the reference element such that x^k(r,s) = [xs,ys]
%
% Written by Timo Lahivaara

x1 = nodes(1);
x2 = nodes(2);
x3 = nodes(3);
y1 = nodes(4);
y2 = nodes(5);
y3 = nodes(6);

% r = a+b*s
a = (2*xs-x2-x3)/(-x1+x2);
b = -(-x1+x3)/(-x1+x2);

s = (2*ys -y2 -y3 -(-y1+y2)*a)/(b*(-y1+y2)-y1+y3);
r = a+b*s;
if abs(-x1+x2) > 0.0000001
    a = (2*xs-x2-x3)/(-x1+x2);
    b = -(-x1+x3)/(-x1+x2);
    
    s = (2*ys -y2 -y3 -(-y1+y2)*a)/(b*(-y1+y2)-y1+y3);
    r = a+b*s;
else
    % check this!!!
    a = (2*ys-y2-y3)/(-y1+y2);
    b = -(-y1+y3)/(-y1+y2);
    s = (2*xs -x2 -x3 -(-x1+x2)*a)/(b*(-x1+x2)-x1+x3);
    r = a+b*s;
end

if abs(-(r+s)*x1 + (1+r)*x2 + (1+s)*x3 - 2*xs) > 0.00001
    disp('!!!!! getrs !!!!!')
end
if abs(-(r+s)*y1 + (1+r)*y2 + (1+s)*y3 - 2*ys ) > 0.00001
    disp('!!!!! getrs !!!!!')
end
  
 % return
%   
%      -(r+s)*x1 + (1+r)*x2 + (1+s)*x3 - 2*xs
%      -(r+s)*y1 + (1+r)*y2 + (1+s)*y3-      2*ys
% %    keyboard
%    

% $$$ r1 = r;
% $$$ s1 = s;
% $$$ 
% $$$    syms r s
% $$$    [r,s] = solve(-(r+s)*x1 + (1+r)*x2 + (1+s)*x3 == 2*xs, ...
% $$$                  -(r+s)*y1 + (1+r)*y2 + (1+s)*y3 == 2*ys, r, s);
% $$$  
% $$$    
% $$$    r = double(r);
% $$$    s = double(s);
% $$$    
% $$$     r-r1
% $$$     s-s1
%    
% %        -(r+s)*x1 + (1+r)*x2 + (1+s)*x3
% %     2*xs
% %     -(r+s)*y1 + (1+r)*y2 + (1+s)*y3
% %      2*ys
%         -(r+s)*x1 + (1+r)*x2 + (1+s)*x3 - 2*xs
%      -(r+s)*y1 + (1+r)*y2 + (1+s)*y3-      2*ys
% 
%      keyboard
%    
