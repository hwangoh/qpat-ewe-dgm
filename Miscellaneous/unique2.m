function [Y P Q]=unique2(X,S)

%[Y P Q]=unique2(X,S);
%
% unique2 Is unique but has an unsort option and option to find unique rows
%          in a matrix. The first occurance of each number (row) is kept. 
%      NOTE: If sort is used unique2 does the same as unique(X,'rows','first').
%
%Inputs:
%   X - Vector or matrix of which the unique elements are to be found 
%   S - Vector of Booleans.
%         S(1) = whether to sort or not. Default S=0
%         S(2) = whether to use matrix sort. Default S=1
%
%Outputs:
%   Y - Column Vector (Matrix) of unique elements (rows)
%   P - Indices of the first Unique elements (rows) in X so that Y=X(P,:)
%        Note if S(1)=1 P is I as seen in documentation for unique
%   Q - If S(1)=0 Indices of non-unique elements (row rows if matrix)
%       If S(1)=1 J as seen in documentation for unique
%
% Hwan Goh, University of Auckland, New Zealand - 31/05/2013
% Code adapted from P. J. Hadwin, University of Auckland, New Zealand
%    27/09/2011 - Original
%    18/01/2012 - Matrices Added


%%Setting Defaults and input checks
if nargin<2, S=[0 1];  end
if ~any(find([0;1]==S(1)))||~(length(S)==2)||~any(find([0;1]==S(2)))
    error('unique2:Scheck','Provided S is unsupported in unique2.')
end
if ~(ndims(X)==2), 
    error('unique2:dimCheck','Provided X must be a vector or 2D matrix.')
end

%%Computation
if ~S(2)||(size(X,1)==1)
    X=X(:);
end

N=size(X,1);
[Y, P, Q]=unique(X,'rows','first');

if ~S(1)
    P=sort(P(:));
    Q=setdiff(1:N,P)';
    Y=X(P,:);
end
