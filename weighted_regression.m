function [ beta ] = weighted_regression( X,y,W )
%Weighted linear regression of y on X with weights W
%   Detailed explanation goes here

[n,p]=size(X);
if ~exist('W')
    W=ones(n,1);
end
beta=(X'*diag(sparse(W))*X)\(X'*diag(sparse(W))*y);


end

