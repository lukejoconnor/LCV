function [ mu ] = weighted_mean( y,W )
%Weighted linear regression of y on X with weights W
%   Detailed explanation goes here

if ~exist('W')
    W=ones(n);
end
mu=sum(y.*W)/sum(W);


end

