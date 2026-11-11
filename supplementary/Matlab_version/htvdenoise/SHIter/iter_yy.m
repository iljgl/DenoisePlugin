function [Lyy,byy]=iter_yy(g,byy,para,hess,mu,gpu)

if nargin < 5 || isempty(mu)
    mu=1;
end
if nargin < 6 || isempty(gpu)
    gpu=0;
end

gyy = back_diff(forward_diff(g,1,2,gpu),1,2,gpu);
signdyy = abs(gyy+byy)-hess/mu;
signdyy(signdyy<0)=0;
signdyy=signdyy.*sign(gyy+byy);
dyy=signdyy;
byy = byy+(gyy-dyy);
Lyy=para*back_diff(forward_diff(dyy-byy,1,2,gpu),1,2,gpu);