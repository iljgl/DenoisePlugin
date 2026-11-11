function [Lxx,bxx]=iter_xx(g,bxx,para,hess,mu,gpu)

if nargin < 5 || isempty(mu)
    mu=1;
end
if nargin < 6 || isempty(gpu)
    gpu=0;
end
gxx = back_diff(forward_diff(g,1,1,gpu),1,1,gpu);
signdxx = abs(gxx+bxx)-hess/mu;
signdxx(signdxx<0)=0;
signdxx=signdxx.*sign(gxx+bxx);
dxx=signdxx;
bxx = bxx+(gxx-dxx);
Lxx = para*back_diff(forward_diff(dxx-bxx,1,1,gpu),1,1,gpu);