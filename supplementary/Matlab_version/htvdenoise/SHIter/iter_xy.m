function [Lxy,bxy]=iter_xy(g,bxy,para,hess,mu,gpu)


if nargin < 5 || isempty(mu)
    mu=1;
end
if nargin < 6 || isempty(gpu)
    gpu=0;
end
gxy = forward_diff(forward_diff(g,1,1,gpu),1,2,gpu);
signdxy = abs(gxy+bxy)-hess/mu;
signdxy(signdxy<0)=0;
signdxy=signdxy.*sign(gxy+bxy);
dxy=signdxy;
bxy = bxy+(gxy-dxy);
Lxy=para* back_diff(back_diff(dxy-bxy,1,2,gpu),1,1,gpu);