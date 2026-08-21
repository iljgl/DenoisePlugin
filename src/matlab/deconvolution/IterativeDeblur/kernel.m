function psf=kernel(pixel, lambda, NA,z,nn,r)

if nargin < 4 || isempty(z)
    z=0;
end
if nargin < 5 || isempty(nn)
    nn=64;
end
if nargin < 6 || isempty(r)
    r=2;
end

psf=Generate_PSF(pixel,lambda,nn,NA,z,r);
psf=psf./sum(sum(psf));
