function [Lsparse,bsparse]=iter_sparse(g,bsparse,para,boost,gpu,mu)

if nargin < 5 || isempty(gpu)
    gpu=cudaAvailable;
end
if nargin < 6 || isempty(mu)
    mu=1;
end
if boost==1
    gsparse=para*g;
else
    gsparse=g;
end

signd=abs(gsparse+bsparse)-1/mu;
signd(signd<0)=0;
signd=signd.*sign(gsparse+bsparse);
d=signd;
bsparse = bsparse+(gsparse-d);
Lsparse=para*(d-bsparse);
end