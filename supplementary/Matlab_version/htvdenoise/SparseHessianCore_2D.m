
function [final_g,reachCrite]=SparseHessianCore_2D(f,fidelity,hess,paral1,mu,iteration,boost,rel,fra,frames,criterion)

% if nargin < 5 || isempty(mu)
%     mu=1;
% end
% if nargin < 6 || isempty(iteration)
%     iteration=100;
% end
% if nargin < 7 || isempty(boost)
%     boost=0;
% end
% if nargin < 8 || isempty(rel)
%     rel=0.01;
% end
reachCrite = 0;
f = single(f);
gpu=0;
%%%%%
% total_time = 0; % 累计总时间
% iteration_data = zeros(4, iteration); % 预分配数据存储矩阵：4行(iter, diff_g, total_time, ttime) × iteration列
%%%%%%
% progressbar('Sparsity reconstruction');

record_max = max(f(:));
f=f./record_max;
fidelityOverMu = single(fidelity/mu);
[sx,sy,sz]=size(f);
sizeg=[sx,sy,sz];
xxfft=operation_xx(sizeg);
yyfft=operation_yy(sizeg);
xyfft=operation_xy(sizeg);

operationfft=hess*xxfft+ hess*yyfft+ 2*hess*xyfft;
normlize = single(fidelityOverMu +paral1^2 +operationfft);
clear xxfft yyfft xyfft  operationfft

bxx = zeros(sizeg,'single');
byy = bxx;
bxy = bxx;
bl1 = bxx;

PixelNum = sx*sy;
near_rel = 3*rel;
g=f;

for iter = 1:iteration
    tic;
    prev_g = g;       % 初始解g=f
    g_update = fidelityOverMu*f;
   
    [Lxx,bxx]=iter_xx(g,bxx,1,gpu);
    g_update = g_update+Lxx;
    
    [Lyy,byy]=iter_yy(g,byy,1,gpu);
    g_update = g_update+Lyy;
    
    [Lxy,bxy]=iter_xy(g,bxy,2,gpu);
    g_update = g_update+Lxy;

    [Lsparse,bl1]=iter_sparse(g,bl1,paral1,boost,gpu);

    g_update = g_update+Lsparse;

    g_update = fftn(g_update);

    g = real(ifftn(g_update./normlize));

    diff_g = norm(g(:) - prev_g(:)) / (norm(prev_g(:)) + eps);

    [OK, ZeroPercent] = ConvergenceCriterion2(g,criterion,PixelNum,record_max,1);

    if(OK)
        ttime = toc;
        disp(['breakByConvergedCriterion2 ',num2str(fra),' | ',num2str(frames) , ',  iter ' num2str(iter) ' | ' num2str(iteration),', rel_change: ', num2str(diff_g),', ZeroPercent: ', num2str(ZeroPercent),', took ' num2str(ttime) ' secs']);
        reachCrite =1;
        break;
    end

    if diff_g<near_rel
        boost=0;
    end
    if diff_g<rel
        disp(['Converged at iteration ', num2str(iter),...
                  ', ZeroPercent: ', num2str(ZeroPercent),...
              ' . Relative change: ', num2str(diff_g),...
              ' < rel: ', num2str(rel)]);
        break;
    end
    ttime = toc;
    disp([num2str(fra),' | ',num2str(frames) , ',  iter ' num2str(iter) ' | ' num2str(iteration),', rel_change: ', num2str(diff_g), ', ZeroPercent: ', num2str(ZeroPercent),', took ' num2str(ttime) ' secs']);
end

g(g<0)=0;
n_g = zeros(size(g), 'single');
n_g = ( g - min( g(:) ) ) / ( max(g(:)) - min(g(:)) );
final_g = uint16(record_max*(n_g));
clear bxx byy bzz bxz bxy byz bl1 f normlize g_update
