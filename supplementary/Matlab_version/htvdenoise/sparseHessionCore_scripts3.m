clear
clc
addpath('./Utils/');
addpath('./SHOperation/');
addpath('./SHIter/');
addpath('./IterativeDeblur/');
%% 
filename = 'EzrA300_img_000000029_Channel2_000.tif';
tif_stack = imreadTiff(filename);
%% 
% tif_stack = tif_stack(500:1000,500:1000,:);
% tif_stack_part = tif_stack;
filename_prefix = filename(1:end-4);%这里暂时移除.tif后缀，增添一些关键信息后再补上;
fidelity=150;
hess=10;
paral1=15;
mu=1;
iteration=100;
boost=1;
rel=0.001;
bit=16;
% tt = imread('\\YangLab_NAS\XiangHui_Ou\HJ_YD_SparseDeconv\YD\MreB motion-20250221 data\100 nM SiftReg\561nm-30mW-OD1-100ms-0000 SiftReg 100in200.tif');
% ttg = SparseHessian_core_2D_optimized(tt,fidelity,hess,paral1,mu,iteration,boost,rel);
% ttg = SparseHessian_core_2D_optimized(ttg,fidelity,hess,paral1,mu,iteration,boost,rel);
[~,~,frames] = size(tif_stack);
% PixelNum=h*w;
ZeroPercentCriterion = 0.995;
% B(B > max_A) = max_A
% max_100 = max(max(tif_stack(:,:,100)));
% bkg=75;
ReachCrite = zeros(1,frames);
% parpool(6);
for fra=1:frames
    
    % max_cur = max(max(tif_stack_part(:,:,fra)));
    % paral12 = single(max_100/max_cur *paral1);
    tmp = tif_stack(:,:,fra);
    % tmp(tmp>minMax)=minMax;
    % tmp = tmp-bkg;
    % tmp(tmp<0)=0;

    [tmp,ReachCrite(fra)]=SparseHessianCore_2D(tmp,fidelity,hess,paral1,mu,iteration,boost,rel,fra,frames,ZeroPercentCriterion);
    while(~ReachCrite(fra))
        [tmp,ReachCrite(fra)]=SparseHessianCore_2D(tmp,fidelity,hess,paral1,mu,iteration,0,rel,fra,frames,ZeroPercentCriterion);  
    end
    tif_stack(:,:,fra)=tmp;
    % while(~Converged2(tmp,criterion,PixelNum,0,0)) % 只获取第一个返回值
    %     tmp = SparseHessian_core_2D_optimized(tmp,fidelity,hess,paral1,mu,iteration,boost,rel,fra,frames,criterion);
    % end
        
    % tif_stack(:,:,fra)=tmp;
    progressbar(fra/frames);
end

fullfilename1 = [filename_prefix,'-fidelity_', num2str(fidelity) ,' hess_',num2str(hess),' paral1_',num2str(paral1),' boost_',num2str(boost),' rel_',num2str(rel),'  zeroPer_',num2str(ZeroPercentCriterion),' 1epochs.tif'];
imwritestack(tif_stack,fullfilename1,bit) ;%(stack, filename,bit)
