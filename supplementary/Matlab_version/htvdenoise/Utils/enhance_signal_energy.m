raw = imread('\\YangLab_NAS\Yu_Yan\Data\AausFP1\2026\ExM-Hela\Project0107004-PostExM_YY496-postExM-2_RAW_ch00.tif');
deal = imread('\\YangLab_NAS\Xianghui_Ou\Yu_Yan\Data\ExM-Hela\ExM-Hela SparseDenoiseDeconv iter10\Project0107004-PostExM_YY496-postExM-2_RAW_ch00-fidelity_150 hess_10 paral1_15 rel_0.01 cons_0 Possion_10 Pixel_0.284 NA_1.3 Lamda_488.tif');
%% 

gxx_raw = back_diff(forward_diff(raw,1,1,gpu),1,1,gpu);% back(forward)和forward(back)差异忽略，具体差异见diff_collections中注释
gyy_raw = back_diff(forward_diff(raw,1,2,gpu),1,2,gpu);
gxy_raw = forward_diff(forward_diff(raw,1,1,gpu),1,2,gpu);
hess_det_raw = gxx_raw.*gyy_raw-gxy_raw.^2;
abs_det = abs(hess_det_raw);
floatsave(abs_det,'abs_det.tif');
A = median(abs_det(:));% 50%百分位数
A2 = prctile(abs_det(:),75);%75%百分位数
%% 
gpu=0;
gxx = back_diff(forward_diff(deal,1,1,gpu),1,1,gpu);
gyy = back_diff(forward_diff(deal,1,2,gpu),1,2,gpu);
gxy = forward_diff(forward_diff(deal,1,1,gpu),1,2,gpu);
hess_det_deal = gxx.*gyy-gxy.^2;
abs_deal = abs(hess_det_deal);
%% 
record_max = max(deal(:));
A = prctile(abs_deal(:),99.5);
mask = abs_deal > A;
imwrite(uint8(255*mask),'mask.tif');

g_restored = deal;
g_restored(mask) = deal(mask) * (sum(raw(mask(:)))/ sum(deal(mask(:))));
g_restored(g_restored > record_max) = record_max; % 确保最大值约束

%%%%%%%%%结束，上面的比较好
%% 
% 不好，还是不如上面的手动99.5百分位数
abs_lambda1_norm = mat2gray(abs_deal);
level = graythresh(abs_lambda1_norm); % 自适应阈值
maska = abs_lambda1_norm > level;  
imwrite(uint8(255*maska),'maska.tif');



%% 

signal_energy = sum(deal(mask(:))); % 信号区域的当前能量
current_total = sum(deal(:)); 
target_total = sum(raw(:)); % target_total_mask = sum(raw(mask(:))); 

noise_energy = current_total - signal_energy; 
k = (target_total - noise_energy) / signal_energy;



g_restored = deal;
g_restored(mask) = deal(mask) * k; % 仅放大信号区域
g_restored(g_restored > record_max) = record_max; % 确保最大值约束

function g_restored = restore_signal_energy(raw, g, gxx, gyy, gxy)
    % 输入：
    %   raw - 原图（用于计算目标总能量）
    %   g   - 去噪后的图（需要恢复能量）
    %   gxx, gyy, gxy - 你已计算的Hessian矩阵分量（逐像素二阶差分）
    % 输出：
    %   g_restored - 能量恢复后的图
    
    % ========== 1. 基础参数与预处理 ==========
    record_max = max(g(:));  % 确保最大值不变
    [sx, sy] = size(g);
    
    % （可选但推荐）对g做轻微高斯平滑，抑制噪声对Hessian的干扰
    sigma_smooth = 1.0;  % 可根据图像分辨率调整（1~2之间）
    g_smooth = imgaussfilt(g, sigma_smooth);
    
    % ========== 2. 重新计算更鲁棒的Hessian矩阵（用平滑后的图，中心差分） ==========
    % 如果你信任自己的gxx/gyy/gxy，可跳过此步，直接用你的输入
    [gxx_smooth, gyy_smooth, gxy_smooth] = compute_hessian(g_smooth);
    
    % ========== 3. 逐像素计算Hessian的最大特征值绝对值（|λ1|） ==========
    % 2x2 Hessian矩阵 [[a, b], [b, c]] 的特征值解析解：
    % λ = [ (a+c) ± sqrt( (a-c)^2 + 4b^2 ) ] / 2
    %% 

    % a = gxx_smooth;
    % b = gxy_smooth;
    % c = gyy_smooth;
    a = gxxs;
    b = gxys;
    c = gyys;
    %% 
    
    trace_H = a + c;          % 迹（λ1 + λ2）
    det_H = a .* c - b.^2;    % 行列式（λ1 * λ2）
    sqrt_term = sqrt( (a - c).^2 + 4 * b.^2 ); % 特征值公式中的根号项
    
    lambda1 = (trace_H + sqrt_term) / 2; % 绝对值最大的特征值（通常λ1 ≥ λ2）
    lambda2 = (trace_H - sqrt_term) / 2;
    abs_lambda1 = abs(lambda1); % 用绝对值衡量曲率强度（正负仅代表明暗变化）
    %% 
    
    % ========== 4. 自适应阈值分割（Otsu法）得到信号mask ==========
    % 将abs_lambda1归一化到[0,1]，用graythresh找Otsu阈值
    abs_lambda1_norm = mat2gray(abs_lambda1);
    level = graythresh(abs_lambda1_norm); % 自适应阈值
    mask = abs_lambda1_norm > level;       % 信号点=1，噪声/背景=0
    imwrite(uint8(255*mask),'masksss.tif');
    %% 
    
    % ========== 5. 自适应能量缩放 ==========
    E_total_target = sum(raw(:));       % 目标总能量（原图的总能量）
    E_current = sum(g(:));              % 当前g的总能量
    E_signal_current = sum(g( mask(:)) ); % 当前信号区域的能量
    E_non_signal = E_current - E_signal_current; % 非信号区域能量（保持不变）
    
    % 计算信号区域需要的缩放因子
    if E_signal_current > 0
        E_signal_target = E_total_target - E_non_signal; % 信号区域需要的目标能量
        scaling_factor = E_signal_target / E_signal_current;
    else
        scaling_factor = 1; % 若无信号区域，不缩放
    end
    %% 
    
    % 生成缩放系数矩阵：信号点=scaling_factor，非信号点=1
    scalar_matrix = ones(sx, sy);
    scalar_matrix(mask) = scaling_factor;
    %% 
    
    % ========== 6. 应用缩放并截断到record_max ==========
    g_restored = g .* uint16(scalar_matrix);
    g_restored(g_restored > record_max) = record_max; % 保证最大值不变
    g_restored(g_restored < 0) = 0;                   % 保证非负（图像像素值）

end

% ========== 辅助函数：计算Hessian矩阵（中心差分，边界对称填充） ==========
function [gxx, gyy, gxy] = compute_hessian(img)
    [sx, sy] = size(img);
    img_pad = padarray(img, [1, 1], 'symmetric'); % 对称填充边界，避免边缘效应
    
    % 二阶导数xx：d²/dx² = I(x+1,y) + I(x-1,y) - 2I(x,y)
    gxx = img_pad(3:end, 2:end-1) + img_pad(1:end-2, 2:end-1) - 2 * img;
    
    % 二阶导数yy：d²/dy² = I(x,y+1) + I(x,y-1) - 2I(x,y)
    gyy = img_pad(2:end-1, 3:end) + img_pad(2:end-1, 1:end-2) - 2 * img;
    
    % 混合导数xy：d²/dxdy = [I(x+1,y+1)-I(x+1,y-1)-I(x-1,y+1)+I(x-1,y-1)] / 4
    gxy = (img_pad(3:end, 3:end) - img_pad(3:end, 1:end-2) ...
         - img_pad(1:end-2, 3:end) + img_pad(1:end-2, 1:end-2)) / 4;
end