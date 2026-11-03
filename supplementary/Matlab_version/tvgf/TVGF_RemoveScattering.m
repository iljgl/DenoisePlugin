function [U, S, a_map, b_map] = TVGF_RemoveScattering(G, S0, WinSize, lambda)
    % 输入：
    %   G 去噪后的引导图,单通道，double，[0,1]
    %   S0 四分之一窗口暗通道初始估计,double，[0,1]
    %   WinRowSize 标量，default:15
    %   lambda 正则系数，default:0.01）
    % 输出：
    %   U - 最终去散射光的清晰图像
    %   S - 变分模型最终的散射光图像
    %   a_map - 逐像素的a系数图
    %   b_map - 逐像素的b系数图

    [Gy, Gx] = size(G);
    halfWin = floor(WinSize / 2);
    
    % N_win = WinRow * WinCol; % 窗口像素总数
    % N_win = WinSize * WinSize;
    % 1：预计算G的梯度1范数 |∇G|1 
    % % 中心差分
    kernelx = zeros(Gy, Gx);
    kernelx(1,end)=-0.5;
    kernelx(1,2)=0.5;
    kernely = zeros(Gy, Gx);
    kernely(2,1) = 0.5;
    kernely(end,1) = -0.5;
    otfx = fftn(kernelx,size(G));
    otfy = fftn(kernely,size(G));
    
    G_x = real(ifftn(fftn(G) .* otfx));
    G_y = real(ifftn(fftn(G).* otfy));

    abs_grad_G = abs(G_x) + abs(G_y);
    % for y = 1:Gy
    %     for x = 1:Gx
    %         % x方向梯度（前向差分，边界用0）
    %         if x < Gx
    %             dx = G(y, x+1) - G(y, x);
    %         else
    %             dx = 0;
    %         end
    %         % y方向梯度（前向差分，边界用0）
    %         if y < Gy
    %             dy = G(y+1, x) - G(y, x);
    %         else
    %             dy = 0;
    %         end
    %         abs_grad_G(y, x) = abs(dx) + abs(dy);
    %     end
    % end
    % abs_grad_G = zeros(Gy, Gx);
    % abs_grad_G(y, x) = abs(G_x) + abs(G_y);

    % 2：初始化输出图
    S = zeros(Gy, Gx);
    a_map = zeros(Gy, Gx);
    b_map = zeros(Gy, Gx);

    % 3：滑动窗口遍历每个像素
    for y = 1:Gy
        for x = 1:Gx
            r_start = max(1, y - halfWin);
            r_end = min(Gy, y + halfWin);
            c_start = max(1, x - halfWin);
            c_end = min(Gx, x + halfWin);

            % Extract window data
            win_G = G(r_start:r_end, c_start:c_end);
            win_S0 = S0(r_start:r_end, c_start:c_end);
            win_grad = abs_grad_G(r_start:r_end, c_start:c_end);
            win_G_vec = win_G(:);
            win_S0_vec = win_S0(:);
            win_grad_vec = win_grad(:);
            current_N = numel(win_G);
            % 5个统计量
            mean_G = mean(win_G_vec);
            mean_S0 = mean(win_S0_vec);
            mean_grad = mean(win_grad_vec);
            % 中心化
            cent_G = win_G_vec - mean_G;
            cent_S0 = win_S0_vec - mean_S0;
            var_G = sum(cent_G.^2) / current_N;
            cov_GS = sum(cent_G .* cent_S0) / current_N;
            % 软阈值T
            T = lambda * mean_grad;
            % 计算a_xy
            if var_G < 1e-8 % 数值稳定性：窗口内G完全平坦
                a_xy = 0;
            else
                if cov_GS > T
                    a_xy = (cov_GS - T) / var_G;
                elseif cov_GS < -T
                    a_xy = (cov_GS + T) / var_G;
                else
                    a_xy = 0;
                end
            end
            % 计算b_xy
            b_xy = mean_S0 - a_xy * mean_G;
            % 计算当前像素的S_xy
            S_xy = a_xy * G(y, x) + b_xy;
            S_xy = max(0, min(1, S_xy)); % 值域截断
            S(y, x) = S_xy;
            a_map(y, x) = a_xy;
            b_map(y, x) = b_xy;
        end
    end

    % 4：中间图V
    denominator = max(1 - S, 1e-8);    % 数值稳定性：避免除以0
    V = (G - S) ./ denominator;
    % 归一化因子alpha
    max_V = max(V(:));
    if max_V < 1e-8
        alpha = 1;
    else
        alpha = 1 / max_V;
    end
    %5最终清晰图像U
    U = alpha * V;
    U = max(0, min(1, U)); % 最终值域截断
end