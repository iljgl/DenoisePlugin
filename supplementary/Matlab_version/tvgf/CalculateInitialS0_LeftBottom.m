function S0 = CalculateInitialS0_LeftBottom(G, WinRow, WinCol)
    % 输入：
    %   G - 去噪后的引导图像（单通道，double类型，[0,1]）
    %   WinRow - 窗口行数（如10）
    %   WinCol - 窗口列数（如10）
    % 输出：
    %   S0 - 散射光初始估计（与G同尺寸）
    [Gy, Gx] = size(G);
    S0 = zeros(Gy, Gx);
    % 预计算窗口半长（用于边界判断）
    halfWinRow = WinRow - 1;
    halfWinCol = WinCol - 1;
    for row = 1:Gy
        for col = 1:Gx
            % 定四分之一窗口，行检测
            if row - halfWinRow >=1
                % 不越界时
                r_start = row - halfWinRow;
                r_end = row;
            else
                % 越界反转
                r_start = row;
                r_end = row+halfWinRow;
                % 窗口大于图片时
                if r_start < 1
                    r_start = 1;
                    r_end = Gy;
                end
            end

            % 列检测
            if col + halfWinCol <= Gx
                % 不越界
                c_start = col;
                c_end = col + halfWinCol;
            else
                % 越界反转
                c_start = col - halfWinCol;
                c_end = col;
                % 窗口大于图片时
                if c_start < 1
                    c_start = 1;
                    c_end = Gx;
                end
            end

            %
            window_G = G(r_start:r_end, c_start:c_end);
            S0(row, col) = min(window_G(:));
        end
    end
end