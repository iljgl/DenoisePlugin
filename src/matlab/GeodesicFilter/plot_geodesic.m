function plot_geodesic(window, prev_vec, dist_vec, points)
% PLOT_GEODESIC 可视化窗口内测地线距离与最短路径
% 输入:
%   window   - W x W 二维窗口灰度图 (与geodesic_filter2输出的window_o一致)
%   prev_vec - W^2 x 1 列优先展平的前驱索引向量 (1-based, 中心点值为0)
%   dist_vec - W^2 x 1 列优先展平的测地线距离向量
%   points   - N x 2 待绘制的点坐标矩阵, 每行格式为 [x, y]
%              x: 列号, y: 行号, 与图像坐标系一致

    % ========== 参数校验与预处理 ==========
    W = size(window, 1);
    if size(window, 2) ~= W
        error('window必须是方阵');
    end
    if length(prev_vec) ~= W*W || length(dist_vec) ~= W*W
        error('prev_vec与dist_vec长度必须等于W^2');
    end
    
    % 向量重塑为方阵 (列优先, 与Matlab reshape行为完全一致)
    dist_mat = reshape(dist_vec, W, W);
    
    N_points = size(points, 1);
    % 生成差异化颜色集
    line_colors = lines(N_points);
    
    half = (W - 1) / 2;
    center_r = half + 1; % 中心点行号
    center_c = half + 1; % 中心点列号

    % ========== 图1: 测地线距离热力图 + 格点数值标注 ==========
    figure('Name', '测地线距离热力图');
    imagesc(dist_mat);
    
    % ===== 修改点：使用金色→紫色渐变 =====
    n_colors = 256;
    gold = [1, 0.84, 0];      % 金色 RGB（距离短）
    purple = [0.5, 0, 0.5];   % 紫色 RGB（距离长）
    custom_cmap = [linspace(gold(1), purple(1), n_colors)', ...
                   linspace(gold(2), purple(2), n_colors)', ...
                   linspace(gold(3), purple(3), n_colors)'];
    colormap(custom_cmap);
    % ===================================
    
    colorbar;
    title('The distance of the geodesic and weights inside the window)');
    hold on;
    
    % 每个像素格标注距离值
    for r = 1 : W
        for c = 1 : W
            text(c, r, sprintf('%.2f', dist_mat(r,c)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'white', ...
                'FontSize', 8, ...
                'FontWeight', 'bold');
        end
    end
    
    axis equal tight;
    set(gca, 'YDir', 'reverse'); % 对齐图像坐标系: 行1在顶部
    xlabel('Col (x)');
    ylabel('Row (y)');
    hold off;

    % ========== 图2: 窗口灰度图 + 多点测地线路径 ==========
    figure('Name', '测地线路径','Position',[100,100,900,900]);
    imagesc(window);
    colormap(gray);
    title('The geodesic path within the window (center point → each target point)');
    hold on;
    
    % 标记中心点
    plot(center_c, center_r, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(center_c + 0.3, center_r, 'center', ...
        'Color', 'red', 'FontWeight', 'bold', 'FontSize', 10);
    
    % 逐个点回溯并绘制路径
    for i = 1 : N_points
        x = points(i, 1); % 输入x = 列号
        y = points(i, 2); % 输入y = 行号
        
        % 边界校验
        if x < 1 || x > W || y < 1 || y > W
            warning('点(%d,%d)超出窗口范围, 已跳过', x, y);
            continue;
        end
        
        % ---- 回溯前驱链条: 从目标点追溯到中心点 ----
        path = []; % 每行存储 [行号, 列号]
        current_idx = (x - 1) * W + y; % 列优先1-based索引
        
        while current_idx ~= 0
            % 索引转行号、列号
            r = mod(current_idx - 1, W) + 1;
            c = floor((current_idx - 1) / W) + 1;
            path = [path; r, c];
            current_idx = prev_vec(current_idx);
        end
        
        % ---- 绘制折线路径 ----
        plot(path(:,2), path(:,1), '.-', ...
            'Color', line_colors(i,:), ...
            'LineWidth', 2, ...
            'MarkerSize', 12);
        
        % ---- 标注该点的测地线距离 ----
        point_dist = dist_vec( (x-1)*W + y );
        text(x , y+ 0.3, sprintf('d=%.2f', point_dist), ...
            'Color', line_colors(i,:), ...
            'FontSize', 10, ...
            'FontWeight', 'bold');
        
        % 标记目标点
        plot(x, y, 's', 'Color', line_colors(i,:), ...
            'MarkerSize', 8, 'LineWidth', 1.5);
    end
    
    hold off;
end