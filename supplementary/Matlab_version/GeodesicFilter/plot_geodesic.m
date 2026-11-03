function plot_geodesic(window, prev_vec, dist_vec, points)
% PLOT_GEODESIC  Visualizes geodesic distances and shortest paths within a local window.
% Inputs:
%   window   - W x W grayscale image patch (consistent with the window_o output of geodesic_filter2).
%   prev_vec - W^2 x 1 column vector of predecessor indices (column-major order, 1-based;
%              the value at the center pixel is 0, indicating the root of the path tree).
%   dist_vec - W^2 x 1 column vector of accumulated geodesic distances from the center pixel
%              (column-major order, consistent with MATLAB reshape convention).
%   points   - N x 2 matrix of query point coordinates; each row is [x, y], where x denotes
%              the column index and y denotes the row index, following the image coordinate system.

    % ===== Parameter validation and preprocessing =====
    W = size(window, 1);
    if size(window, 2) ~= W
        error('window must be a square matrix');
    end
    if length(prev_vec) ~= W*W || length(dist_vec) ~= W*W
        error('prev_vec and dist_vec must each have length W^2');
    end
    
    % Reshape the distance vector into a W x W matrix (column-major order,
    % consistent with MATLAB's default reshape behavior).
    dist_mat = reshape(dist_vec, W, W);
    
    N_points = size(points, 1);
    % Generate a set of perceptually distinct colors for the N_paths paths.
    line_colors = lines(N_points);
    
    half = (W - 1) / 2;
    center_r = half + 1; % Row index of the window center (root of the shortest-path tree).
    center_c = half + 1; % Column index of the window center.

    % ===== Figure 1: Geodesic distance heatmap with per-pixel numerical annotation =====
    figure('Name', 'Geodesic distance heat map');
    imagesc(dist_mat);
    
    % Custom colormap: gold-to-purple gradient encoding short-to-long geodesic distances.
    n_colors = 256;
    gold = [1, 0.84, 0];      % Gold RGB (short distance).
    purple = [0.5, 0, 0.5];   % Purple RGB (long distance).
    custom_cmap = [linspace(gold(1), purple(1), n_colors)', ...
                   linspace(gold(2), purple(2), n_colors)', ...
                   linspace(gold(3), purple(3), n_colors)'];
    colormap(custom_cmap);
    % =====================================================
    
    colorbar;
    title('The distance of the geodesic and weights inside the window)');
    hold on;
    
    % Annotate each grid cell with its scalar geodesic distance value.
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
    set(gca, 'YDir', 'reverse'); % Flip y-axis to match image coordinates (row 1 at top).
    xlabel('Col (x)');
    ylabel('Row (y)');
    hold off;

    % ===== Figure 2: Window intensity image with multi-point geodesic paths overlaid =====
    figure('Name', 'Geodesic path','Position',[100,100,900,900]);
    imagesc(window);
    colormap(gray);
    title('The geodesic path within the window (center point → each target point)');
    hold on;
    
    % Mark the center (root) pixel.
    plot(center_c, center_r, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(center_c + 0.3, center_r, 'center', ...
        'Color', 'red', 'FontWeight', 'bold', 'FontSize', 10);
    
    % Backtrack the predecessor chain and render the shortest path for each query point.
    for i = 1 : N_points
        x = points(i, 1); % x: column index.
        y = points(i, 2); % y: row index.
        
        % Boundary check: skip query points lying outside the window extent.
        if x < 1 || x > W || y < 1 || y > W
            warning('Point (%d,%d) is out of the window bounds; skipped.', x, y);
            continue;
        end
        
        % ---- Backtrack the predecessor chain from the target pixel to the center ----
        path = []; % Each row stores [row, col].
        current_idx = (x - 1) * W + y; % Column-major 1-based linear index.
        
        while current_idx ~= 0
            % Convert linear index back to row and column subscripts.
            r = mod(current_idx - 1, W) + 1;
            c = floor((current_idx - 1) / W) + 1;
            path = [path; r, c];
            current_idx = prev_vec(current_idx);
        end
        
        % ---- Render the polyline path ----
        plot(path(:,2), path(:,1), '.-', ...
            'Color', line_colors(i,:), ...
            'LineWidth', 2, ...
            'MarkerSize', 12);
        
        % ---- Annotate the accumulated geodesic distance at the target point ----
        point_dist = dist_vec( (x-1)*W + y );
        text(x , y+ 0.3, sprintf('d=%.2f', point_dist), ...
            'Color', line_colors(i,:), ...
            'FontSize', 10, ...
            'FontWeight', 'bold');
        
        % Mark the target (query) point with a square marker.
        plot(x, y, 's', 'Color', line_colors(i,:), ...
            'MarkerSize', 8, 'LineWidth', 1.5);
    end
    
    hold off;
end