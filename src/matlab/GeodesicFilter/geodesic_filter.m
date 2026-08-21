function [dist_o,window_o,prev_o,filtered] = geodesic_filter(img, point,W, alpha, sigma)
% GEODESIC_FILTER implements geodesic filtering (based on Dijkstra's shortest path algorithm)
%   Input: img - grayscale image (2D double, any size)
%         W - window size (odd number, default 11)
%         alpha - gray difference weight (default 1)
%         sigma - standard deviation of Gaussian kernel (default 10)
%   Output: filtered - filtered image

    if nargin < 5, sigma = 2; end
    if nargin < 4, alpha = 10; end
    if nargin < 3, W = 7; end
    if nargin < 2, point = floor(size(img)/2); end
    if mod(W,2)==0, error('The window size W must be an odd number'); end
    if sigma<=0, error('sigma>0'); end
    if alpha<=0, error('alpha>0'); end

    % Symmetric Padding
    half = floor(W/2);
    padded = padarray(img, [half half], 'symmetric');
    filtered = zeros(size(img));
    
    offsets = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1]; % 8 directions
    Eu_spatial_dists = sqrt(sum(offsets.^2, 2));%Pre-computed Euclidean distance
    
    %Flatten the window into a vector, with node numbers 1.. W^2
    Nnodes = W*W;
    prev_o = zeros(Nnodes, 1);
    window_o=zeros(Nnodes,1);
    dist_o=zeros(Nnodes,1);
    %Traverse each pixel (only processing the original image area) 
    for i = 1:size(img,1)
        for j = 1:size(img,2)
            %The starting index of the current window in filling the image
            r0 = i + half;
            c0 = j + half;
            %Extract window (W x W) 
            window = padded(r0-half:r0+half, c0-half:c0+half);
            
            %1.Calculate the geodesic distance from all pixels within the calculation window to the central pixel.  ----------

            %The index of the central pixel in the flattened image (MATLAB:Column-Major)
            center_idx = (half+1) + (half)*W;
            %Initialization distance
            dist = inf(Nnodes, 1);
            dist(center_idx) = 0;
            visited = false(Nnodes, 1);
            prev = zeros(Nnodes, 1);
            %Priority Queue (Simple implementation using an array and linear search)
            Q = center_idx;
            
            %Calculate the row and column coordinates (within the window) for each node in advance
            [rows, cols] = ndgrid(1:W, 1:W);
            rows = rows(:); cols = cols(:);%3x3:(1,1),(2,1),(3,1),(1,2),(2,2),(3,2),(1,3),(2,3),(3,3)
            
            while ~isempty(Q)
                %Select the node with the smallest distance
                [~, idx] = min(dist(Q));
                u = Q(idx);
                Q(idx) = [];
                if visited(u), continue; end
                visited(u) = true;
                %The row and column of the current node
                ur = rows(u); uc = cols(u);
                %Check 8 adjacent areas
                for k = 1:size(offsets,1)
                    vr = ur + offsets(k,1);
                    vc = uc + offsets(k,2);
                    %Boundary Check
                    if vr < 1 || vr > W || vc < 1 || vc > W
                        continue;
                    end
                    %Column-Major
                    v = (vc-1)*W + vr;
                    if visited(v)
                        continue;
                    end
                    %Calculate edge weights     
                    intensity_diff = window(ur,uc) - window(vr,vc);
                    w = sqrt(Eu_spatial_dists(k)^2 + alpha^2 * intensity_diff^2);
                    %update geodesic distance
                    new_dist = dist(u) + w;
                    if new_dist < dist(v)
                        dist(v) = new_dist;
                        prev(v) = u;
                        %If v is not in Q, then add it.
                        if ~any(Q == v)
                            Q = [Q; v];
                        end
                    end
                end
            end
            if i==point(1)&&j==point(2)
                prev_o=prev;
                window_o=window;
                dist_o=dist;end
            %"dist" represents the geodesic distance.
            %2.Calculate the Gaussian weights and perform weighted averaging  ----------
            weights = exp(-dist.^2 / (2*sigma^2));%The smaller the sigma value is, the smaller the weight corresponding to the large distance will be.
            %normalization
            weights = weights / sum(weights);
            %weighted average
            filtered(i,j) = sum(weights .* window(:));
        end
    end
end