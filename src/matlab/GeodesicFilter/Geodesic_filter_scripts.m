
clear
clc
img=imread('EzrA300_img_000000029_Channel2_000-1-t.tif');
% img=imread('raw_htv_p-0030.tif');
f=img;%接受任何类型的输入并在处理前统一变为double且仅最大值归一化为1,最小值不要求归为0.
f = double(f);
record_max = max(f(:));
f=f./record_max;%记录最大值

%% 
W=7;
alpha=10;
sigma=1;
[dist,window,prev,filtered] = geodesic_filter(f,[24,57], W, alpha, sigma);

o=uint16(filtered*record_max);
imwrite(o,'EzrA300_img_000000029_Channel2_000-1-t_w7a10s1_revised2.tif')
% imwrite(uint16(window*record_max),'EzrA300_img_000000029_Channel2_000-1-t_w7a10s1_revised2window.tif')
imwrite(window,'EzrA300_img_000000029_Channel2_000-1-t_w7a10s1_revised2window.tif')
imwrite(reshape(prev,[7,7]),'EzrA300_img_000000029_Channel2_000-1-t_w7a10s1_revised2window_prev.tif')

%% 

% dist_o=reshape(dist,[W,W]);
points = [2,5; 7,7; 1,4;4,1;7,1];
weights = exp(-dist.^2 / (2*sigma^2));
weights = weights / sum(weights);
plot_geodesic(window, prev, dist, points);
plot_geodesic(window, prev, weights, points);
%% 