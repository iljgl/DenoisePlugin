clear
clc
addpath('./Utils/');
addpath('./IterativeDeblur/');
%% 
% filename = 'img_000000029_Channel2_000-380.tif';
filename = 'EzrA300_img_000000029_Channel2_000-1-t.tif';
tif_stack = imreadTiff(filename);

filename_prefix = filename(1:end-4);
lambda=0.561; %μm 0.488
NA = 1.49;
z=0;
n=64;
r=2;%Cut Airy disk radius at 2μm
gpu=0;
iteration=10;
pixel = 0.0433; % μm
psf=kernel(pixel, lambda, NA,z,n,r);
% floatsave(psf,fullfile('',['PSF_double Pixel_',num2str(pixel),'NA_1.49', ' Lamda_',num2str(lambda),'.tif']));
imwritestack(uint16(65535*psf),fullfile('',['PSF_uint16 Pixel_',num2str(pixel),'NA_1.49', ' Lamda_',num2str(lambda),'.tif']),16) 

[~,~,frames] = size(tif_stack);
rule=2;
metho = ["Possion","Gauss"];

for num = 1: frames
    f = tif_stack(:,:,num);
    f = double(f);
    record_max = max(f(:));
    data_de=Iterative_deblur(f,psf,iteration,rule,gpu);
    n_d = ( data_de - min( data_de(:) ) ) / ( max(data_de(:)) - min(data_de(:)) );
    tif_stack(:,:,num) = uint16(record_max*(n_d));
    progressbar(num/frames);

end
fullfilename1 = [filename_prefix,'-DeconvIter_', num2str(iteration),' Pixel_',num2str(pixel), ' NA_',num2str(NA), ' Lamda_',num2str(lambda),char(metho(rule)),' 2.tif'];
imwritestack(tif_stack,fullfilename1,16) ;%(stack, filename,bit)
