function floatsave(tifimage, filename)
i = Tiff(filename, 'w');
tagstruct.ImageLength = size(tifimage, 1);
tagstruct.ImageWidth = size(tifimage, 2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
% tagstruct.Photometric = Tiff.Photometric.BlackIsZero;  % 灰度图像
tagstruct.BitsPerSample = 32;%  % single 对应 32 位，double 对应 64 位
% tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.SamplesPerPixel = 1;  % 单通道（灰度图）
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
% 关键：指定数据类型为浮点
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;  % IEEE 浮点格式

if isa(tifimage, 'double')
    tagstruct.BitsPerSample = 64;  % double 对应 64 位
end
i.setTag(tagstruct)
i.write(tifimage);
% i.writeDirectory();
i.close();
end
% function myimwritetifstack(currentStack,output_filename)
% 
% % keys = [8, 16];  
% % values = [8, 16];
% % Bits = dictionary(keys, values);
% [height, width, numFrames] = size(currentStack);
% bits = 8;
% % 设置TIFF标签结构
% tagstruct = struct();
% tagstruct.ImageLength = height;
% tagstruct.ImageWidth = width;
% tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
% tagstruct.BitsPerSample = bits; % 根据实际数据类型调整
% tagstruct.SamplesPerPixel = 1;
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% tagstruct.Compression = Tiff.Compression.None; % 无压缩
% tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
% % 创建Tiff对象
% t = Tiff(output_filename, 'w');
% for round = 1:numFrames
%     if round==1
%         t.setTag(tagstruct);
%         t.write(uint8(currentStack(:, :, round)));
%     else
%         t.writeDirectory(); % 新建一个目录（页）
%         t.setTag(tagstruct);
%         t.write(uint8(currentStack(:, :, round)));
%     end
% end
% % 关闭Tiff对象
% t.close();
% end